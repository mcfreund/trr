library(here)
library(tidyr)
library(dplyr)
library(data.table)
library(tibble)
library(mfutils)
library(doParallel)
library(foreach)
library(lme4)
library(brms)
library(ggsegSchaefer)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))
# source(here("code", "inferential", "_plotting.R"))


########################### Constants ###########################

n_core_brm <- 4  # Number of cores for parallelization within brm()
n_roi_used <- -1  # Number of rois to look at (set as -1 to use all)
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
waves <- waves[do_waves]
n_cores <- 20
tasks <- "Stroop"
mle_output <- "univariate_linear_model.csv"
bayes_output <- "uv_bayes_MCMC_coefs.rds"

# Make sure we don't use more cores than available
stopifnot(n_core_brm <= n_cores)

# Atlas
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"

# ROIs
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
}


############################## Utilities #############################

# For extracting group-level effects from lme models
pull_fixef <- function(x, nms = c("term", "b", "se", "tstat")) {
  res <- coef(summary(x))
  res <- cbind(rownames(res), data.table::data.table(res))
  names(res) <- nms
  res
}

# For extracting population-level and group-level effects from brms models
pull_bayes_ef <- function(mdl) {

  # Convert effect matrix to tibble and add the grouping variable
  ef2tibble <- function(x, group = NA) {
    as_tibble(x, rownames = "Term") %>%
      add_column(Grouping = group, .before = "Term")
  }

  # Fixed effects
  res_fix <- ef2tibble(fixef(mdl))

  # Random effects and the residual
  vcov_mdl <- VarCorr(mdl)
  res_rnd <- bind_rows(lapply(names(vcov_mdl), function(nm) {
    ef2tibble(vcov_mdl[[nm]][["sd"]], group = nm)
  }))
  res_rnd[res_rnd$Grouping == "residual__",  "Term"] <- "Residual"
  res_rnd[res_rnd$Grouping == "residual__",  "Grouping"] <- NA

  # Data
  dat <- mdl$data[[mdl$formula$formula[[2]]]]
  res_dat <- as_tibble(list(Term = "Data_sd", Estimate = sd(dat)))

  # Combine them as the output
  bind_rows(list(res_fix, res_rnd, res_dat))

}


########################### Load data ###########################

# read trial-wise univariate means:
means <- read_results(
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) paste0("trial-means_", atlas_nm, "-", roi_col, "_resid-errts.csv"),
  read_fun = fread,
  n_cores = n_cores
)

# wrangle into wide-form data.table:
means <- lapply(means, function(x) as.data.table(t(x), keep.rownames = "trialnum"))
means <- rbindlist(means, idcol = "wave_task_session_subject")
names(means)[!names(means) %in% c("trialnum", "wave_task_session_subject")] <- rois
means <- separate(means, "wave_task_session_subject", into = c("wave", "task", "session", "subj"))
means$trialnum <- as.numeric(gsub("trial_", "", means$trialnum))

# read behavioral data:
behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_alltasks.csv"), na.strings = c("", "NA"))
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves & subj %in% subjs]

# bind
stopifnot(nrow(behav) == nrow(means))
d <- merge(means, behav, by = c("wave", "task", "session", "subj", "trialnum"), allow.cartesian = TRUE)


######################## Summary stat ##########################

d_long <- d %>% melt(id.vars = names(behav), variable.name = "region")  ## wide to long
d_sum_subj <- d_long[, .(value = mean(value)), by = .(subj, wave, session, region, hilo_all)]
d_sum_group <- d_sum_subj %>%
  pivot_wider(id_cols = c("subj", "wave", "session", "region"), names_from = "hilo_all", values_from = "value") %>%
  mutate(hivlo = hi - lo) %>%
  group_by(wave, session, region) %>%
  summarize(
    stat = t.test(hivlo)$statistic,
    p = t.test(hivlo)$p.value
    )  # Splitting the table by "hilo_all" into "hi" and "lo", then compute t across subjects

# correct p-values across all regions (but not waves and sessions):
d_sum_group <- d_sum_group %>%
  group_by(wave, session) %>%
  mutate(p_fdr = p.adjust(p, "fdr"))

# # Plot:
# d_temp <- d_sum_group %>%
#   filter(session == "baseline", wave == "wave1")
# tmp <- brain_plot(d_temp, stat_term = "stat",
#   fig_title = "t-statistics, summary stat method, stroop, baseline, wave1")
# print(tmp)


################ Prepare data for hierarchical linear modeling ############

# Select part of the ROIs
if (n_roi_used > 0) {
  rois <- rois[1:n_roi_used]
}

d <- as_tibble(d[d$session == "baseline",]) %>%
  mutate(hilo_all = factor(hilo_all, c("lo", "hi"), ordered = TRUE),
    wave = factor(wave, c("wave1", "wave2"), ordered = TRUE)) %>%
  mutate(hiloc = ifelse(hilo_all == "hi", 0.5, -0.5)) %>%
  mutate(hiloc1 = hiloc * (wave == "wave1"), hiloc2 = hiloc * (wave == "wave2"),
    w1_vs_mw1 = ifelse(wave == "wave1", 1, 0), w2_vs_mw2 = ifelse(wave == "wave2", 1, 0)) %>%
  select(!hiloc)

# Set contrasts as the difference between two levels (twice the difference between level 2 and the mean)
contrasts(d$wave) <- matrix(c(-0.5, 0.5), nrow = 2, dimnames = list(c("wave1", "wave2"), "w2_vs_w1"))
contrasts(d$hilo_all) <- matrix(c(-0.5, 0.5), nrow = 2, dimnames = list(c("lo", "hi"), "hi_vs_lo"))
d


############################### lme4 modeling ##############################

# # fit model to single ROI:
# fit_full <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave * hilo_all + (wave * hilo_all | subj), d)
# summary(fit_full)
# fit_reduced <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave + hilo_all + (wave + hilo_all | subj), d)
# fit_nowave <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ hilo_all + (hilo_all | subj), d)
# anova(fit_nowave, fit_reduced, fit_full)


# fit all models
formulas <- paste0("`", rois, "` ~ wave * hilo_all + (0 + w1_vs_mw1 + w2_vs_mw2 | subj)",
  " + (0 + hiloc1 + hiloc2 | subj)")
fits <- mclapply(formulas, function(x) lmer(as.formula(x), d), mc.cores = n_cores)
names(fits) <- rois
b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")

# Saving
write.csv(b, here("out", "spatial", mle_output))


############################# Fit Bayesian model #########################

# Fix naming problems with brms() formulas
input_for_bayes <- d %>%
  filter(if_all(starts_with("17Networks"), ~ !is.na(.x))) %>%
  setNames(gsub("17Networks", "Networks", names(.)))
rois_bayes <- gsub("17Networks", "Networks", rois)

# Fit a Bayesian model
bayes_model <- as.formula(paste0(rois_bayes[[1]],
  " ~ wave * hilo_all + (0 + w1_vs_mw1 + w2_vs_mw2 | subj)",
  " + (0 + hiloc1 + hiloc2 | subj)"))
get_prior(bayes_model, input_for_bayes)
fit_bayes <- brm(bayes_model, input_for_bayes, cores = n_core_brm)

# Fit all models
formulas_bayes <- paste0(rois_bayes,
  " ~ wave * hilo_all + (0 + w1_vs_mw1 + w2_vs_mw2 | subj)",
  " + (0 + hiloc1 + hiloc2 | subj)")
fits_bayes <- mclapply(formulas_bayes, function(x) tryCatch(
    brm(as.formula(x), input_for_bayes, cores = n_core_brm),
    error = function(e) {
      print()
      print(e)
      print()
      return(NA)
    }
  ),
  mc.cores = min(length(formulas_bayes), n_cores %/% n_core_brm))
names(fits_bayes) <- rois  # Note: need to get back the "17" now!

# # Save the summary of coefficients
# b_bayes <- bind_rows(lapply(fits_bayes, pull_bayes_ef), .id = "region")
# write.csv(b_bayes, here("out", "spatial", bayes_output))

# Save the coefficients from all MCMC draws
b_bayes <- lapply(fits_bayes, as.data.frame)
saveRDS(b_bayes, here("out", "spatial", bayes_output))