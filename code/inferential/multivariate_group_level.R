library(here)
library(readr)
library(tidyr)
library(dplyr)
library(data.table)
library(tibble)
library(brms)
library(mfutils)
library(doParallel)
library(foreach)
library(lme4)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))
source(here("code", "inferential", "_plotting.R"))

# Constants
n_roi_used <- -1  # Number of rois to look at (set as -1 to use all)
subjs <- subjs_wave12_complete
do_waves <- c(1, 2)
n_cores <- 20
tasks <- "Stroop"
sessions <- "baseline"
fname <- "projections__stroop__rda_lambda_100__n_resamples100.csv"
mle_output <- "multivariate_linear_model.csv"
bayes_output <- "multivariate_bayesian_model.csv"

## for extracting group-level effects from lme models:
pull_fixef <- function(x, nms = c("term", "b", "se", "tstat")) {
  res <- coef(summary(x))
  res <- cbind(rownames(res), data.table::data.table(res))
  names(res) <- nms
  res
}

## For extracting random effects and their explained variances from brms models
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

## Read multivariate projection
mv_proj <- read_csv(here("out", "spatial", fname))
if (n_roi_used > 0) {
  rois <- unique(mv_proj$roi)[1:n_roi_used]
}
mv_proj_wide <- mv_proj %>%
  filter(roi %in% .env$rois) %>%
  pivot_wider(id_cols = c(trial, subj, task, wave, variable),
    names_from = roi, values_from = value) %>%
  rename(hilo_all = variable)
mv_proj_wide

# Fix naming problems with brms() formulas
input_for_bayes <- mv_proj_wide %>%
  setNames(gsub("17Networks", "Networks", names(.)))
rois_bayes <- gsub("17Networks", "Networks", rois)


## fit model to single ROI:

# nowave_model <- as.formula(paste0("`", rois[[1]], "` ~ hilo_all + (hilo_all | subj)"))
# reduced_model <- as.formula(paste0("`", rois[[1]], "` ~ wave + hilo_all + (wave + hilo_all | subj)"))
# full_model <- as.formula(paste0("`", rois[[1]], "` ~ wave * hilo_all + (wave * hilo_all | subj)"))

# fit_nowave <- lmer(nowave_model, mv_proj_wide)
# fit_reduced <- lmer(reduced_model, mv_proj_wide)
# fit_full <- lmer(full_model, mv_proj_wide)

# summary(fit_full)
# anova(fit_nowave, fit_reduced, fit_full)


# Fit a Bayesian model
bayes_model <- as.formula(paste0("`", rois_bayes[[1]], "` ~ wave + hilo_all + (wave + hilo_all | subj)"))
get_prior(bayes_model, input_for_bayes)
fit_bayes <- brm(bayes_model, input_for_bayes, cores = n_cores)

## fit many models (all ROIs)

# formulas <- paste0("`", rois, "` ~ wave + hilo_all + (wave + hilo_all | subj)")
# fits <- mclapply(formulas, function(x) lmer(as.formula(x), mv_proj_wide), mc.cores = n_cores)
# names(fits) <- rois
# b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")

formulas_bayes <- paste0(rois_bayes, " ~ wave + hilo_all + (wave + hilo_all | subj)")
fits_bayes <- mclapply(formulas_bayes, function(x) brm(as.formula(x), input_for_bayes),
  mc.cores = n_cores)
names(fits_bayes) <- rois  # Note: need to get back the "17" now!
b_bayes <- bind_rows(lapply(fits_bayes, pull_bayes_ef), .id = "region")

## Plotting and saving
# tmp <- brain_plot(b, eff_term = "term", eff = "hilo_alllo", lim = c(-12, 3), direct = -1,
#   fig_title = "t-statistics for hi-lo, multivariate method, stroop, baseline")
# print(tmp)
# write.csv(b, here("out", "spatial", mle_output))

write.csv(b_bayes, here("out", "spatial", bayes_output))