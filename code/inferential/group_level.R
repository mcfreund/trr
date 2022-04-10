library(here)
library(tidyr)
library(dplyr)
library(data.table)
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
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
waves <- waves[do_waves]
n_cores <- 20
tasks <- "Stroop"
output_fname <- "univariate_linear_model.csv"


## for extracting group-level effects from models:
pull_fixef <- function(x, nms = c("term", "b", "se", "tstat")) {
  res <- coef(summary(x))
  res <- cbind(rownames(res), data.table::data.table(res))
  names(res) <- nms
  res
}


## read trial-wise univariate means:
means <- read_results(
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) paste0("trial-means_", atlas_nm, "-", roi_col, "_resid-errts.csv"),
  read_fun = fread,
  n_cores = n_cores
)

## wrangle into wide-form data.table:
means <- lapply(means, function(x) as.data.table(t(x), keep.rownames = "trialnum"))
means <- rbindlist(means, idcol = "wave_task_session_subject")
names(means)[!names(means) %in% c("trialnum", "wave_task_session_subject")] <- rois
means <- separate(means, "wave_task_session_subject", into = c("wave", "task", "session", "subj"))
means$trialnum <- as.numeric(gsub("trial_", "", means$trialnum))


## read trial data:
behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_alltasks.csv"), na.strings = c("", "NA"))
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves & subj %in% subjs]

## bind
stopifnot(nrow(behav) == nrow(means))
d <- merge(means, behav, by = c("wave", "task", "session", "subj", "trialnum"), allow.cartesian = TRUE)

d_long <- d %>% melt(id.vars = names(behav), variable.name = "region")  ## wide to long


## fit models ----


## summary stat

## compute:
d_sum_subj <- d_long[, .(value = mean(value)), by = .(subj, wave, session, region, hilo_all)]
d_sum_group <- d_sum_subj %>%
  pivot_wider(id_cols = c("subj", "wave", "session", "region"), names_from = "hilo_all", values_from = "value") %>%
  mutate(hivlo = hi - lo) %>%
  group_by(wave, session, region) %>%
  summarize(
    stat = t.test(hivlo)$statistic,
    p = t.test(hivlo)$p.value
    )  # Splitting the table by "hilo_all" into "hi" and "lo", then compute t across subjects

## correct p-values across all regions (but not waves and sessions):
d_sum_group <- d_sum_group %>%
  group_by(wave, session) %>%
  mutate(p_fdr = p.adjust(p, "fdr"))


## Plot:
d_temp <- d_sum_group %>%
  filter(session == "baseline", wave == "wave1")
tmp <- brain_plot(d_temp, stat_term = "stat",
  fig_title = "t-statistics, summary stat method, stroop, baseline, wave1")
print(tmp)


## hierarchical linear model (lme4)
d$hilo_all <- factor(d$hilo_all, levels = c("lo", "hi"))  ## ensure intercept is coded as 'lo'


## fit model to single ROI:

fit_full <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave * hilo_all + (wave * hilo_all | subj),
  d[session == "baseline"])
summary(fit_full)
fit_reduced <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave + hilo_all + (wave + hilo_all | subj),
  d[session == "baseline"])
fit_nowave <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ hilo_all + (hilo_all | subj),
  d[session == "baseline"])
anova(fit_nowave, fit_reduced, fit_full)


## fit many models (all 400 ROIs)

formulas <- paste0("`", rois, "` ~ wave + hilo_all + (wave + hilo_all | subj)")
fits <- mclapply(formulas, function(x) lmer(as.formula(x), d[session == "baseline"]), mc.cores = n_cores)
names(fits) <- rois
b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")
b
write.csv(b, here("out", "spatial", output_fname))