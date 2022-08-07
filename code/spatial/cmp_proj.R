# Compare the output of different spatial models
#
# Author: Ruiqi Chen
#
# 08/05/2022 Update: now we can use this script to compare the output of two different
# spatial models (correlation between them, AUC, and roughly-estimated TRR).
# Currently we use it to compare the same model with or without divisive normalization.
#
# Note that wo_divnorm seems to include some NAs (probably due to all-NA trial?)
#
# This is a simple script to compare the predicted value and AUC of "schafer_full" and
# "schafer_diag" classifiers, i.e., whether the covariances matters. The script computes
# the correlation between the predicted value from both classifiers for every (subj *
# wave * region), and plot the grand mean over the subjects for each wave and region.
# Trials with NA for any classifiers will be discarded before computing correlation.
#
# Besides, the script computes the grand mean of AUC over subjects and plot it for each
# classifier, wave and region.
#
# The plotting relies on the brain_plot() function from ./code/inferential/_plotting.R

library(tidyverse)
library(patchwork)
library(data.table)
library(here)

source(here("code", "inferential", "_plotting.R"))
source(here("code", "_constants.R"))
core32_roi <- rois[core32]

val_term <- "value.ridge"
auc_term <- "auc_ridge"
do_session <- "baseline"
fnames <- c(
  w_divnorm = paste0("projections__stroop__rda__n_resamples100",
    "__divnorm_run__divnorm_vertex__cv_allsess.csv"),
  wo_divnorm = "projections__stroop__rda__n_resamples100__cv_allsess.csv"
)
models <- names(fnames)

dat <- lapply(fnames, function(x) as_tibble(fread(here("out", "spatial", x))))
dat <- bind_rows(dat, .id = "model") %>%
  filter(test_session == .env$do_session) %>%
  rename(region = roi) %>%
  na.omit()

# AUC
auc_summary <- dat %>%
  group_by(subj, region, wave, model) %>%
  summarize(auc = mean(.data[[auc_term]])) %>%
  group_by(region, wave, model) %>%
  summarize(m_auc = mean(auc), sem_auc = sd(auc) / sqrt(n())) %>%
  ungroup()

# Plotting AUC
clim <- range(auc_summary$m_auc)
figs <- NULL
for (wv in c("wave1", "wave2")) {
  for (mdl in models) {
    tmp_dat <- auc_summary %>%
      select(region, m_auc, wave, model) %>%
      filter(wave == .env$wv, model == .env$mdl)
    fig <- brain_plot(tmp_dat, stat_term = "m_auc", lim = clim,
      fig_title = paste0(wv, ", ", mdl))
    if (is.null(figs)) figs <- fig else figs <- figs + fig
  }
}
figs + plot_layout(nrow = 2) + plot_annotation(title = "Mean AUC across subjects")
ggsave(here("out", "spatial", "mean_AUC.png"))

# Roughly estimated TRR
trr_summary <- dat %>%
  group_by(subj, wave, region, model, variable) %>%
  summarize(response = mean(.data[[val_term]])) %>%
  pivot_wider(names_from = variable, values_from = response,
    id_cols = c(subj, wave, region, model)) %>%
  mutate(hilo = hi - lo) %>%
  pivot_wider(names_from = wave, values_from = hilo,
    id_cols = c(subj, region, model)) %>%
  group_by(region, model) %>%
  summarize(trr = cor(wave1, wave2)) %>%
  select(region, model, trr) %>%
  ungroup()

# Difference between two models
trr_diff <- get_diff_dat(trr_summary, "model", models[[1]], "trr")

# Plotting TRR and difference between TRRs
clim <- range(trr_summary$trr)
figs <- NULL
for (mdl in models) {
  fig <- brain_plot(trr_summary %>% filter(model == .env$mdl), stat_term = "trr",
    lim = clim, fig_title = mdl)
  if (is.null(figs)) figs <- fig else figs <- figs + fig
}
clim <- range(trr_diff$trr)
for (mdl in unique(trr_diff$model)) {
  figs <- figs + brain_plot(trr_diff %>% filter(model == .env$mdl), stat_term = "trr",
    lim = clim, fig_title = mdl)
}
figs + plot_annotation(title = "Estimated TRR for each model")
ggsave(here("out", "spatial", "estimated_TRR.png"))

# Distribution plot
clim <- range(trr_summary$trr)
figs <- ggplot(bind_rows(trr_summary, trr_diff)) +
  aes(x = model, y = trr, fill = model) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(size = 0.02, height = 0, width = 0.1) +
  ylim(clim[[1]], clim[[2]]) +
  labs(x = "All regions", y = "Estimation") +
  theme(axis.text.x = element_blank())
figs <- figs +
  ggplot(bind_rows(trr_summary, trr_diff) %>% filter(region %in% .env$core32_roi)) +
  aes(x = model, y = trr, fill = model) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(size = 0.02, height = 0, width = 0.1) +
  ylim(clim[[1]], clim[[2]]) +
  labs(x = "Core32 regions", y = NULL) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(title = "Distribution of estimated TRR") +
  theme(axis.text.x = element_blank())
figs
ggsave(here("out", "spatial", "dist_estimated_TRR.png"))

# Correlation
corr_summary <- dat %>%
  pivot_wider(names_from = model, values_from = .env$val_term,
    id_cols = c(subj, wave, region, trial)) %>%
  group_by(subj, wave, region) %>%
  summarize(corr = cor(.env$models[[1]], .env$models[[2]])) %>%
  group_by(wave, region) %>%
  summarize(m_corr = mean(corr), sem_corr = sd(corr) / sqrt(n())) %>%
  ungroup()

# Plotting correlation
clim <- range(corr_summary$m_corr)
f1 <- brain_plot(filter(corr_summary, wave == "wave1"), stat_term = "m_corr", lim = clim,
  fig_title = "wave1")
f2 <- brain_plot(filter(corr_summary, wave == "wave2"), stat_term = "m_corr", lim = clim,
  fig_title = "wave2")
f1 + f2 + plot_annotation(
  title = "Mean correlation between the scores from both model")
ggsave(here("out", "spatial", "mean_corr.png"))
