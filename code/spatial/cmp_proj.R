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

val_terms <- c(ridge = "value.ridge", rda = "value.rda", uv = "uv")
auc_term <- "auc_ridge"
do_session <- "baseline"
old_subjs_only <- TRUE  # Whether to use previous 18 subjects only
fnames <- c(
  w_divnorm = paste0("projections__stroop__rda__n_resamples100",
    "__divnorm_run__divnorm_vertex__cv_allsess.csv"),
  wo_divnorm = "projections__stroop__rda__n_resamples100__cv_allsess.csv"
)
models <- names(fnames)
responses <- names(val_terms)
core32_roi <- rois[core32]

dat <- lapply(fnames, function(x) as_tibble(fread(here("out", "spatial", x))))
dat <- bind_rows(dat, .id = "model") %>%
  filter(test_session == .env$do_session) %>%
  rename(region = roi) %>%
  na.omit()

# Subjects
if (old_subjs_only) {
  other_dat <- dat %>% filter(!subj %in% .env$subjs_old_complete)
  dat <- dat %>% filter(subj %in% .env$subjs_old_complete)
}
n_subj <- length(unique(dat$subj))

# # AUC
# auc_summary <- dat %>%
#   group_by(subj, region, wave, model) %>%
#   summarize(auc = mean(.data[[auc_term]])) %>%
#   group_by(region, wave, model) %>%
#   summarize(m_auc = mean(auc), sem_auc = sd(auc) / sqrt(n())) %>%
#   ungroup()

# # Plotting AUC
# clim <- range(auc_summary$m_auc)
# figs <- NULL
# for (wv in c("wave1", "wave2")) {
#   for (mdl in models) {
#     tmp_dat <- auc_summary %>%
#       select(region, m_auc, wave, model) %>%
#       filter(wave == .env$wv, model == .env$mdl)
#     fig <- brain_plot(tmp_dat, stat_term = "m_auc", lim = clim,
#       fig_title = paste0(wv, ", ", mdl))
#     if (is.null(figs)) figs <- fig else figs <- figs + fig
#   }
# }
# figs + plot_layout(nrow = 2) +
#   plot_annotation(title = paste("Mean", auc_term, "across subjects"))
# ggsave(here("out", "spatial", "mean_AUC.png"))

# Roughly estimated TRR
trr_summary <- bind_rows(lapply(val_terms, function(x) {
  dat %>%
    group_by(subj, wave, region, model, variable) %>%
    summarize(response = mean(.data[[x]])) %>%
    pivot_wider(names_from = variable, values_from = response,
      id_cols = c(subj, wave, region, model)) %>%
    mutate(hilo = hi - lo) %>%
    pivot_wider(names_from = wave, values_from = hilo,
      id_cols = c(subj, region, model)) %>%
    group_by(region, model) %>%
    summarize(trr = cor(wave1, wave2)) %>%
    select(region, model, trr) %>%
    ungroup()
}), .id = "response")
trr_all_subjs <- bind_rows(lapply(val_terms, function(x) {
  dat %>%
    bind_rows(other_dat) %>%
    group_by(subj, wave, region, model, variable) %>%
    summarize(response = mean(.data[[x]])) %>%
    pivot_wider(names_from = variable, values_from = response,
      id_cols = c(subj, wave, region, model)) %>%
    mutate(hilo = hi - lo) %>%
    pivot_wider(names_from = wave, values_from = hilo,
      id_cols = c(subj, region, model)) %>%
    group_by(region, model) %>%
    summarize(trr = cor(wave1, wave2)) %>%
    select(region, model, trr) %>%
    ungroup()
}), .id = "response")

# TRR for each response type, more or fewer subjects, with or without normalization
trr_full <- bind_rows(list(subset18 = trr_summary, full = trr_all_subjs),
  .id = "subjects")
clim <- max(range(trr_full$trr))
figs <- NULL
for (res in responses) {
  tmp <- trr_full %>%
    filter(response == .env$res) %>%
    filter(region %in% core32_roi) %>%
    pivot_wider(id_cols = c(region, subjects), names_from = model, values_from = trr)
  fig <- ggplot(tmp) +
    aes(x = w_divnorm, y = wo_divnorm, color = subjects) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "with divisive normalization", y = "without divisive normalization",
      title = res) +
    xlim(-clim, clim) +
    ylim(-clim, clim)
  if (is.null(figs)) figs <- fig else figs <- figs + fig
}
figs <- figs + plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(title = "TRR estimated by different methods over core32 regions")
figs
ggsave(here("out", "spatial", "TRR_res_core32.png"), width = 15, height = 6)

# Difference between ridge and uv for both models
trr_res_diff <- get_diff_dat(trr_summary, "response", "uv", "trr") %>%
  filter(response == "ridge - uv") %>%
  filter(region %in% core32_roi) %>%
  pivot_wider(id_cols = region, names_from = model, values_from = trr)
trr_res_diff_all_subjs <- get_diff_dat(trr_all_subjs, "response", "uv", "trr") %>%
  filter(response == "ridge - uv") %>%
  filter(region %in% core32_roi) %>%
  pivot_wider(id_cols = region, names_from = model, values_from = trr)
all_trr <- bind_rows(list(subset18 = trr_res_diff, full = trr_res_diff_all_subjs),
  .id = "subjects")
clim <- max(range(all_trr[, models]))
fig <- ggplot(all_trr, aes(x = w_divnorm, y = wo_divnorm, color = subjects)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "with divisive normalization", y = "without divisive normalization",
    title = "Difference between TRRs calculated by ridge and uv for core32 region") +
  xlim(-clim, clim) +
  ylim(-clim, clim)
fig
ggsave(here("out", "spatial", "TRR_ridge_vs_uv_core32.png"), width = 8, height = 6)

# Difference between two models
trr_diff <- get_diff_dat(trr_summary, "model", models[[1]], "trr")
trr_full <- bind_rows(trr_summary, trr_diff)
all_models <- unique(trr_full$model)

# # Plotting TRRs
# clim <- range(trr_summary$trr)
# figs <- NULL
# for (mdl in models) {
#   for (res in responses) {
#     fig <- brain_plot(trr_summary %>% filter(model == .env$mdl, response == .env$res),
#       stat_term = "trr", lim = clim, fig_title = paste0(mdl, ", ", res))
#     if (is.null(figs)) figs <- fig else figs <- figs + fig
#   }
# }
# figs + plot_layout(nrow = length(models), guides = "collect") +
#   plot_annotation(title = "Estimated TRR for each model")
# ggsave(here("out", "spatial", "estimated_TRR.png"))

# Distribution plots

# clim <- range(trr_summary$trr)
# figs <- NULL
# for (res in responses) {
#   fig <- ggplot(trr_full %>% filter(response == .env$res)) +
#     aes(x = model, y = trr, fill = model) +
#     geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#     # geom_jitter(size = 0.02, height = 0, width = 0.1) +
#     ylim(clim[[1]], clim[[2]]) +
#     labs(x = res, y = NULL) +
#     theme(axis.text.x = element_blank())
#   if (res == responses[[1]]) figs <- fig + labs(y = "TRR") else figs <- figs + fig
# }
# figs + plot_layout(nrow = 1, guides = "collect") +
#   plot_annotation(title = "Distribution of TRR or difference between TRRs over all parcels")
# ggsave(here("out", "spatial", "res_dist_TRR.png"))

clim <- range(trr_summary$trr)
figs <- NULL
for (mdl in all_models) {
  fig <- ggplot(trr_full %>% filter(model == .env$mdl)) +
    aes(x = response, y = trr, fill = response) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    # geom_jitter(size = 0.02, height = 0, width = 0.1) +
    ylim(clim[[1]], clim[[2]]) +
    labs(x = mdl, y = NULL) +
    theme(axis.text.x = element_blank())
  if (mdl == models[[1]]) figs <- fig + labs(y = "TRR") else figs <- figs + fig
}
figs + plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(title = "Distribution of TRR or difference between TRRs over all parcels")
ggsave(here("out", "spatial", paste0("mdl_dist_TRR_", n_subj, "_subj.png")),
  width = 12, height = 9)

clim <- range(trr_summary$trr)
figs <- NULL
for (mdl in all_models) {
  fig <- ggplot(trr_full %>% filter(model == .env$mdl, region %in% core32_roi)) +
    aes(x = response, y = trr, fill = response) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    # geom_jitter(size = 0.02, height = 0, width = 0.1) +
    ylim(clim[[1]], clim[[2]]) +
    labs(x = mdl, y = NULL) +
    theme(axis.text.x = element_blank())
  if (mdl == models[[1]]) figs <- fig + labs(y = "TRR") else figs <- figs + fig
}
figs + plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(title = "Distribution of TRR or difference between TRRs over core32 parcels")
ggsave(here("out", "spatial", paste0("mdl_dist_TRR_core32_", n_subj, "_subj.png")),
  width = 12, height = 9)

# # Correlation
# corr_summary <- dat %>%
#   pivot_wider(names_from = model, values_from = .env$val_term,
#     id_cols = c(subj, wave, region, trial)) %>%
#   group_by(subj, wave, region) %>%
#   summarize(corr = cor(.env$models[[1]], .env$models[[2]])) %>%
#   group_by(wave, region) %>%
#   summarize(m_corr = mean(corr), sem_corr = sd(corr) / sqrt(n())) %>%
#   ungroup()

# # Plotting correlation
# clim <- range(corr_summary$m_corr)
# f1 <- brain_plot(filter(corr_summary, wave == "wave1"), stat_term = "m_corr", lim = clim,
#   fig_title = "wave1")
# f2 <- brain_plot(filter(corr_summary, wave == "wave2"), stat_term = "m_corr", lim = clim,
#   fig_title = "wave2")
# f1 + f2 + plot_annotation(
#   title = "Mean correlation between the scores from both model")
# ggsave(here("out", "spatial", "mean_corr.png"))
