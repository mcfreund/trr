# Compare the value and AUC for "schafer_full" and "schafer_diag"
#
# Author: Ruiqi Chen
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
library(here)

source(here("code", "inferential", "_plotting.R"))

full_dat_fname <- "projections__stroop__schafer_full__n_resamples100.csv"
diag_dat_fname <- "projections__stroop__schafer_diag__n_resamples100.csv"

full_dat <- read_csv(here("out", "spatial", full_dat_fname))
diag_dat <- read_csv(here("out", "spatial", diag_dat_fname))

new_dat <- bind_rows(list(full = full_dat, diag = diag_dat), .id = "model") %>%
  pivot_wider(names_from = model, values_from = c(value, auc)) %>%
  rename(region = roi)

# The NAs
new_dat %>% filter(is.na(value_full) | is.na(value_diag)) %>% arrange(subj, wave, trial, region)
new_dat <- new_dat %>% filter(!is.na(value_full) & !is.na(value_diag))

# Correlation
corr_summary <- new_dat %>%
  group_by(subj, wave, task, region) %>%
  summarize(corr = cor(value_full, value_diag)) %>%
  group_by(wave, task, region) %>%
  summarize(m_corr = mean(corr), sem_corr = sd(corr) / sqrt(n())) %>%
  ungroup()

# Plotting correlation
clim <- c(0.1, 0.7)
f1 <- brain_plot(filter(corr_summary, wave == "wave1"), stat_term = "m_corr", lim = clim,
  fig_title = "wave1")
f2 <- brain_plot(filter(corr_summary, wave == "wave2"), stat_term = "m_corr", lim = clim,
  fig_title = "wave2")
f1 + f2 + plot_annotation(
  title = "Mean correlation between the scores from full & diagonal model")
ggsave(here("out", "spatial", "mean_corr.png"))

# AUC
auc_summary <- new_dat %>%
  group_by(subj, wave, task, region) %>%
  summarize(auc_full = mean(auc_full), auc_diag = mean(auc_diag)) %>%
  group_by(wave, task, region) %>%
  summarize(m_auc_full = mean(auc_full), sem_auc_full = sd(auc_full) / sqrt(n()),
    m_auc_diag = mean(auc_diag), sem_auc_diag = sd(auc_diag) / sqrt(n())) %>%
  ungroup()

# Plotting AUC
clim <- c(0.45, 0.65)
f1 <- brain_plot(filter(auc_summary, wave == "wave1"), stat_term = "m_auc_full", lim = clim,
  fig_title = "wave1, full covariance matrix")
f2 <- brain_plot(filter(auc_summary, wave == "wave1"), stat_term = "m_auc_diag", lim = clim,
  fig_title = "wave1, diagonal covariance matrix")
f3 <- brain_plot(filter(auc_summary, wave == "wave2"), stat_term = "m_auc_full", lim = clim,
  fig_title = "wave2, full covariance matrix")
f4 <- brain_plot(filter(auc_summary, wave == "wave2"), stat_term = "m_auc_diag", lim = clim,
  fig_title = "wave2, diagonal covariance matrix")
(f1 + f2) / (f3 + f4) + plot_annotation(title = "Mean AUC across subjects")
ggsave(here("out", "spatial", "mean_AUC.png"))
