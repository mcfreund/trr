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

rda_dat_fname <- "projections__stroop__rda_lambda_100__n_resamples100.csv"
full_dat_fname <- "projections__stroop__schafer_full__n_resamples100.csv"
diag_dat_fname <- "projections__stroop__schafer_diag__n_resamples100.csv"

full_dat <- read_csv(here("out", "spatial", full_dat_fname))
diag_dat <- read_csv(here("out", "spatial", diag_dat_fname))
rda_dat <- read_csv(here("out", "spatial", rda_dat_fname))

new_dat <- bind_rows(list(full = full_dat, diag = diag_dat, rda = rda_dat), .id = "model") %>%
  rename(region = roi)

# The NAs

new_dat %>% filter(is.na(value)) %>% arrange(subj, wave, trial, region)
new_dat %>% filter(is.na(value)) %>% 
  unique %>% 
  select(subj, wave, model) %>%
  table
ids <- new_dat %>% 
  filter(is.na(value)) %>% 
  mutate(id = paste0(subj, "_", wave, "_", region, "_", trial)) %>%
  pull(id) %>%
  unique
new_dat <- new_dat %>% 
  mutate(id = paste0(subj, "_", wave, "_", region, "_", trial)) %>%
  filter(!id %in% ids)


auc_summary <- new_dat %>%
  group_by(model, subj, wave, task, region) %>%
  summarize(auc = mean(auc)) %>%
  group_by(model, subj, task, region) %>%
  summarize(auc = mean(auc)) %>%
  group_by(model, task, region) %>%
  summarize(auc = mean(auc))

auc_summary %>%
  pivot_wider(names_from = "model", values_from = "auc") %>%
  ggplot(aes(diag, full)) +
  geom_abline() +
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) +
  geom_point()

auc_summary %>%
  pivot_wider(names_from = "model", values_from = "auc") %>%
  ggplot(aes(diag, rda)) +
  geom_abline() +
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) +
  geom_point()

auc_summary %>%
  pivot_wider(names_from = "model", values_from = "auc") %>%
  ggplot(aes(full, rda)) +
  geom_abline() +
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) +
  geom_point()

# Correlation

cor_summary <- new_dat %>%
  select(region, variable, trial, subj, task, wave, value, model) %>%
  pivot_wider(names_from = "model", values_from = "value") %>%
  group_by(subj, wave, task, region) %>%
  summarize(r_full_diag = cor(full, diag), r_full_rda = cor(full, rda), r_diag_rda = cor(diag, rda))

cor_summary %>%
  pivot_longer(cols = c("r_full_diag", "r_full_rda", "r_diag_rda")) %>%
  group_by(region, name) %>%
  summarize(r = tanh(mean(atanh(value)))) %>%
  
  ggplot(aes(name, r)) + 
  geom_boxplot()



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
  summarize(auc_full = mean(auc_full), auc_diag = mean(auc_diag), auc_rda = mean(auc_rda)) %>%
  group_by(wave, task, region) %>%
  summarize(m_auc_full = mean(auc_full), sem_auc_full = sd(auc_full) / sqrt(n()),
    m_auc_diag = mean(auc_diag), sem_auc_diag = sd(auc_diag) / sqrt(n())) %>%
  ungroup()

# Plotting AUC
clim <- c(0.45, 0.65)
f1 <- brain_plot(filter(auc_summary, wave == "wave1"), stat_term = "m_auc_diag", lim = clim,
  fig_title = "variance only, wave1")
f2 <- brain_plot(filter(auc_summary, wave == "wave1"), stat_term = "m_auc_full", lim = clim,
  fig_title = "variance and covariance, wave1")
f3 <- brain_plot(filter(auc_summary, wave == "wave2"), stat_term = "m_auc_diag", lim = clim,
  fig_title = "variance only, wave2")
f4 <- brain_plot(filter(auc_summary, wave == "wave2"), stat_term = "m_auc_full", lim = clim,
  fig_title = "variance and covariance, wave2")
(f1 + f2) / (f3 + f4) + plot_annotation(title = "Mean AUC across subjects")
ggsave(here("out", "spatial", "mean_AUC.png"))
