library(here)
library(readr)
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
n_roi_used <- -1  # Number of rois to look at (set as -1 to use all)
subjs <- subjs_wave12_complete
do_waves <- c(1, 2)
n_cores <- 20
tasks <- "Stroop"
sessions <- "baseline"
fname <- "projections__stroop__rda_lambda_100__n_resamples100.csv"
output_fname <- "multivariate_linear_model.csv"


## for extracting group-level effects from models:
pull_fixef <- function(x, nms = c("term", "b", "se", "tstat")) {
  res <- coef(summary(x))
  res <- cbind(rownames(res), data.table::data.table(res))
  names(res) <- nms
  res
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


## fit model to single ROI:

nowave_model <- as.formula(paste0("`", rois[[1]], "` ~ hilo_all + (hilo_all | subj)"))
reduced_model <- as.formula(paste0("`", rois[[1]], "` ~ wave + hilo_all + (wave + hilo_all | subj)"))
full_model <- as.formula(paste0("`", rois[[1]], "` ~ wave * hilo_all + (wave * hilo_all | subj)"))

fit_nowave <- lmer(nowave_model, mv_proj_wide)
fit_reduced <- lmer(reduced_model, mv_proj_wide)
fit_full <- lmer(full_model, mv_proj_wide)

summary(fit_full)
anova(fit_nowave, fit_reduced, fit_full)

## fit many models (all ROIs)
formulas <- paste0("`", rois, "` ~ wave + hilo_all + (wave + hilo_all | subj)")
fits <- mclapply(formulas, function(x) lmer(as.formula(x), mv_proj_wide), mc.cores = n_cores)
names(fits) <- rois
b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")

## Plotting and saving
tmp <- brain_plot(b, eff_term = "term", eff = "hilo_alllo", lim = c(-12, 3), direct = -1,
  fig_title = "t-statistics for hi-lo, multivariate method, stroop, baseline")
print(tmp)
write.csv(b, here("out", "spatial", output_fname))