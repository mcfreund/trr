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

atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
n_roi_used <- -1  # Number of rois to look at (set as -1 to use all)
subjs <- subjs_wave12_complete
do_waves <- c(1, 2)
n_cores <- 20
tasks <- "Stroop"
sessions <- "baseline"
fname <- "projections__stroop__rda_lambda_100__n_resamples100.csv"
output_fname <- "multivariate_linear_model.csv"

theme_set(theme_bw(base_size = 12))
theme_surface <- list(
    theme(
        axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_text(size = 7),
        legend.background = element_blank(), legend.text = element_text(size = 7), legend.direction = "horizontal",
        legend.key.height = unit(1 / 4, "cm"), legend.key.width = unit(1 / 3, "cm")
    )
)

## for extracting group-level effects from models:
pull_fixef <- function(x, nms = c("term", "b", "se", "tstat")) {
    res <- coef(summary(x))
    res <- cbind(rownames(res), data.table::data.table(res))
    names(res) <- nms
    res
}

## load data ----

waves <- waves[do_waves]
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
    rois <- get(atlas_nm)$key[[roi_col]]
    atlas <- schaefer17_400
    atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
    stop("not configured for atlas")
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
write.csv(b, here("out", "spatial", output_fname))

## Plotting
b %>%
    filter(term == "hilo_alllo") %>%
    ggplot() +
    geom_brain(aes(fill = tstat), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        limits = c(-12, 3),
        direction = -1,
        option = "magma", na.value = "grey",
        breaks = scales::extended_breaks(4)
        ) +
    theme_surface +
    labs(title = "t-statistics for hi-lo, multivariate method, stroop, baseline", fill = NULL)
