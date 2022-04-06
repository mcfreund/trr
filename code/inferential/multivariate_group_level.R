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
subjs <- subjs_wave12_complete
do_waves <- c(1, 2)
n_cores <- 20
tasks <- "Stroop"
sessions <- "baseline"
fname <- "projections__stroop__rda_lambda_100__n_resamples100.csv"

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
rois <- unique(mv_proj$roi)
mv_proj_wide <- mv_proj %>%
    pivot_wider(id_cols = c("trial", "subj", "task", "wave", "variable"),
        names_from = "roi", values_from = "value") %>%
    rename("hilo_all" = "variable")
mv_proj_wide


## fit model to single ROI:

fit_full <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave * hilo_all + (wave * hilo_all | subj),
  mv_proj_wide)
summary(fit_full)
fit_reduced <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave + hilo_all + (wave + hilo_all | subj),
  mv_proj_wide)
fit_nowave <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ hilo_all + (hilo_all | subj),
  mv_proj_wide)
anova(fit_nowave, fit_reduced, fit_full)

## fit many models (all 10 ROIs)
formulas <- paste0("`", rois, "` ~ wave + hilo_all + (wave + hilo_all | subj)")
fits <- mclapply(formulas, function(x) lmer(as.formula(x), mv_proj_wide), mc.cores = n_cores)
names(fits) <- rois
b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")
