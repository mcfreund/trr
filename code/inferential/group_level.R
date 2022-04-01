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

atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
n_cores <- 20
tasks <- "Stroop"

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
    )

## correct p-values across all regions:
d_sum_group <- d_sum_group %>%
  group_by(wave, session) %>%
  mutate(p_fdr = p.adjust(p, "fdr"))


## plot:
d_sum_group %>%
  filter(session == "baseline", wave == "wave1") %>%
  ggplot() +
  geom_brain(aes(fill = stat), atlas = atlas, position = position_brain(side ~ hemi)) +
  scale_fill_viridis_c(
    option = "magma", na.value = "grey",
    breaks = scales::extended_breaks(4)
    ) +
  theme_surface +
  labs(title = "t-statistics, summary stat method, stroop, baseline, wave1", fill = NULL)



## hierarchical linear model (lme4)

d$hilo_all <- factor(d$hilo_all, levels = c("lo", "hi"))  ## ensure intercept is coded as 'lo'


## fit model to single ROI:

fit <- lmer(`17Networks_LH_VisCent_ExStr_1` ~ wave + hilo_all + (wave + hilo_all | subj), d[session == "baseline"])
summary(fit)
coef(summary(fit))


## fit many models (all ROIs)

formulas <- paste0("`", rois, "` ~ wave + hilo_all + (wave + hilo_all | subj)")
fits <- mclapply(formulas, function(x) lmer(as.formula(x), d[session == "baseline"]), mc.cores = n_cores)
names(fits) <- rois
b <- rbindlist(lapply(fits, pull_fixef), idcol = "region")
b