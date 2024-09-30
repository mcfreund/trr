library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(mfutils)

# ROIs
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"

# Atlas
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
}

# Theme
theme_set(theme_bw(base_size = 12))
theme_surface <- list(
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_text(size = 7),
    legend.background = element_blank(), legend.text = element_text(size = 7), legend.direction = "horizontal",
    legend.key.height = unit(1 / 4, "cm"), legend.key.width = unit(1 / 3, "cm")
  )
)

mv_brm_fname <- "multivariate_bayesian_model.csv"  # Effects extracted by pull_bayes_ef()
uv_brm_fname <- "univariate_bayesian_model.csv"

# Summarize sampling statistics
vec2sum <- function(dat, term_name = NA, group_name = NA, alpha = .05) {
    ci <- c(alpha / 2, 1 - alpha / 2)
    as_tibble(
        list(
            Term = term_name, Grouping = group_name,
            Estimate = mean(dat),
            `Est.Error` = sd(dat),
            tstat = mean(dat) / sd(dat),
            CI_L = quantile(dat, probs = ci[[1]]),
            CI_U = quantile(dat, probs = ci[[2]]))
            )
}


summarize_samples <- function(file_name) {

    fpath <- here("out", "spatial", file_name)
    # Read file, drop the failed models and uninteresting coefficients
    samples <- readRDS(fpath)
    #print(sum(is.na(samples)))
    #samples <- samples[!is.na(samples)]
    samples <- lapply(samples, function(x) x[, !grepl("^r_subj|^lp__", names(x))])

    samples_summary <- lapply(samples, function(x) {

        # Calculate some "effect sizes" by mutate()
        x <- as_tibble(x) %>%
            mutate(
                hiloc1_subj_norm = sd_subj__hiloc1^2 / sigma^2,
                hiloc2_subj_norm = sd_subj__hiloc2^2 / sigma^2
            )

        # Summarizing
        bind_rows(lapply(as.list(x), vec2sum), .id = "Term") %>%
            mutate(Grouping = ifelse(grepl("subj", Term), "subj", NA))

    })

    samples_summary

}


## read and summarize:

mv_mcmc <- "mv_bayes_MCMC_coefs_fda.rds"  # Coefficients from every MCMC sample
uv_mcmc <- "uv_bayes_MCMC_coefs.rds"


mv_summary <- bind_rows(summarize_samples(mv_mcmc), .id = "region")
uv_summary <- bind_rows(summarize_samples(uv_mcmc), .id = "region")
dat_summary <- bind_rows(list(mv = mv_summary, uv = uv_summary), .id = "spatial_model")

rois <- dat_summary %>%
  filter(spatial_model == "mv", Grouping == "subj", Term == "cor_subj__hiloc1__hiloc2", CI_L > 0.5) %>%
  pull(region)

noise <- fread(here("out", "spatial", "noisepca__stroop.csv"))


noise[, pc := 1:.N, by = c("roi", "session", "level", "subj", "task", "wave")]
noise_sum <-
  noise[,
  .(proj = mean(proj), eigval = mean(eigval)), by = c("roi", "subj", "pc", "wave", "session")
  ]  ## ave over condition, wave, session
noise_sum <- noise_sum[, .(pev = sum(proj) / sum(eigval), proj1 = proj[1]), by = c("roi", "subj", "wave", "session")]
noise_sum <- noise_sum[, .(pev = mean(pev), proj1 = mean(proj1)), by = c("roi", "wave", "session")]
noise_sum <- noise_sum[session == "baseline", .(pev = mean(pev), proj1 = mean(proj1)), by = "roi"]
noise_sum <- noise_sum %>% rename(region = roi)

#noise_sum %>%
trr <- dat_summary %>%
  filter(Grouping == "subj", Term %in% c("hiloc1_subj_norm", "hiloc2_subj_norm")) %>%
  group_by(spatial_model, region) %>%
  summarize(Estimate = mean(Estimate)) %>%
  pivot_wider(names_from = "spatial_model", values_from = c("Estimate"))

# trr <- dat_summary %>%
#   filter(Grouping == "subj", Term == "cor_subj__hiloc1__hiloc2") %>%
#   select(spatial_model, region, Estimate) %>%
#   pivot_wider(names_from = "spatial_model", values_from = c("Estimate"))

trr <- merge(trr, noise_sum, by = "region")


#rois <- dat_summary %>% filter(tstat > 2 & Term == "b_hilo_allhi_vs_lo" & spatial_model == "uv") %>% pull(region) %>% unique
#trr <- trr %>% filter(region %in% rois)
cor(trr$mv, trr$uv, method = "spearman")
cor(trr$uv, trr$pev, method = "spearman")
cor(trr$mv, trr$pev, method = "spearman")
cor(trr$mv / (trr$uv + trr$mv), trr$pev, method = "spearman")
cor(trr$uv, trr$proj1, method = "spearman")
cor(trr$mv, trr$proj1, method = "spearman")
cor(trr$mv / (trr$uv + trr$mv), trr$proj1, method = "spearman")

pairs(cbind(diff = trr$mv - trr$uv, ratiodiff = ((trr$mv) / (trr$mv + trr$uv)), trr[, -1]))


trr %>%
  filter((mv - uv) > 0.05)

trr %>%
  filter(proj1 < 0.8)


pairs(cbind(trr[, -1], diff = trr$mv - trr$uv))


dat_summary %>%
  filter(Grouping == "subj", Term == "cor_subj__hiloc1__hiloc2", spatial_model == "mv") %>%
  slice_max(CI_L, n = 5) %>%
  pull(region)

dat_summary %>%
  filter(Grouping == "subj", Term == "cor_subj__hiloc1__hiloc2", spatial_model == "uv") %>%
  slice_max(CI_L, n = 5) %>%
  pull(region)

rois_min <- trr %>% filter(mv > uv, mv > 0.5) %>% slice_min(mv - uv, n = 10) %>% pull(region)
rois_max <- trr %>% filter(mv > uv, mv > 0.5) %>% slice_max(mv - uv, n = 10) %>% pull(region)
noise_bas <- noise[session == "baseline" & roi %in% c(rois_min, rois_max)]
noise_bas[roi %in% rois_min, roi_type := "min_diff"]
noise_bas[roi %in% rois_max, roi_type := "max_diff"]

noise_bas[, .(proj = mean(proj*eigval), eigval = mean(eigval)), by = c("subj", "pc", "roi_type")] %>%
  filter(pc < 20) %>%
  ggplot(aes(pc, proj, color = roi_type)) +
  geom_line(aes(group = interaction(subj, roi_type)), alpha = 0.1) +
  stat_summary(fun = mean, geom = "line")

ggsave("")
