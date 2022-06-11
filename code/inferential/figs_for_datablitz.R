library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(mfutils)
library(lemon)

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
theme_set(theme_bw(base_size = 18))
theme_surface <- list(
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_blank(),
    legend.background = element_blank(), legend.text = element_text(size = 12), legend.direction = "horizontal",
    legend.key.height = unit(1 / 3, "cm"), legend.key.width = unit(1, "cm")
  )
)

roi_color <- "#f804c8"
not_roi_color <- "black"
plot_sizes <- list(
    scatter_width = 5,
    scatter_height = 5,
    brain_width = 10,
    brain_height = 5
)

## functions:

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


dat_summary %>%
    filter(Term %in% c("hiloc1_subj_norm", "hiloc2_subj_norm")) %>%
    select(spatial_model:Estimate) %>%
    pivot_wider(id_cols = c("spatial_model", "region"), names_from = "Term", values_from = "Estimate") %>%
    group_by(spatial_model) %>%
    summarize(cor(hiloc1_subj_norm, hiloc2_subj_norm))

## subject / trial variance ratio ----

## average across runs:
dat_ratio <- dat_summary %>%
    filter(Term %in% c("hiloc1_subj_norm", "hiloc2_subj_norm")) %>%
    group_by(spatial_model, region) %>%
    summarize_if(is.numeric, mean)

## scatterplot:

p_ratio_scatter <- 
    dat_ratio %>%
    pivot_wider(id_cols = c("region"), names_from = "spatial_model", values_from = c("Estimate", "CI_L", "CI_U")) %>%
    ggplot(aes(Estimate_uv, Estimate_mv)) +
    geom_abline() +
    geom_point(
        shape = 21, size = 3.5, color = "white", fill = not_roi_color
        ) +
    coord_capped_cart(left = "both") +
    scale_x_continuous(breaks = c(0, 0.01), limits = c(0, 0.05)) +
    scale_y_continuous(breaks = c(0, 0.01), limits = c(0, 0.16)) +
    theme(
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line()
        ) +
    labs(x = "Univariate: var(subject)/var(trial)", y = "Multivariate: var(subject)/var(trial)")

ggsave(
    here("out", "datablitz_scatter_ratio.pdf"),
    dev = "pdf", 
    width = plot_sizes$scatter_width,
    height = plot_sizes$scatter_height
    )

p_ratio_scatter_color <- 
    dat_ratio %>%
    pivot_wider(id_cols = c("region"), names_from = "spatial_model", values_from = c("Estimate", "CI_L", "CI_U")) %>%
    mutate(is_roi = grepl("DorsAttn|ContA|ContB", region)) %>%
    ggplot(aes(Estimate_uv, Estimate_mv)) +
    geom_abline() +
    geom_point(
        aes(fill = is_roi),
        shape = 21, size = 3.5, color = "white"
        ) +
    scale_fill_manual(values = c("TRUE" = roi_color, "FALSE" = not_roi_color)) +
    annotate(
        geom = "text", label = "DAN, FPN", x = 0.025, y = 0.13, color = roi_color, hjust = 0,
        size = 10
        ) +
    coord_capped_cart(left = "both") +
    scale_x_continuous(breaks = c(0, 0.01), limits = c(0, 0.05)) +
    scale_y_continuous(breaks = c(0, 0.01), limits = c(0, 0.16)) +
    theme(
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line()
        ) +
    labs(x = "Univariate: var(subject)/var(trial)", y = "Multivariate: var(subject)/var(trial)")

ggsave(
    here("out", "datablitz_scatter_ratio_color.pdf"), 
    dev = "pdf",
    width = plot_sizes$scatter_width,
    height = plot_sizes$scatter_height
    )

## brains:

p_ratio_mv_brain <-
    dat_ratio %>%
    filter(spatial_model == "mv") %>%
    ggplot() +
    geom_brain(aes(fill = Estimate), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        option = "magma", na.value = "grey",
        limits = c(0, 0.16)
        ) +
    theme_surface +
    labs(title = "Multivariate")

p_ratio_uv_brain <-
    dat_ratio %>%
    filter(spatial_model == "uv") %>%
    ggplot() +
    geom_brain(aes(fill = Estimate), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        option = "magma", na.value = "grey",
        limits = c(0, 0.16)
        ) +
    theme_surface +
    labs(title = "Univariate")

(p_ratio_uv_brain + p_ratio_mv_brain) + theme(plot.title = element_text(hjust = 0.5))

ggsave(
    here("out", "datablitz_brain_ratio.pdf"),
    dev = "pdf",
    width = plot_sizes$brain_width,
    height = plot_sizes$brain_height
    )


## test--retest correlation ----

## extract:

dat_cor <- dat_summary %>%
    filter(Term %in% c("cor_subj__hiloc1__hiloc2")) %>%
    group_by(spatial_model, region) %>%
    summarize_if(is.numeric, mean)

dat_cor %>%
    pivot_wider(id_cols = c("region"), names_from = "spatial_model", values_from = c("Estimate", "CI_L", "CI_U")) %>%
    mutate(is_roi = grepl("DorsAttn|ContA|ContB", region)) %>%

    ggplot(aes(Estimate_uv, Estimate_mv)) +
    geom_abline() +
    geom_point(
        aes(fill = is_roi),
        shape = 21, size = 3, color = "white"
        ) +
    annotate(
        geom = "text", label = "DAN, FPN", x = 0.3, y = -0.3, color = roi_color, hjust = 0,
        size = 10
    ) +
    scale_fill_manual(values = c("TRUE" = roi_color, "FALSE" = not_roi_color)) +

    coord_capped_cart(left = "both") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.55, 1)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.55, 1)) +
    theme(
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line()
        ) +

    theme(legend.position = "none", panel.grid = element_blank()) +
    labs(x = "Univariate: cor(test, retest)", y = "Multivariate: cor(test, retest)")

ggsave(
    here("out", "datablitz_scatter_trr.pdf"),
    dev = "pdf",
    width = plot_sizes$scatter_width,
    height = plot_sizes$scatter_height
    )

dat_cor %>%
    pivot_wider(id_cols = c("region"), names_from = "spatial_model", values_from = c("Estimate", "CI_L", "CI_U")) %>%
    mutate(is_roi = grepl("DorsAttn|ContA|ContB", region)) %>%

    ggplot(aes(Estimate_uv, Estimate_mv)) +
    geom_abline() +
    geom_point(
        aes(fill = is_roi, alpha = ifelse(Estimate_uv > 0.7 | Estimate_mv > 0.7, 1, 0.2)),
        shape = 21, size = 3, color = "white"
        ) +
    annotate(
        geom = "text", label = "DAN, FPN", x = 0.3, y = -0.3, color = roi_color, hjust = 0,
        size = 10
    ) +
    scale_fill_manual(values = c("TRUE" = "#f804c8", "FALSE" = not_roi_color)) +
    scale_alpha_identity() +
    coord_capped_cart(left = "both") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.55, 1)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.55, 1)) +
    theme(
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line()
        ) +

    theme(legend.position = "none", panel.grid = element_blank()) +
    labs(x = "Univariate: cor(test, retest)", y = "Multivariate: cor(test, retest)")

ggsave(
    here("out", "datablitz_scatter_trr_thresh.pdf"),
    dev = "pdf",
    width = plot_sizes$scatter_width,
    height = plot_sizes$scatter_height
    )

## brains:

dat_cor %>%
    filter(spatial_model == "mv") %>%
    mutate(
        #Estimate = ifelse(Estimate > 0.7, Estimate, NA)
        #Estimate = ifelse(CI_L > 0, Estimate, NA)
        ) %>%
    ggplot() +
    geom_brain(aes(fill = Estimate), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        option = "magma", na.value = "grey",
#        breaks = scales::extended_breaks(4)
        ) +
    theme_surface


p_trr_mv_brain <-
    dat_cor %>%
    filter(spatial_model == "mv") %>%
    ggplot() +
    geom_brain(aes(fill = Estimate), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        option = "magma", na.value = "grey",
        limits = c(-0.55, 1)
        ) +
    theme_surface +
    labs(title = "Multivariate")

p_trr_uv_brain <-
    dat_cor %>%
    filter(spatial_model == "uv") %>%
    ggplot() +
    geom_brain(aes(fill = Estimate), atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
        option = "magma", na.value = "grey",
        limits = c(-0.55, 1)
        ) +
    theme_surface +
    labs(title = "Univariate")

(p_trr_uv_brain + p_trr_mv_brain) + theme(plot.title = element_text(hjust = 0.5))

ggsave(
    here("out", "datablitz_brain_trr.pdf"),
    dev = "pdf",
    width = plot_sizes$brain_width,
    height = plot_sizes$brain_height
    )
