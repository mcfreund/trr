library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)

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

# Plotting function
brain_plot <- function(df, eff_term = NULL, eff = NULL, stat_term = "tstat",
                    lim = NULL, direct = 1, fig_title = "Example figure", savename = NULL) {
  if (!is.null(eff_term)) {
    df <- df %>% filter(.data[[eff_term]] %in% .env$eff)
    # df <- df %>% group_by(.data[[eff_term]])  # Will fail due to error in brain_join()
  }
  fig <- df %>%
    ggplot() +
    geom_brain(aes(fill = .data[[stat_term]]),
      atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
      limits = lim,
      direction = direct,
      option = "magma", na.value = "grey",
      breaks = scales::extended_breaks(4)
      ) +
    theme_surface +
    labs(title = fig_title, fill = NULL)

  # Saving or displaying
  if (is.null(savename)) {
    message()
    message("Note: you need to use 'print(brain_plot(...))'",
      " instead of 'brain_plot(...)' to display the result.")
    message()
  } else {
    ggsave(savename, plot = fig)
  }
  return(fig)
}


# Main function when running this script directly
if (sys.nframe() == 0) {

  library(here)
  library(readr)
  library(mfutils)

  # Data
  mv_brm_fname <- "multivariate_bayesian_model.csv"
  uv_brm_fname <- "univariate_bayesian_model.csv"

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


  ############################ lme4 Results #################################

  # fname1 <- "multivariate_linear_model.csv"
  # b1 <- read_csv(here("out", "spatial", fname1)) %>%
  #   mutate(term = ifelse(term == "hilo_alllo", "hilo_allhi", term)) %>%
  #   mutate(b = ifelse(term == "hilo_allhi", -b, b),
  #     tstat = ifelse(term == "hilo_allhi", -tstat, tstat))
  # f1 <- brain_plot(b1, eff_term = "term", eff = "hilo_allhi", lim = c(-6, 12),
  #   fig_title = "multivariate")

  # fname2 <- "univariate_linear_model.csv"
  # b2 <- read_csv(here("out", "spatial", fname2))
  # f2 <- brain_plot(b2, eff_term = "term", eff = "hilo_allhi", lim = c(-6, 12),
  #   fig_title = "univariate")

  # f3 <- brain_plot(b1, eff_term = "term", eff = "wavewave2", lim = c(-2, 2),
  #   fig_title = "multivariate")
  # f4 <- brain_plot(b2, eff_term = "term", eff = "wavewave2", lim = c(-2, 2),
  #   fig_title = "univariate")
  # f3 + f4 + plot_annotation(title = "t-statistics for wave effect in Stroop baseline")
  # ggsave(here("out", "spatial", "wave.png"))


  ############################ brms effects #################################

  prep_dat <- function(fpath) {

    tib <- read_csv(fpath) %>%
      select(!`...1`)

    # Change lo-hi main effect to hi-lo
    tib <- tib %>%
      mutate(nQ975 = -`Q2.5`, nQ25 = -`Q97.5`) %>%
      mutate(Estimate = if_else(Term == "hilo_alllo" & is.na(Grouping), -Estimate, Estimate)) %>%
      mutate(`Q2.5` = if_else(Term == "hilo_alllo" & is.na(Grouping), nQ25, `Q2.5`)) %>%
      mutate(`Q97.5` = if_else(Term == "hilo_alllo" & is.na(Grouping), nQ975, `Q97.5`)) %>%
      mutate(Term = ifelse(Term == "hilo_alllo", "hilo_allhi", Term)) %>%
      select(!c(nQ25, nQ975))

    # Add normalized effects: (hilo|subj) / res, (wave|subj) / res, res/data, abs(wave)/data, hilo/data
    norm_tib <- tib %>%
      unite(Term_Grouping, c(Term, Grouping), sep = "|", na.rm = TRUE) %>%
      pivot_wider(id_cols = c(region, Term_Grouping), names_from = Term_Grouping, values_from = Estimate) %>%
      mutate(`wavewave2_norm|subj` = `wavewave2|subj` / Residual,
        `hilo_allhi_norm|subj` = `hilo_allhi|subj` / Residual,
        Residual_norm = Residual / Data_sd,
        wavewave2_abs_norm = abs(wavewave2) / Data_sd,
        hilo_allhi_norm = hilo_allhi / Data_sd) %>%
      select(region, contains("norm")) %>%
      pivot_longer(!region, names_to = "Term_Grouping", values_to = "Estimate") %>%
      separate(Term_Grouping, c("Term", "Grouping"), sep = "\\|", fill = "right")

    tib <- bind_rows(tib, norm_tib) %>%
      arrange(region, Grouping, Term)

    return(tib)
  }
  uv_dat <- prep_dat(here("out", "spatial", uv_brm_fname))
  mv_dat <- prep_dat(here("out", "spatial", mv_brm_fname))

  clim <- c(0, 0.4)
  f_uv_hilo <- brain_plot(filter(uv_dat, Grouping == "subj"), eff_term = "Term", eff = "hilo_allhi_norm",
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f_mv_hilo <- brain_plot(filter(mv_dat, Grouping == "subj"), eff_term = "Term", eff = "hilo_allhi_norm",
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f_uv_hilo + f_mv_hilo + plot_annotation(
    title = "Effect of (hi_lo|subj) relative to trial-level error in Stroop baseline")
  ggsave(here("out", "spatial", "hilo.png"))

  clim <- c(0, 0.07)
  f_uv_trr <- brain_plot(filter(uv_dat, is.na(Grouping)), eff_term = "Term", eff = "wavewave2_abs_norm",
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f_mv_trr <- brain_plot(filter(mv_dat, is.na(Grouping)), eff_term = "Term", eff = "wavewave2_abs_norm",
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f_uv_trr + f_mv_trr + plot_annotation(
    title = "Magnitude of the fixed effect of wave / magnitude of input sd(y) in Stroop baseline")
  ggsave(here("out", "spatial", "wave.png"))
}