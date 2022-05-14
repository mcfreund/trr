# Plotting function & Generate figures
#
# Author: Ruiqi Chen
#
# When being sourced, this script provides a function `brain_plot()` that can plot
# some statistics for each parcel over the brain. Please refer to the comments above
# the definition of `brain_plot()` for its usage.
#
# When executed directly, this script plots several effects using `brain_plot()`.
# There are two functions to process the input: `prep_dat_csv()` and `prep_dat_rds()`,
# corresponding to the two ways to save results in "group_level.R" or "multi...level.R".
#
# The proprocessing function accepts a full path to the input file and returns a
# tibble. The codes below use this tibble and `brain_plot()` to make the plots.

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

# brain_plot(): Plot a column in a tibble onto the brain
#
# Inputs:
#  - `df`: a tibble, the names of the parcels should be in `df$region`
#  - `stat_term`: a string, name of the column in `df` to plot
#  - `eff_term` and `eff`: `eff_term` is a string, `eff` can be a string or a sequence of string,
#      if specified, `df` will be filtered by `df[[eff_term]] %in% eff`
#  - `lim`: a sequence of length 2, the limit of colormap to use (by default automatically determined)
#  - `direct`: can set as -1 to reverse the colormap
#  - `fig_title`: a string, title of the plot
#  - `savename`: NULL (don't save) or a string, full path to save the figure
#
# Output:
#    `fig`: A `ggplot2()` figure plotted with `geom_brain()` from package `ggseg`.
#
brain_plot <- function(df, stat_term = "tstat", eff_term = NULL, eff = NULL,
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

  # Data
  mv_brm_fname <- "multivariate_bayesian_model.csv"  # Effects extracted by pull_bayes_ef()
  uv_brm_fname <- "univariate_bayesian_model.csv"
  mv_mcmc <- "mv_bayes_MCMC_coefs.rds"  # Coefficients from every MCMC sample
  uv_mcmc <- "uv_bayes_MCMC_coefs.rds"


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

  # Summarize sampling statistics
  vec2sum <- function(dat, term_name = NA, group_name = NA, alpha = .05) {
    ci <- c(alpha / 2, 1 - alpha / 2)
    as_tibble(list(Term = term_name, Grouping = group_name,
      Estimate = mean(dat), `Est.Error` = sd(dat), tstat = mean(dat) / sd(dat),
      CI_L = quantile(dat, probs = ci[[1]]),
      CI_U = quantile(dat, probs = ci[[2]])))
  }

  # Preprocess the MCMC coefficients saved in .rds file
  prep_dat_rds <- function(fpath) {

    # Read file, drop the failed models and uninteresting coefficients
    samples <- readRDS(fpath)
    samples <- samples[!is.na(samples)]
    samples <- lapply(samples, function(x) x[, !grepl("^r_subj|^lp__", names(x))])

    # Summarize the estimations
    new_samples <- lapply(samples, function(x) {

      # Calculate some "effect sizes" by mutate()
      x <- as_tibble(x) %>%
        mutate(
          hiloc1_subj_norm = sd_subj__hiloc1 / sigma,
          hiloc2_subj_norm = sd_subj__hiloc2 / sigma
        )

      # Summarizing
      bind_rows(lapply(as.list(x), vec2sum), .id = "Term") %>%
        mutate(Grouping = ifelse(grepl("subj", Term), "subj", NA))

    })

    # Bind into a large tibble
    bind_rows(new_samples, .id = "region")
  }
  uv_mdl <- prep_dat_rds(here("out", "spatial", uv_mcmc))
  mv_mdl <- prep_dat_rds(here("out", "spatial", mv_mcmc))

  # Random effect of hilo relative to trial-level error (the residual)
  clim <- c(0, 0.5)
  f_uv_hilo1 <- brain_plot(filter(uv_mdl, Term == "hiloc1_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate, wave1")
  f_mv_hilo1 <- brain_plot(filter(mv_mdl, Term == "hiloc1_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate, wave1")
  f_uv_hilo2 <- brain_plot(filter(uv_mdl, Term == "hiloc2_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate, wave2")
  f_mv_hilo2 <- brain_plot(filter(mv_mdl, Term == "hiloc2_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate, wave2")
  (f_uv_hilo1 + f_mv_hilo1) / (f_uv_hilo2 + f_mv_hilo2) + plot_annotation(
    title = "Effect of (hi_lo|subj) relative to trial-level error in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_norm.png"))

  # Test-retest reliability
  clim <- c(-0.6, 1)
  f_uv_trr <- brain_plot(filter(uv_mdl, Term == "cor_subj__hiloc1__hiloc2"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f_mv_trr <- brain_plot(filter(mv_mdl, Term == "cor_subj__hiloc1__hiloc2"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f_uv_trr + f_mv_trr + plot_annotation(
    title = "Test-retest correlation of hi-lo contrast in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_trr.png"))

  # Preprocess the effects saved in a csv by pull_bayes_ef in the training code
  prep_dat_csv <- function(fpath) {

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
      mutate(
        `wavewave2_norm|subj` = `wavewave2|subj` / Residual,
        `hilo_allhi_norm|subj` = `hilo_allhi|subj` / Residual,
        Residual_norm = Residual / Data_sd,
        wavewave2_abs_norm = abs(wavewave2) / Data_sd,
        hilo_allhi_norm = abs(hilo_allhi) / Data_sd
      ) %>%
      select(region, contains("norm")) %>%
      pivot_longer(!region, names_to = "Term_Grouping", values_to = "Estimate") %>%
      separate(Term_Grouping, c("Term", "Grouping"), sep = "\\|", fill = "right")

    tib <- bind_rows(tib, norm_tib) %>%
      arrange(region, Grouping, Term)

    return(tib)
  }
  uv_dat <- prep_dat_csv(here("out", "spatial", uv_brm_fname))
  mv_dat <- prep_dat_csv(here("out", "spatial", mv_brm_fname))

  # Compute the t-statistics (by default all fixed effects including residual)
  b2t <- function(x, eff = NULL, grp = NA) {
    if (!is.null(eff)) x <- filter(x, Term %in% .env$eff)
    if (!is.null(grp)) x <- filter(x, Grouping %in% .env$grp)
    x %>%
      filter(!is.na(`Est.Error`)) %>%
      mutate(tstat = Estimate / `Est.Error`)
  }

  # Fixed effect of hilo (t-statistics)
  clim <- NULL  # Different scale, since mv's main effect should be generally >= 0
  f1 <- brain_plot(b2t(uv_dat, "hilo_allhi"), lim = clim, fig_title = "univariate")
  f2 <- brain_plot(b2t(mv_dat, "hilo_allhi"), lim = clim, fig_title = "multivariate")
  f1 + f2 + plot_annotation(
    title = "t statistics for the fixed effect of hi-lo in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_fixed_t.png"))

  # Fixed effect of wave (t-statistics)
  clim <- c(0, 2.6)
  f1 <- brain_plot(mutate(b2t(uv_dat, "wavewave2"), tstat = abs(tstat)), lim = clim, fig_title = "univariate")
  f2 <- brain_plot(mutate(b2t(mv_dat, "wavewave2"), tstat = abs(tstat)), lim = clim, fig_title = "multivariate")
  f1 + f2 + plot_annotation(
    title = "Magnitude of t statistics for the fixed effect of wave in Stroop baseline")
  ggsave(here("out", "spatial", "wave_fixed_t_abs.png"))

  # Trial-level error relative to the scale of uni- or multi-variate statistics
  clim <- c(0.97, 1.01)
  f1 <- brain_plot(filter(uv_dat, is.na(Grouping), Term == "Residual_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f2 <- brain_plot(filter(mv_dat, is.na(Grouping), Term == "Residual_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f1 + f2 + plot_annotation(
    title = "Trial-level error relative to the scale of the input sd(y) in Stroop baseline")
  ggsave(here("out", "spatial", "error_norm.png"))

  # Magnitude of the fixed effect of hi-lo contrast relative to the input scale
  clim <- c(0, 0.5)
  f1 <- brain_plot(filter(uv_dat, is.na(Grouping), Term == "hilo_allhi_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f2 <- brain_plot(filter(mv_dat, is.na(Grouping), Term == "hilo_allhi_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f1 + f2 + plot_annotation(
    title = "Fixed effect of hi-lo relative to the scale of the input sd(y) in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_fixed_norm.png"))

  # Magnitude of the fixed effect of wave relative to the input scale
  clim <- c(0, 0.11)
  f1 <- brain_plot(filter(uv_dat, is.na(Grouping), Term == "wavewave2_abs_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f2 <- brain_plot(filter(mv_dat, is.na(Grouping), Term == "wavewave2_abs_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f1 + f2 + plot_annotation(
    title = "Magnitude of the fixed effect of wave relative to the scale of the input sd(y) in Stroop baseline")
  ggsave(here("out", "spatial", "wave_fixed_norm.png"))

  # Random effect of hilo relative to trial-level error (the residual)
  clim <- c(0, 0.35)
  f_uv_hilo <- brain_plot(filter(uv_dat, Grouping == "subj", Term == "hilo_allhi_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f_mv_hilo <- brain_plot(filter(mv_dat, Grouping == "subj", Term == "hilo_allhi_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f_uv_hilo + f_mv_hilo + plot_annotation(
    title = "Effect of (hi_lo|subj) relative to trial-level error in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_norm.png"))

  # Random effect of wave relative to trial-level error (the residual)
  clim <- c(0, 0.11)
  f_uv_wave <- brain_plot(filter(uv_dat, Grouping == "subj", Term == "wavewave2_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "univariate")
  f_mv_wave <- brain_plot(filter(mv_dat, Grouping == "subj", Term == "wavewave2_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "multivariate")
  f_uv_wave + f_mv_wave + plot_annotation(
    title = "Effect of (wave|subj) relative to trial-level error in Stroop baseline")
  ggsave(here("out", "spatial", "wave_norm.png"))
}