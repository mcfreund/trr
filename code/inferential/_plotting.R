# Plotting function & Generate figures
#
# Author: Ruiqi Chen
#
# To-do: consider using `posterior()` package to summarize results
#
# 07/28/2022 update: now using `prep_mdl()` and `mdl2sum()` to prepare data.
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
library(brms)

# Get file naming function
source(here::here("code", "_funs.R"))

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

# **brain_plot(): Plot a column in a tibble onto the brain**
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

# **vec2sum() - Summarize sampling statistics**
#
# Inputs:
# - `dat`: a vector or a dataframe column
# - `term_name`, `group_name`: names for `Term` and `Grouping`
# - `alpha`: determining the percentile
# - `one-sided`: whether to use `alpha` or `alpha`/2 probability. Default is TRUE!
#
# Output: a tibble()
vec2sum <- function(dat, term_name = NA, group_name = NA, alpha = .05, one_sided = TRUE) {
  if (one_sided) ci <- c(alpha, 1 - alpha) else ci <- c(alpha / 2, 1 - alpha / 2)
  res <- tibble(
    Term = term_name, Estimate = mean(dat), `Est.Error` = sd(dat),
    `CI.Lower` = quantile(dat, probs = ci[[1]]),
    `CI.Upper` = quantile(dat, probs = ci[[2]]),
    `Q.Lower` = ci[[1]], `Q.Upper` = ci[[2]], Grouping = group_name
  )
}

#  **get_trr() - Calculate test-retest reliability**
#
# The reliability is calculated by the correlation of estimated hi-lo contrast
# for each subject in each posterior sample, then summarized across all samples.
get_trr <- function(mdl, alpha = .05) {
  mu <- ranef(mdl, summary = FALSE)$subj
  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
    mu_stroop1 <- mu[, , "hilo_wave1"]
    mu_stroop2 <- mu[, , "hilo_wave2"]
  } else {
    mu_stroop1 <- mu[, , "hi_wave1"] - mu[, , "lo_wave1"]
    mu_stroop2 <- mu[, , "hi_wave2"] - mu[, , "lo_wave2"]
  }
  trr <- rep(0, dim(mu)[1])
  for (ii in seq_along(trr)) trr[ii] <- cor(mu_stroop1[ii, ], mu_stroop2[ii, ])
  vec2sum(trr, term_name = "TRR", group_name = "subj", alpha = alpha)
}

# **get_norm_sd() - Calculate the relative contribution of subject-level variations**
#
# Result is the ratio between sd(hi-lo) over subjects and sd(response), estimated over
# posterior samples.
get_norm_sd <- function(mdl, alpha = .05) {
  data_nm <- toString(mdl$formula[[1]][[2]])
  denom <- sd(mdl$data[[data_nm]])  # "Scale" of all input data
  mu <- ranef(mdl, summary = FALSE)$subj
  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
    mu_stroop1 <- mu[, , "hilo_wave1"]
    mu_stroop2 <- mu[, , "hilo_wave2"]
  } else {
    mu_stroop1 <- mu[, , "hi_wave1"] - mu[, , "lo_wave1"] # 4000*27
    mu_stroop2 <- mu[, , "hi_wave2"] - mu[, , "lo_wave2"]
  }
  mu_stroop <- (mu_stroop1 + mu_stroop2) / 2
  sd_hilos_norm <- list(
    norm_sd_hilo_wave1 = apply(mu_stroop1, 1, sd) / denom,
    norm_sd_hilo_wave2 = apply(mu_stroop2, 1, sd) / denom,
    norm_sd_hilo = apply(mu_stroop, 1, sd) / denom
  )
  bind_rows(lapply(sd_hilos_norm, vec2sum, group_name = "subj",
    alpha = alpha), .id = "Term")
}

# **mdl2sum() - Summarize statistics of interest from a model**
#
# Input:
# - a `brmsfit()` model
# - `model_name`: "full", "no_lscov", or "no_lscov_symm"
# - `response_name`: "rda", "uv", or "ridge"
# - `session`: "baseline", "reactive", "proactive"
# - alpha level (by default 0.05)
#
# Output: a tibble with the following terms:
# - `model`, `response`, `session`
# - `Term`, `Grouping`, `Estimate`, `Est.Error`
# - `CI.Lower` and `CI.Upper`: lower and upper bound of the credible interval defined by `alpha`.
# - `Q.Lower` and `Q.upper`: indicating the percentile for `CI.Lower` and `CI.Upper`
mdl2sum <- function(mdl, roi_val = NA, roi_term = "region", alpha = .05) {

  # Loo and WAIC
  res <- bind_rows(as_tibble(mdl$criteria$loo$estimates, rownames = "Term"),
                  as_tibble(mdl$criteria$waic$estimates, rownames = "Term")) %>%
    rename(`Est.Error` = SE)

  # Bayes R2
  r2 <- as_tibble(rstantools::bayes_R2(mdl, probs = c(alpha, 1 - alpha))) %>%
    mutate(Term = "Bayes_R2", `Q.Lower` = alpha, `Q.Upper` = 1 - alpha)
  # Replace Qx.xx with easy-to-handle names
  ci_term_ind <- grep("^Q[0-9.]+$", names(r2))
  stopifnot(length(ci_term_ind) == 2)  # Throw an error if no unique match found
  names(r2)[ci_term_ind] <- c("CI.Lower", "CI.Upper")

  res <- bind_rows(res, r2)

  # Population-level high-low contrast
  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {  # "no_lscov_symm"
    hypos <- c(hilo_wave1 = "hilo_wave1 > 0", hilo_wave2 = "hilo_wave2 > 0")
  } else {
    hypos <- c(hilo_wave1 = "hi_wave1 - lo_wave1 > 0",
      hilo_wave2 = "hi_wave2 - lo_wave2 > 0")
  }
  hypo_res <- as_tibble(hypothesis(mdl, hypos, alpha = alpha)$hypothesis) %>%
    mutate(Term = names(hypos), `Q.Lower` = alpha, `Q.Upper` = 1 - alpha) %>%
    select(-c(`Evid.Ratio`, `Post.Prob`, Star, Hypothesis))
  res <- bind_rows(res, hypo_res)

  # Test-retest reliability
  res <- bind_rows(res, get_trr(mdl, alpha = alpha))

  # Subject-level (normalized) variation of high-low contrast
  res <- bind_rows(res, get_norm_sd(mdl, alpha = alpha))

  # ROI term
  res[[roi_term]] <- roi_val

  res
}

# **prep_mdl() - Extract statistics from the saved models**
#
# Inputs: `model_names`, `response_names`, `sessions`, three vectors (see estimate_reliability.R)
# Output: `res`, a list named by "modelname__responsename__session"
prep_mdl <- function(model_names = c("full"), response_names = c("rda"),
  sessions = c("baseline"), in_path = here::here("out", "inferential", atlas_nm)) {
  res <- mfutils::enlist(mfutils::combo_paste(model_names, "__", response_names,
    "__", sessions))
  for (model_name in model_names) {
    for (response_name in response_names) {
      for (session in sessions) {
        model_info <- get_model_info(model_name, response_name, session)
        rds_name <- paste0(model_info$model_prefix, ".rds")
        files <- list.files(in_path, pattern = rds_name, full.names = T,
          recursive = T)
        curr_res <- bind_rows(lapply(files, function(f) {
          mdl <- readRDS(f)
          mdl_region <- basename(dirname(f))
          mdl2sum(mdl, roi_val = mdl_region)
        }))
        res[[paste0(model_name, "__", response_name, "__", session)]] <- curr_res
      }
    }
  }
  res
}


# Main function when running this script directly
if (sys.nframe() == 0) {

  library(here)

  # Data
  mv_mcmc <- "mv_bayes_MCMC_coefs_fda.rds"  # Coefficients from every MCMC sample
  uv_mcmc <- "uv_bayes_MCMC_coefs.rds"
  mv_full <- "mv_bayes_MCMC_coefs_schafer_full.rds"  # For "schafer_full"
  mv_diag <- "mv_bayes_MCMC_coefs_schafer_diag.rds"  # For "schafer_diag"

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


  ############# Comparing "schafer_full" and "schafer_diag" ###############

  full_mdl <- prep_dat_rds(here("out", "spatial", mv_full))
  diag_mdl <- prep_dat_rds(here("out", "spatial", mv_diag))
  fda_mdl <- prep_dat_rds(here("out", "spatial", mv_mcmc))

  # Random effect of hilo relative to trial-level error (the residual)
  clim <- c(0, 0.5)
  f_uv_hilo1 <- brain_plot(filter(diag_mdl, Term == "hiloc1_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "variance only, wave1")
  f_mv_hilo1 <- brain_plot(filter(full_mdl, Term == "hiloc1_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "variance and covariance, wave1")
  f_fda_hilo1 <- brain_plot(filter(fda_mdl, Term == "hiloc1_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "fda() model, wave1")
  f_uv_hilo2 <- brain_plot(filter(diag_mdl, Term == "hiloc2_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "variance only, wave2")
  f_mv_hilo2 <- brain_plot(filter(full_mdl, Term == "hiloc2_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "variance and covariance, wave2")
  f_fda_hilo2 <- brain_plot(filter(fda_mdl, Term == "hiloc2_subj_norm"),
    stat_term = "Estimate", lim = clim, fig_title = "fda() model, wave2")
  (f_uv_hilo1 + f_mv_hilo1 + f_fda_hilo1) / (f_uv_hilo2 + f_mv_hilo2 + f_fda_hilo2) + plot_annotation(
    title = "Effect of (hi_lo|subj) relative to trial-level error in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_norm_mv_cmp.png"))

  # Test-retest reliability
  clim <- c(-0.6, 1)
  f_uv_trr <- brain_plot(filter(diag_mdl, Term == "cor_subj__hiloc1__hiloc2"),
    stat_term = "Estimate", lim = clim, fig_title = "variance only")
  f_mv_trr <- brain_plot(filter(full_mdl, Term == "cor_subj__hiloc1__hiloc2"),
    stat_term = "Estimate", lim = clim, fig_title = "variance and covariance")
  f_fda_trr <- brain_plot(filter(fda_mdl, Term == "cor_subj__hiloc1__hiloc2"),
    stat_term = "Estimate", lim = clim, fig_title = "fda() model")
  f_uv_trr + f_mv_trr + f_fda_trr + plot_annotation(
    title = "Test-retest correlation of hi-lo contrast in Stroop baseline")
  ggsave(here("out", "spatial", "hilo_trr_mv_cmp.png"))


  ################ Comparing "fda" and univarite methods ################

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

}