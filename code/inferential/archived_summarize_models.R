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
library(here)

source(here("code", "_constants.R"))

# ROIs
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
roi_idx <- core32

# Atlas
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
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

# Preprocess the MCMC coefficients saved in .rds file
prep_dat_rds <- function(fpath) {

  # Read file, keep ROIs and drop failed models
  samples <- readRDS(fpath)
  samples <- samples[roi_idx]
  samples <- samples[!is.na(samples)]

  # Calculate TRR from samples and drop uninteresting terms
  hilo1names <- grep("^r_subj\\[.*hiloc1\\]$", names(samples[[1]]), value = TRUE)
  hilo2names <- grep("^r_subj\\[.*hiloc2\\]$", names(samples[[1]]), value = TRUE)
  samples <- lapply(samples, function(x) {
    x$TRR <- NA
    hilo1 <- t(x[, hilo1names])
    hilo2 <- t(x[, hilo2names])
    for (i in seq_len(dim(x)[1])) x[i, "TRR"] <- cor(hilo1[, i], hilo2[, i])
    x[, !grepl("^r_subj|^lp__", names(x))]
  })

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

# Data
inputs <- list(
  fix_std__fda__baseline = "mv_bayes_MCMC_coefs_fda.rds",
  fix_std__uv__baseline = "uv_bayes_MCMC_coefs.rds",
  fix_std__lda_full__baseline = "mv_bayes_MCMC_coefs_schafer_full.rds",
  fix_std__lda_diag__baseline = "mv_bayes_MCMC_coefs_schafer_diag.rds"
) # input is the coefficients from every MCMC sample (`as.data.frame(fitted_brms_model)`)

res <- lapply(inputs, function(f) {
  prep_dat_rds(here("out", "spatial", "archive", f))
})
saveRDS(res, here("out", "spatial", "archive", "fix_std_mdl_stats.rds"))