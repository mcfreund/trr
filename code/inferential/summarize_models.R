# **Summarizing the statistics from hierarchical models**
#
# Author: Ruiqi Chen
# Version: 07/29/2022
#
# To-do:
# - Investigate the posterior distribution of TRR
# - Parallelize `prep_mdl()`

library(here)
library(tidyr)
library(dplyr)
library(ggsegSchaefer)
library(mfutils)
library(brms)
library(posterior)

# Get file naming function
source(here("code", "_funs.R"))

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

# **rvar_map() - Maximum A Posterior Estimation**
#
# Inputs:
# - x: an `rvar`
# - n: number of points for density calculation, by default 100
#
# Ouput: a numeric value, the peak position of the density distribution
rvar_map <- function(x, n = 100) {
  segs <- seq(min(x), max(x), length.out = n)
  d <- density(x, segs)
  segs[[which.max(d)]]
}

# **get_summary() - Summarize sampling statistics**
#
# Inputs:
# - `dat`: a vector; or a dataframe/matrix with N columns (variables)
# - `term_name`, `group_name`: names for `Term` and `Grouping`, can be a length-N vector
# - `alpha`: probability to calculate percentiles, by default 0.05
# - `one-sided`: whether to use `alpha` or `alpha`/2 probability. Default is TRUE!
# - `n_chains`: number of chains, by default 1.
#
# Output: a tibble() with N rows and following columns:
# - `Term`: name of the variable
# - `Estimate`, `Est.Error`: mean and sd
# - `MAP`: maximum a posterior estimation
# - `rhat`, `ess_bulk`, `ess_tail`: convergence diagnostic statistics
# - `CI.Lower`, `CI.Upper`, `Q.Lower`, `Q.Upper`: percentiles and corresponding probabilities
# - `Grouping`: name of the groups
get_summary <- function(dat, term_name = NA, group_name = NA, alpha = .05,
  one_sided = TRUE, n_chains = 1) {

  # Transform to draws data frame
  if (length(dim(dat)) < 3) {
    x <- as.data.frame(dat)
    if (!is.na(term_name)) names(x) <- term_name
    n_iter <- dim(x)[1] / n_chains
    if (dim(x)[1] %% n_chains) stop("Number of samples is not divisible by number of chains!")
    x[[".chain"]] <- rep(1:n_chains, each = n_iter)
  }
  x <- as_draws_rvars(x)

  # Credible Interval
  if (one_sided) ci <- c(alpha, 1 - alpha) else ci <- c(alpha / 2, 1 - alpha / 2)
  ci_res <- summary(x, "quantile2", .args = list(probs = ci)) %>%
    mutate(variable = NULL, `Q.Lower` = ci[[1]], `Q.Upper` = ci[[2]])
  names(ci_res)[1:2] <- c("CI.Lower", "CI.Upper")

  # MAP
  map_res <- unlist(lapply(x, rvar_map))

  res <- summary(x, Estimate = mean, `Est.Error` = sd, default_convergence_measures()) %>%
    bind_cols(ci_res) %>%
    mutate(MAP = map_res, Grouping = group_name) %>%
    relocate(MAP, .after = `Est.Error`) %>%
    rename(Term = variable)
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
  get_summary(trr, term_name = "TRR", group_name = "subj",
    alpha = alpha, n_chains = dim(mdl$fit)[2])
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
  bind_rows(lapply(sd_hilos_norm, get_summary, group_name = "subj",
    alpha = alpha, n_chains = dim(mdl$fit)[2]), .id = "Term")
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
  if ("b_hilo_wave1" %in% variables(mdl)) {
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
  sessions = c("baseline"), in_path = here("out", "inferential", atlas_nm)) {

  res <- mfutils::enlist(mfutils::combo_paste(model_names, "__", response_names,
    "__", sessions))
  for (model_name in model_names) {
    for (response_name in response_names) {
      for (session in sessions) {
        model_info <- get_model_info(model_name, response_name, session)
        rds_name <- paste0(model_info$model_prefix, ".rds")
        files <- list.files(in_path, pattern = rds_name, full.names = T,
          recursive = T)
        if (length(files) == 0) next
        curr_res <- bind_rows(lapply(files, function(f) {
          mdl <- readRDS(f)
          mdl_region <- basename(dirname(f))
          print(paste("Now working on", model_name, response_name,
            session, mdl_region, "..."))
          mdl2sum(mdl, roi_val = mdl_region)
        }))
        res[[paste0(model_name, "__", response_name, "__", session)]] <- curr_res
      }
    }
  }
  res

}


########################### Main function ##############################

res <- prep_mdl(model_names = c("no_lscov_symm"),
  response_names = c("uv", "rda"), in_path = here("out", "inferential", atlas_nm))
saveRDS(res, file = here("out", "inferential", atlas_nm, "Schaefer400_stats.rds"))