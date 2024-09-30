## Adapted from figs_ms_2023-05-31.rmd, removing some irrelavant parts


## ----setup---------


library(here)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
library(brms)
library(posterior)
library(colorspace)
library(knitr)
library(mfutils)
library(ggsegSchaefer)
library(purrr)
library(furrr)
library(patchwork)

source(here("code", "_funs.R"))
source(here("code", "_constants.R"))
source(here("code", "inferential", "_plotting.R"))

## other constants, paths, etc...

n_cores <- 30
plan(multicore, workers = n_cores)

responses <- c("rda", "uv")
models <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma", "summarystat")
models_hbm <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma")
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
rois <- mfutils::schaefer2018_17_400_fsaverage5$key[[roi_col]]

path_out <- file.path("/data/nil-external/ccp/chenr/trr", "out")
model_info <- expand.grid(model_nm = "no_lscov_symm", response = responses, roi_nm = rois, session = sessions)


## --------summarystat TRR--------

## read data for summarystat model, subset relevant rows/cols, get subject-level stats:

d_summarystat <-
  fread(file.path(path_out, "spatial", "projections__stroop__rda__n_resamples100__demean_run__cv_allsess.csv")) %>%
  rename(ridge = value.ridge, rda = value.rda, session = test_session)
cols_keep <- c("roi", "uv", "rda", "variable", "trial", "subj", "wave", "session")
d_summarystat <- na.omit(d_summarystat[, ..cols_keep])
s_summarystat_subj <- d_summarystat[,
  lapply(.SD, mean),
  by = c("variable", "subj", "wave", "roi", "session"),
  .SDcols = responses
  ]

## get summarystat model ests:
s_summarystat_subj_trr <- s_summarystat_subj %>%
  pivot_longer(cols = all_of(responses), names_to = "response") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(hilo = hi - lo) %>%
  pivot_wider(id_cols = c("subj", "roi", "response", "session"), names_from = wave, values_from = hilo) %>%
  rename(roi_nm = roi) %>%
  mutate(model_nm = "summarystat")

## get summarystat model ests:
r_summarystat <- s_summarystat_subj_trr %>%
  group_by(roi_nm, model_nm, response, session) %>%
  summarize(r = cor(wave1, wave2)) %>%
  mutate(network = get_network(roi_nm))





## ---- t-statistics and ROI ------

## get z-values:
smrz <- function(x) {
  res <- data.table(
    est_mean = colMeans(x),
    est_median = apply(x, 2, median),
    est_sd = apply(x, 2, sd)
  )
  res$param <- colnames(x)
  res
}

nms <- c("hilo_wave1", "hilo_wave2")

## `div` method: extract posterior samples, compute z-scores per sample, then summarize over samples
s_pop_z_div <- future_pmap_bind(
  CJ(model_nm = "no_lscov_symm", response = "uv", roi_nm = rois, session = "baseline"),
  read_summarize_hbm,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(x, hyp) {
    sd_subj <- VarCorr(x, summary = FALSE)$subj$sd[, nms]
    sd_trial <- exp(fixef(x, summary = FALSE)[, paste0("sigma_", nms)])
    m <- fixef(x, summary = FALSE)[, nms]
    l <- list(pop_mean = m, sd_subj = sd_subj, sd_trial = sd_trial, z_subj = m / sd_subj, z_trial = m / sd_trial)
    rbindlist(lapply(l, smrz), idcol = "term")
    }
  )

## `div` method: extract posterior samples, summarize over samples, then compute z-scores with summary stats
s_pop_z_sum <- future_pmap_bind(
  CJ(model_nm = "no_lscov_symm", response = "uv", roi_nm = rois, session = "baseline"),
  read_summarize_hbm,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(x, hyp) {
    m <- colMeans(fixef(x, summary = FALSE)[, nms])
    sd_subj <- exp(colMeans(log(VarCorr(x, summary = FALSE)$subj$sd[, nms])))
    sd_trial <- exp(colMeans(fixef(x, summary = FALSE)[, paste0("sigma_", nms)]))
    res <- data.table(rbind(m, sd_subj, sd_trial, m / sd_subj, m / sd_trial))
    res$term <- c("pop_mean", "sd_subj", "sd_trial", "z_subj", "z_trial")
    melt(res, value.name = "est_mean", id.vars = "term", variable.name = "param")
    }
  )
s_pop_z <- full_join(
  s_pop_z_div %>% mutate(term = paste0(term, "_div")),
  s_pop_z_sum %>% mutate(term = paste0(term, "_sum"))
)

## add t-stat
s_pop_z <- full_join(
  s_pop_z[term == "pop_mean_div"] %>%
    mutate(est_mean = est_mean / est_sd, term = "pop_tstat_div", est_median = NA, est_sd = NA),
  s_pop_z
)

## average over waves:
s_pop_z_mean <- s_pop_z[, .(est_mean = mean(est_mean)), by = c("roi_nm", "term")]


## t_value
s_pop_t_brain <- s_pop_z_mean[term == "pop_tstat_div"] %>%
  arrange(-est_mean) %>%
  mutate(is_roi = row_number() <= 40) %>%
  rename(region = "roi_nm", value = est_mean)





## ------------HBM TRR------------


## get HBM posteriors:
s_hbm_posterior <- future_pmap_bind(
  model_info,
  read_summarize_hbm,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(mdl) {
    mu <- ranef(mdl, summary = FALSE)$subj  ## for conditional TRR
    S <- VarCorr(mdl, summary = FALSE)$subj  ## for marginal TRR
    if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
      mu_stroop1 <- mu[, , "hilo_wave1"]
      mu_stroop2 <- mu[, , "hilo_wave2"]
      trr <- S$cor[, "hilo_wave1", "hilo_wave2"]
    } else {
      mu_stroop1 <- mu[, , "hi_wave1"] - mu[, , "lo_wave1"]
      mu_stroop2 <- mu[, , "hi_wave2"] - mu[, , "lo_wave2"]
      trr <- NA  ## marginal ranefs not yet implemented for this parameterization
    }
    trr_conditional <- rep(0, dim(mu)[1])
    for (ii in seq_along(trr_conditional)) trr_conditional[ii] <- cor(mu_stroop1[ii, ], mu_stroop2[ii, ])
    stopifnot(length(trr_conditional) == length(trr))
    data.table(trr = trr, trr_conditional = trr_conditional, resample = seq_along(trr))
  }
)

s_hbm_posterior[, network := get_network(roi_nm, TRUE)]

## summarize HBM posterior

s_hbm_posterior_sum <-
  s_hbm_posterior[,
  .(
    trr_map = max_aposteriori(trr),
    trr_mean = mean(trr),
    trr_05 = quantile(trr, 0.05),
    trr_sd = sd(trr)
  ),
  by = c("model_nm", "roi_nm", "response", "network", "session")
  ]
s_hbm_posterior_sum <- melt(
  s_hbm_posterior_sum,
  id.vars = c("model_nm", "roi_nm", "response", "network", "session")
  )



##------Variability ratio -----------

## get variability ratio:

s_hbm_subj_ratio <- future_pmap_bind(
  model_info,
  read_summarize_hbm,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(x) {
    ranefs <- ranef(x, summary = FALSE)$subj
    if ("sd_subj__hilo_wave1" %in% variables(x)) {
      stroop1 <- ranefs[, , "hilo_wave1"]
      stroop2 <- ranefs[, , "hilo_wave2"]
      if ("b_sigma_Intercept" %in% variables(x))  {
        hyp <- "exp(sigma_Intercept) = 0"
      } else {
        hyp <- "exp((sigma_mean_wave1 + sigma_mean_wave2) / 2) = 0"
      }
    } else {
      stroop1 <- ranefs[, , "hi_wave1"] - ranefs[, , "lo_wave1"]
      stroop2 <- ranefs[, , "hi_wave2"] - ranefs[, , "lo_wave2"]
      hyp <- "exp((sigma_hi_wave1 + sigma_hi_wave2 + sigma_lo_wave1 + sigma_lo_wave2) / 4) = 0"
    }
    sigma_stroop1 <- sqrt(Var(stroop1, 1))
    sigma_stroop2 <- sqrt(Var(stroop2, 1))
    sigma_stroop <- (sigma_stroop1 + sigma_stroop2)/2
    sigma_resid <- hypothesis(x, hyp)$samples$H1
    ratio <- sigma_resid / sigma_stroop

    data.table(
      sigma_resid = sigma_resid, 
      sigma_stroop = sigma_stroop, 
      ratio = sigma_resid / sigma_stroop,
      resample = seq_along(ratio)
      )

  }
)

s_hbm_subj_ratio_sum <- s_hbm_subj_ratio[,
  .(ratio = max_aposteriori(ratio)),
  by = c("model_nm", "roi_nm", "response", "session")
  ] %>%
  rename(region = roi_nm, value = ratio) %>% 
  mutate(variable = "ratio", network = get_network(region, TRUE))


##----Save data--------
##
## Save data

# Populational Stroop effect for UV
uv_popt <- s_pop_t_brain %>%
  select(!term) %>%
  rename(uv_pop_tstat_div = value)

# ICC
ICC <- r_summarystat %>%
  ungroup() %>%
  filter(session == "baseline") %>%
  mutate(response = recode(response, rda = "mv")) %>%
  rename(region = roi_nm, ICC = r) %>%
  select(!c(model_nm, network, session)) %>%
  pivot_wider(names_from = response, values_from = ICC) %>%
  rename(mv_ICC = mv, uv_ICC = uv)

# HBM TRR & Variability ratio
HBM <- s_hbm_posterior_sum %>%
  rename(region = roi_nm) %>%
  rbind(s_hbm_subj_ratio_sum) %>%
  filter(session == "baseline", model_nm == "no_lscov_symm") %>%
  mutate(response = recode(response, rda = "mv")) %>%
  select(!c(session, model_nm, network)) %>%
  pivot_wider(names_from = c(response, variable), values_from = value)

# Combine
plt_dat <- uv_popt %>%
  inner_join(ICC, by = "region") %>%
  inner_join(HBM, by = "region") %>%
  print()

# Save
fwrite(plt_dat, file.path(path_out, "inferential", "plt_dat.csv"))