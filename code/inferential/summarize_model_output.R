library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(data.table)
library(brms)
library(posterior)
library(mfutils)

source(here("code", "inferential", "_parameters_reliability.R"))
source(here("code", "inferential", "_utils_reliability.R"))
source(here("code", "_atlas.R"))
source(here("code", "_constants.R"))
source(here("code", "_paths.R"))


## misc parameters ----

n_cores <- 8
plan(multicore, workers = n_cores)

path_summary <- here("out", "inferential", atlas_nm)

model_info <- rbind(
  ## all models for core32 parcels:
  CJ(model_nm = models_hbm, response = responses, roi_nm = core32_nms, session = "baseline"),
  ## only winning model (no_lscov_symm) for rest of parcels:
  CJ(model_nm = "no_lscov_symm", response = responses, roi_nm = rois[!rois %in% core32_nms], session = "baseline")
)

## flags for saving different files:

do_summary_stats <- FALSE
do_hbm_posteriors <- TRUE
do_hbm_criteria <- FALSE

## posterior summary functions to use:

posterior_stats <- list(
  popef = extract_popef,
  trr = extract_trr,
  ratio = extract_ratio
)


## Summary statistic reliability models ----

if (do_summary_stats) {

  ## read trial-level data and compute subject-wise means

  d_summarystat <-
    fread(file.path(path_out, "spatial", "projections__stroop__rda__n_resamples100__demean_run__cv_allsess.csv")) %>%
    ## format for downstream analysis:
    rename(ridge = value.ridge, rda = value.rda, session = test_session, region = roi) %>%
    mutate(model_nm = "summarystat")
  cols_keep <- c("region", responses, "variable", "trial", "subj", "wave", "session", "model_nm")
  d_summarystat <- na.omit(d_summarystat[, ..cols_keep])
  d_summarystat_subj <-
    d_summarystat[, lapply(.SD, mean),
                  by = c("variable", "subj", "wave", "region", "session", "model_nm"),
                  .SDcols = responses]

  ## use means to compute population-level stroop coefficient

  d_summarystat_popef <- d_summarystat_subj %>%
    ## compute Stroop contrast
    pivot_longer(cols = all_of(responses), names_to = "response") %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(hilo = hi - lo) %>%
    group_by(region, session, model_nm, response) %>%
    summarize(
      est_mean = mean(hilo),
      est_sd = sd(hilo)
    )
  
  ## use means to compute summary-statistic reliability estimate

  d_summarystat_trr <- d_summarystat_subj %>%
    ## compute Stroop contrast
    pivot_longer(cols = all_of(responses), names_to = "response") %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(hilo = hi - lo) %>%
    ## compute correlation across waves
    pivot_wider(
      id_cols = c("subj", "region", "response", "session", "model_nm"),
      names_from = wave, values_from = hilo
    ) %>%
    group_by(region, response, session, model_nm) %>%
    summarize(r = cor(wave1, wave2))

  ## save

  fwrite(d_summarystat_subj, file.path(path_summary, "summarystat_subj_means.csv"))
  fwrite(d_summarystat_popef, file.path(path_summary, "summarystat_popef.csv"))
  fwrite(d_summarystat_trr, file.path(path_summary, "summarystat_trr.csv"))

}


## HBM model fit criteria ----

if (do_hbm_criteria) {

  criteria <- future_pmap_bind(
    model_info,
    read_summarize_hbm,
    base_path = file.path(path_out, "inferential"),
    atlas_nm = atlas_nm,
    sum_fun = pull_criteria
  )
  criteria <- reformat_for_downstream_analyses(criteria)
  fwrite(criteria, file.path(path_summary, "fit_criteria.csv"))

}


## HBM posteriors ----

if (do_hbm_posteriors) {

  for (statistic_i in seq_along(posterior_stats)) {

    samples <- future_pmap_bind(
      model_info,
      read_summarize_hbm,
      base_path = file.path(path_out, "inferential"),
      atlas_nm = atlas_nm,
      sum_fun = posterior_stats[[statistic_i]]
    )
    ## minor formatting changes for downstream analyses:
    samples <- reformat_for_downstream_analyses(samples)
    
    ## compute summary statistics on posterior distributions:
    summaries <- summarize_posteriors(samples)

    ## save to disk:
    saveRDS(samples, file.path(path_summary, paste0("posterior_samples_", names(posterior_stats)[statistic_i], ".RDS")))
    saveRDS(
      summaries,
      file.path(path_summary, paste0("posterior_summaries_", names(posterior_stats)[statistic_i], ".RDS"))
    )

    ## restore memory:
    rm(samples, summaries)
    gc()

  }

}
