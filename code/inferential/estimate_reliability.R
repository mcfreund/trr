# Hierarchical modeling for test-retest reliability
#
# 07/23/2022 update: now using only `brms` and fitting three different models.
#
# Author: Ruiqi Chen
#
# This script reads in the statistics and builds a hierarchical model for it.
# Note that currently the contrasts are set as (-0.5, 0.5) for ("lo", "hi") and
# (-0.5, 0.5) for ("wave1", "wave2"), so the fixed effect represents the difference
# between "hi" - "lo" or "wave2" - "wave1".
#
# We train the model using both maximum-likelihood method (via package `lme4`)
# and Bayesian MCMC method (via package `brms`). The `lme4` models mostly can't
# converge so we mainly use the results from `brms`. Note that the MCMC is very slow
# (~10h with 20 cores) and the parallelization often causes some jobs to fail.
# You need to check `is.na(fits_bayes)` and re-train the models that fail, which
# typically all finish in the second run.
#
# We save the fixed effects from the `lme4` models in a .csv file. For `brms` model,
# there are two ways:
#
# 1. Extract the effects of interest from the models by the function `pull_bayes_ef()`
#   in the script. Currently this function will save the point estimate, estimation
#   error, 2.5 and 97.5 percentile of the following "effects":
#  - all fixed effects (slope)
#  - all random effects & the error (standard deviation)
#  - "normalized" random effect of (hilo|subj)/error and (wave|subj)/error
#  - "normalized" fixed effect of abs(wave)/sd(y) and error/sd(y), where y is the
#     response (the input statistics we model).
# The extracted effects will be saved in a .csv file similar to that for `lme4` models.
#
# 2. Save all MCMC samples (with the estimated coefficients from each sample) from all
#   models as a list of 400 dataframes (of size 4000 * ~90) into an .rds file.
#
# The second way uses much more disk space (~1GB) but is more flexible so that you don't
# need to re-train the models (for ~10h) if you want to use another term that is not
# calculated by `pull_bayes_ef()`.

library(here)
library(tidyverse)
library(mfutils)
library(doParallel)
library(foreach)
library(brms)
library(parallel)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))
# source(here("code", "inferential", "_plotting.R"))


################### Command line parameters ######################

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 4)

# Help
if (length(args) == 0) {
  print(paste(
    "Usage: Rscript estimate_reliability.R",
    "[response_name = rda / uv] [n_cores (default = 4)]",
    "[roi_idx_first (default = 1)] [roi_idx_last (default = 2)]"
  ))
  q()
}

# Response name
response_names <- c(args[[1]])  # "rda" or "uv"
stopifnot(response_names[[1]] %in% c("rda", "uv"))

# Number of total cores
if (length(args) >= 2) n_cores <- strtoi(args[[2]]) else n_cores <- 4

# Subset of ROIs to use
if (length(args) >= 3) roi_idx_first <- strtoi(args[[3]]) else roi_idx_first <- 1
if (length(args) >= 4) roi_idx_last <- strtoi(args[[4]]) else roi_idx_last <- 2
stopifnot(roi_idx_first <= roi_idx_last)
roi_idx <- roi_idx_first:roi_idx_last


########################## Constants ##############################

n_core_brm <- 4  # Number of cores for parallelization within brm()
file_refit <- "on_change"  ## "on_change", "never", "always", see help(brm)
tasks <- "Stroop"
sessions <- c("baseline")
vterm <- switch(response_names[[1]],
  rda = "value.rda",
  ridge = "value.ridge",
  uv = "uv"
)
model_names <- c("no_lscov_symm")


## Input from ./code/spatial/multi...task.R:
fname <- here("out", "spatial",
  "projections__stroop__rda__n_resamples100__demean_run__cv_allsess.csv"
)

# Output path
out_path <- here("out", "inferential", "schaefer2018_17_400_fsaverage5")


# ROIs
rois <- readRDS(here("in", "rois.RDS"))
rois <- rois[roi_idx]

# Make sure we don't use more cores than available
stopifnot(n_core_brm <= n_cores)


###################### Load and prepare data #########################

# Select part of the ROIs
# Caution: old version of readr will parse the subject column incorrectly
#   if you don't supply col_types, resulting in missing subjects!
d <- read_csv(
  fname, col_types = list(
    test_session = "c",
    roi = "c",
    value.ridge = "d",
    value.rda = "d",
    variable = "c",
    trial = "i",
    auc_ridge = "d",
    auc_rda = "d",
    uv = "d",
    subj = "c",
    task = "c",
    wave = "c"
  )
)
stop_for_problems(d)
d <- d %>% filter(roi %in% .env$rois) %>% na.omit()


# Convert to wide form and create regressors

d_wide <-
  d %>%
  pivot_wider(id_cols = c(trial, subj, test_session, task, wave, variable),
    names_from = roi, values_from = .env$vterm) %>%
  mutate(
    ## coerce to factor:
    variable  = factor(variable, c("lo", "hi"), ordered = TRUE),
    wave = factor(wave, c("wave1", "wave2"), ordered = TRUE),
    ## create contrast coded numeric:
    variable_c = ifelse(variable == "hi", 0.5, -0.5),
    ## create dummy codes for chen2021 model:
    mean_wave1 = as.numeric(wave == "wave1"),
    mean_wave2 = as.numeric(wave == "wave2"),
    hilo_wave1 = variable_c * mean_wave1,
    hilo_wave2 = variable_c * mean_wave2,
    ## create dummy codes for 'flat' model:
    lo_wave1 = (variable == "lo") * mean_wave1,
    hi_wave1 = (variable == "hi") * mean_wave1,
    lo_wave2 = (variable == "lo") * mean_wave2,
    hi_wave2 = (variable == "hi") * mean_wave2
  )

# Fix naming problems with brms() formulas
input_for_bayes <- d_wide %>%
  filter(if_all(starts_with("17Networks"), ~ !is.na(.x))) %>%
  setNames(gsub("17Networks", "Networks", names(.)))
rois_bayes <- gsub("17Networks", "Networks", rois)

# Save some memory
rm(d, d_wide)


################## Fit Bayesian hierarchical models ###############

needs_refit <- enlist(combo_paste(model_names, "__", response_names, "__", sessions))

# CHPC only
session <- "baseline"
model_name <- "no_lscov_symm"
response_name <- response_names[[1]]

input_for_bayes_i <- input_for_bayes %>% filter(test_session == session)

## get formulas and out file prefix
model_info <- get_model_info(model_name = model_name, response_name = response_name, session = session)
formulas <-
  lapply(
    paste0(rois_bayes, model_info[["formula_string"]]),
    function(x) bf(as.formula(x), model_info[["formula_sigma"]])
  )
names(formulas) <- rois

## check:
print(get_prior(brmsformula(formulas[[1]]), input_for_bayes, family = student()))
if (FALSE) {
  t0 <- Sys.time()
  fit_check <- brm(
    formulas[[1]],
    input_for_bayes_i,
    cores = n_core_brm,
    family = student(),
    file = here("out", "spatial", "tmp")  # To get a sense of the size
    )
  t1 <- Sys.time()
  add_criterion(
    fit_check,
    criterion = c("loo", "waic", "bayes_R2"),
    file = here("out", "spatial", "tmp")
  )
  t2 <- Sys.time()
  if (TRUE) {
    cat(paste("Elapsed time to train a", model_name, "model: "))
    print(t1 - t0)
    cat("Elapsed time for evaluating: ")
    print(t2 - t1)
    cat("Total time: ")
    print(t2 - t0)
  }
}

print(paste0(" --------------- starting ", model_name, " ", response_name, " ", session, " --------------- "))

fit_model <- function(x) {
  out_subdir <- file.path(out_path, x)
  if (!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)
  out_file_name <- file.path(out_subdir, model_info$model_prefix)
  fit <- brm(
    formula = formulas[[x]],
    data = input_for_bayes_i,
    family = student(),
    cores = n_core_brm,
    save_model = out_file_name,
    file = out_file_name,
    file_refit = file_refit
  )
  add_criterion(
    fit,
    criterion = c("loo", "waic", "bayes_R2"),
    file = out_file_name
  )
}

fits_bayes <- mclapply(
  names(formulas),
  function(x) {
    tryCatch(expr = fit_model(x), error = function(e) NA)
  },
  mc.cores = min(length(formulas), n_cores %/% n_core_brm)
)

print(paste0(length(rois[is.na(fits_bayes)]), " needs rerun"))
needs_refit[[paste0(model_name, "__", response_name, "__", session)]] <- rois[is.na(fits_bayes)]

print(needs_refit)


# file_names <-
#   paste0(
#     file.path(out_path, rois[10],
#     c(
#       get_model_info(model_name = "full", response_name = "mv", session = "baseline")$model_prefix,
#       get_model_info(model_name = "no_lscov", response_name = "mv", session = "baseline")$model_prefix,
#       get_model_info(model_name = "no_lscov_symm", response_name = "mv", session = "baseline")$model_prefix
#       )
#     ),
#     "baseline.rds"
#     )
# fits <- lapply(file_names, readRDS)
# fits