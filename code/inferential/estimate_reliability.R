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
library(ggsegSchaefer)
library(doParallel)
library(foreach)
library(brms)
library(parallel)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))
# source(here("code", "inferential", "_plotting.R"))


########################## Constants ##############################

n_core_brm <- 4  # Number of cores for parallelization within brm()
roi_idx <- core32  ## which ROIs to use
n_cores <- 16
file_refit <- "on_change"  ## "never", "always", see help(brm)
tasks <- "Stroop"
sessions <- c("baseline", "proactive", "reactive")
response_names <- c("rda")  # Note: currently must be c("rda"), c("ridge") or c("uv")
vterm <- switch(response_names[[1]],
  rda = "value.rda",
  ridge = "value.ridge",
  uv = "uv"
)
model_names <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma")


## Input from ./code/spatial/multi...task.R:
fname <- here("out", "spatial",
  "projections__stroop__rda__n_resamples100__demean_run__cv_allsess.csv"
)

# Output path
out_path <- here("out", "inferential", "wo_divnorm")

# Atlas
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"

# ROIs
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
}
rois <- rois[roi_idx]

# Make sure we don't use more cores than available
stopifnot(n_core_brm <= n_cores)


###################### Load and prepare data #########################

# Select part of the ROIs
d <- read_csv(fname) %>% filter(roi %in% .env$rois) %>% na.omit()


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



################## Fit Bayesian hierarchical models ###############

needs_refit <- enlist(combo_paste(model_names, "__", response_names, "__", sessions))

# Debugging
if (FALSE) {
  session <- "baseline"
  model_name <- "fixed_sigma"
  response_name <- response_names[[1]]
}

for (session in sessions) {
  for (model_name in model_names) {
    for (response_name in response_names) {

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
      time_start <- Sys.time()
      print(time_start)

      fits_bayes <- mclapply(
        names(formulas),

        function(x) {
          tryCatch(

            expr = {
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
            },

            error = function(e) {
              print("------------------------")
              print("ERROR:")
              print(e)
              print("------------------------")
              return(NA)
            }

          )
        },

        mc.cores = min(length(formulas), n_cores %/% n_core_brm)
      )

      print(Sys.time() - time_start)
      print(paste0(length(rois[is.na(fits_bayes)]), " needs rerun"))
      needs_refit[[paste0(model_name, "__", response_name, "__", session)]] <- rois[is.na(fits_bayes)]

    }
  }
}
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