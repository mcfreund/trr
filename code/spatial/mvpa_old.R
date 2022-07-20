# Compute trial-level multivariate statistics
#
# Author: Michael Freund
#
# 05/13/2022 updated by Ruiqi Chen:
#
# This script trains a linear classifier for each subject * task * region using the
# proactive & reactive sessions as training set and tests it on the baseline session.
# The value predicted for each trial is saved along with the classifier's AUC.
#
# The variable `classifier` specifies the type of classifier to use. "fda" will use
# the `fda()` function from the `mda` package for shrinkage-based LDA. "schafer_full"
# uses the `lda_schafer()` function from package `sparsediscrim` for LDA, where the
# covariance matrix is estimated by the method in Schafer and Strimmer (2005).
# "schafer_diag" also uses `lda_schafer()` but with parameter `lambda = 1` so that
# the off-diagonal elements of the estimated covariance matrix are forced to be 0.
#
# For "fda" classifier, the obtained statistic is the prediction for "variates";
# for "schafer_full" and "schafer_diag", the statistic is the log of the ratio
# between predicted posterior probability for high control over low control trial
# type for each trial in the baseline condition.
#
# To-do: the training returns a few NAs for "schafer_full" and "schafer_diag" and
# the reason is not clear yet.

library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(mda)
library(sparsediscrim)
library(pROC)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))


## input vars ----

variable <- "hilo_all"
classes <- c("lo", "hi")  ## -, +
tasks <- "Stroop"
train <- c("proactive", "reactive")
test <- c("baseline")
classifier <- "schafer_diag"  # "schafer_full", "schafer_diag" (ignoring covariances) or "fda"
shrinkage_factor <- 100  # Only used for "fda" classfier
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
n_cores <- 12
n_resamples <- 100

fname <- ifelse(grep("^schafer", classifier),
  paste0("projections__stroop__", classifier, "__n_resamples", n_resamples, ".csv"),
  paste0("projections__stroop__rda_lambda_", shrinkage_factor, "__n_resamples", n_resamples, ".csv")
)  # Output filename (under ./out/spatial/)


## execute ----

atlas <- get(atlas_nm)
rois <- unique(atlas$key[[roi_col]])
waves <- waves[do_waves]
n_classes <- length(classes)

## read trial-wise coefficients:
alltrials <- read_results(
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) "errts_trials_target_epoch.RDS",
  read_fun = readRDS,
  n_cores = n_cores
)

## read trial data:
## behav$hilo indicates high-demand (incongruent) vs low-demand (congruent) for only bias trialtypes
## behav$hilo_all indicates high-demand (incongruent) vs low-demand (congruent) for both bias and pc50 trialtypes
behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_alltasks.csv"), na.strings = c("", "NA"))
cols <- c("subj", "wave", "task", "session", "trialtype", variable, "trialnum")
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves, ..cols]

# # for dev/interactive use:
# subj_i <- 3
# task_i <- 1
# wave_i <- 1
# session_i <- 1
# roi_i <- 1

cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
allres <-
  foreach(task_i = seq_along(tasks), .inorder = FALSE, .combine = "rbind") %:%
  foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "rbind") %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "rbind") %dopar% {

    tryCatch({

    task_val <- tasks[task_i]
    subj_val <- subjs[subj_i]
    wave_val <- waves[wave_i]

    trials <- alltrials[grepl(paste0(wave_val, ".*", task_val, ".*", subj_val), names(alltrials))]
    stopifnot(length(trials) == 3)
    trials <- trials[sort(names(trials))]  ## sorts baseline, proactive, reactive
    names(trials) <- c("baseline", "proactive", "reactive")

    ## get good trials, class lables, and regress nuisance variance
    data_clean <- enlist(sessions)
    trial_idxs <- enlist(sessions)
    for (session_i in seq_along(trials)) {

      session_val <- sessions[session_i]

      ## subset trials (fmri) and behav (metadata) to single subject, session, wave, task
      trials_session <- trials[[session_val]]
      behav_session <- behav[subj == subj_val & wave == wave_val & session == session_val]
      behav_session <- behav_session[sort(trialnum)]
      stopifnot(nrow(behav_session) == nrow(trials_session))
      y <- behav_session[[variable]]  ## classes / outcome vector

      ## identify run1 vs run2 trials
      is_run1 <- seq_len(nrow(trials_session)) < (n_trialspr[paste0(task_val, "_", session_val)] + 1)

      ## identify bad/uninteresting trials
      ok_trials_idx <- which(rowSums(is.na(trials_session)) == 0)  ## identify non-missing/non-censored trials
      ttype_of_interest_idx <- which(!is.na(behav_session[[variable]]))
      idx <- intersect(ok_trials_idx, ttype_of_interest_idx)
      trial_idxs[[session_val]] <- idx

      ## subset
      trials_session_good <- trials_session[idx, ]
      y_good <- y[idx]
      is_run1_good <- is_run1[idx]
      run_labels <- ifelse(is_run1_good, "run1", "run2")

      ## regress nuisance variance:
      ## each vertex may have different mean in different scanning runs, independent of the different trialtypes
      ## here we estimate those means and remove them from each vertex
      X <- cbind(indicator(y_good) * is_run1_good, indicator(y_good) * !is_run1_good)
      mu <- coef(.lm.fit(x = X, y = trials_session_good))  ## yeilds class means per run
      mu_bar <- average(t(mu), rep(c("run1", "run2"), each = n_classes))  ## mean of class means per run (vertex by run)
      mu_bar <- tcrossprod(indicator(run_labels), mu_bar)  ## expand to match dims of data
      trials_session_good_c <- trials_session_good - mu_bar  ## center
      trials_session_good_c <- t(trials_session_good_c)
      ## store class labels as colnames:
      colnames(trials_session_good_c) <- paste0(y_good, "__", run_labels, "__", session_val)

      data_clean[[session_i]] <- trials_session_good_c

    }

    ## build set of trial indices for stratified resampling:
    ids_train <- unlist(lapply(data_clean, colnames)[train], use.names = FALSE)
    set.seed(0)  ## OK to use same seed across subjs/sessions etc (b/c trial positions randomized)
    resampled_idx <- resample_idx(ids_train, n_resamples)  ## rows: n_resamples; cols: n_trial in each training model
    y_train <- abind(strsplit(colnames(resampled_idx), "__"), along = 0)[, 1]  ## extract class labels
    y_train <- factor(y_train, levels = classes)
    ids_test <- unlist(lapply(data_clean, colnames)[test], use.names = FALSE)
    y_test <- abind(strsplit(ids_test, "__"), along = 0)[, 1]  ## extract class labels
    y_test <- factor(y_test, levels = classes)

    ## segment images:
    data_clean_rois <- lapply(data_clean, parcellate, atlas, col_roi = roi_col)

    ## loop over ROIs and train/test models:
    projs <- enlist(rois)
    for (roi_i in seq_along(rois)) {

      print(paste("Now processing", task_val, "task for subject", subj_val, wave_val,
        "in region", roi_i, ":", rois[[roi_i]]))

      ## extract roi_i for each session
      data_clean_roi <- lapply(data_clean_rois, "[[", rois[roi_i])  ## list of matrices of vertex by trial

      ## remove vertices with no BOLD variance
      is_good_vertex_session <- lapply(data_clean_roi, function(x) !is_equal(Var(x, 1), 0))
      is_good_vertex <- Reduce("&", is_good_vertex_session)  ## intersection
      stopifnot(mean(is_good_vertex) > 3/4)  ## stop if less than 3/4 vertices in ROI do not have signal
      data_clean_roi_goodverts <- lapply(data_clean_roi, function(x) x[is_good_vertex, ])

      ## data
      d_test <- t(do.call(cbind, data_clean_roi_goodverts[test]))
      d_train <- t(do.call(cbind, data_clean_roi_goodverts[train]))

      ## Add column names to data for lda_schafer()
      colnames(d_test) <- paste0("v", seq_len(ncol(d_test)))
      colnames(d_train) <- paste0("v", seq_len(ncol(d_train)))

      ## train
      fits <- lapply(
        seq_len(nrow(resampled_idx)),
        function(resample_i) {
          .idx <- resampled_idx[resample_i, ]
          if (classifier == "schafer_full") {
            return(lda_schafer(d_train[.idx, ], y_train))
          } else if (classifier == "schafer_diag") {
            return(lda_schafer(d_train[.idx, ], y_train, lambda = 1))
          } else {
            return(fda(y_train ~ d_train[.idx, ], method = gen.ridge, lambda = shrinkage_factor))
          }
        }
        #mc.cores = n_cores
      )
      ## test
      projs_i <- vapply(
        fits,
        function(.x, .newdata) if (grep("^schafer", classifier)) {
          tmp <- predict(.x, newdata = .newdata, type = "prob")  # Note: using "score" seems to be wrong
          return(as.numeric(t(log(tmp[classes[[2]]] + 1e-6) - log(tmp[classes[[1]]] + 1e-6))))
        } else {
          return(predict(.x, newdata = .newdata, type = "variates"))
        },
        numeric(nrow(d_test)),
        .newdata = d_test
      )
      #image(projs_i)
      #proj_bar <- apply(projs_i, 1, median)
      proj_bar <- rowMeans(projs_i)

      ## summarize performance with auc
      rocobj <- roc(y_test, proj_bar, direction = "<", levels = classes)

      projs[[roi_i]] <- data.table(value = proj_bar, variable = y_test, trial = trial_idxs[[test]], auc = rocobj$auc)

    }

    res <- rbindlist(projs, idcol = "roi")
    res[, ":=" (subj = subj_val, task = task_val, wave = wave_val)]  ## add subj/sess info

    res  ## return

    }, error = function(e) {
      print("")
      print(e)
      print("")
      return(NA)
    })
  }
stopCluster(cl)

fwrite(allres, here("out", "spatial", fname))