# Compute trial-level multivariate statistics
#
# Author: Michael Freund
#
#
# 07/21/2022 update:
#
# Problem: ldf() use the `.covpooled` of lda models, which cannot be estimated when data contain NAs
# See https://github.com/cran/klaR/blob/master/R/rda.R#L526-L531.
#
# Two reason:
#
# 1. Bad trials, e.g., Subject 7 (448347) wave 2 Stroop baseline trial 93 & 94.
#  This should be fixed after the update of `resampleing_trials.r`.
#
# 2. Divide by 0 error when normalizing vertices with no fluctuation.
#  This should be fixed by the update in `sdev <- sqrt(Var(eps, 2)) + 1e-6`.
#
# 07/20/2022 update:
#
# Use `klaR`` instead of the (problematic) `sparsediscrim`` for classification.
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
# The training returns a few NAs for "schafer_full" and "schafer_diag" and
# the reason is not clear yet.

library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(mda)
library(klaR)
library(pROC)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))


## input vars ----

variable <- "hilo_all"
divnorm_vertex <- TRUE
divnorm_trial <- FALSE
demean_trial <- FALSE
classes <- c("lo", "hi")  ## -, +
tasks <- "Stroop"
classifier <- "rda"  ## or "ridge"
shrinkage_factor_ridge <- 100
shrinkage_factor_rda <- 0.25
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
glm_nm <- "null_2rpm"
resid_type <- "errts"
n_cores <- 18
do_waves <- c(1, 2)
subjs <- switch(toString(do_waves),
  "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all
)
input_fname <- here("in", "behav",
  paste0("behavior-and-events_wave", do_waves[1], do_waves[2], "_alltasks.csv")
)
file_name_resamples <- here("out", "spatial",
  paste0("trialidx_stroop_congruency_wave", do_waves[1], do_waves[2], ".RDS")
)

## atlas info and other constants
atlas <- get(atlas_nm)
rois <- unique(atlas$key[[roi_col]])
waves <- waves[do_waves]
n_classes <- length(classes)

## read resampled trial indices
resamples <- readRDS(file_name_resamples)
n_resamples <- nrow(resamples[[1]])
print(noquote(paste0("num resamples: ", n_resamples)))

## out file name
file_name <- paste0(
  "projections__stroop__", classifier, "__n_resamples", n_resamples,
  switch(divnorm_vertex + 1, "", "__divnorm_vertex"),
  switch(divnorm_trial + 1, "", "__divnorm_trial"),
  switch(demean_trial + 1, "", "__demean_trial"),
  "__cv_allsess_wave", do_waves[1], do_waves[2], ".csv"
)

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
behav <- fread(input_fname, na.strings = c("", "NA"))
cols <- c("subj", "wave", "task", "session", "trialtype", variable, "trialnum")
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves, ..cols]


## utilities ----

## convenience function for use with klaR::rda()
ldf <- function(object, newdata, class_names = c("hi", "lo")) {
  ## class_names: positive class first
  r <- object$regularization
  if (r["lambda"] < 1) stop("Not configured for QDA.")
  sigmahat <- object$covpooled
  scaled_identity <- mean(diag(sigmahat)) * diag(nrow(sigmahat))
  sigmahat_reg <- (1 - r["gamma"]) * sigmahat + r["gamma"] * scaled_identity
  w <- solve(sigmahat_reg) %*% (object$means[, class_names] %*% rbind(1, -1))
  w <- w / sqrt(sum(w^2))  ## scale to unit length
  newdata %*% w
}


## for dev/interactive use:
if (FALSE) {
  subj_i <- which(subjs == "448347")
  task_i <- 1
  wave_i <- 2
  session_i <- 1
  roi_i <- 1
}


## execute ----

cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
allres <-
  foreach(task_i = seq_along(tasks), .inorder = FALSE, .combine = "c") %:%
  foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "c") %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %dopar% {

    tryCatch(
      {

      task_val <- tasks[task_i]
      subj_val <- subjs[subj_i]
      wave_val <- waves[wave_i]

      trials <- alltrials[grepl(paste0(wave_val, ".*", task_val, ".*", subj_val), names(alltrials))]
      stopifnot(length(trials) == 3)
      trials <- trials[sort(names(trials))]  ## sorts baseline, proactive, reactive
      names(trials) <- c("baseline", "proactive", "reactive")

      ## get good trials, class lables, and regress nuisance variance
      data_clean <- mfutils::enlist(sessions)
      trial_idxs <- mfutils::enlist(sessions)
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
        run_labels <- ifelse(is_run1, "run1", "run2")

        ## regress nuisance variance:
        ## each vertex may have different mean in different scanning runs, independent of the different trialtypes
        ## here we estimate those means and remove them from each vertex
        X <- cbind(indicator(y_good) * is_run1_good, indicator(y_good) * !is_run1_good)
        mu <- coef(.lm.fit(x = X, y = trials_session_good))  ## yeilds class means per run
        ## mean of class means per run (vertex by run):
        mu_bar <- average(t(mu), rep(c("run1", "run2"), each = n_classes))
        mu_bar <- tcrossprod(indicator(run_labels), mu_bar)  ## expand to match dims of (non-subsetted) data
        trials_session_c <- trials_session - mu_bar  ## center
        trials_session_c <- t(trials_session_c)
        ## store class labels as colnames:
        #colnames(trials_session_c) <- paste0(y, "__", run_labels, "__", session_val)
        colnames(trials_session_c) <- y

        if (divnorm_vertex) {
          ## divisive normalize each vertex by residual sdev over trials? (univariate prewhiten)
          eps <- resid(.lm.fit(x = indicator(behav_session$trialtype[idx]), y = t(trials_session_c[, idx])))
          sdev <- sqrt(Var(eps, 2)) + 1e-6  # To prevent NaN for vectices with no BOLD
          trials_session_c <- sweep(trials_session_c, 1, sdev, "/")
        }

        data_clean[[session_i]] <- trials_session_c  # Note: remove trials by resampled idx, not here

      }

      ## segment images:
      data_clean_rois <- lapply(data_clean, parcellate, atlas, col_roi = roi_col)

      ## preprocess each parcel:
      if (divnorm_trial | demean_trial) {
        for (session_i in seq_along(trials)) {
          data_clean_rois[[session_i]] <-
            lapply(data_clean_rois[[session_i]], scale, center = demean_trial, scale = divnorm_trial)
        }
      }


      projs_all <- mfutils::enlist(sessions)
      for (test in sessions) {

        train <- setdiff(sessions, test)
        #resamples_test <- resamples[[paste0(subj_val, "__", wave_val, "__", test)]]
        resamples_train <- resamples[paste0(subj_val, "__", wave_val, "__", train)]
        names(resamples_train) <- train

        ## loop over ROIs and train/test models:
        projs <- mfutils::enlist(rois)
        for (roi_i in seq_along(rois)) {

          print(paste("Testing on", task_val, test, "for subject", subj_val, wave_val,
            "region", roi_i, ":", rois[[roi_i]]))

          ## extract roi_i for each session
          data_clean_roi <- lapply(data_clean_rois, "[[", rois[roi_i])  ## list of matrices of vertex by trial

          ## remove vertices with no BOLD variance
          is_good_vertex_session <- lapply(data_clean_roi, function(x) !is_equal(Var(x, 1, na.rm = TRUE), 0))
          is_good_vertex <- Reduce("&", is_good_vertex_session)  ## intersection
          stopifnot(mean(is_good_vertex) > 1/4)  ## stop if more than 1/4 vertices in ROI do not have signal
          data_clean_roi_goodverts <- lapply(data_clean_roi, function(x) x[is_good_vertex, ])

          ## data
          ## Add column names to data for lda_schafer()
          d_test <- t(do.call(cbind, data_clean_roi_goodverts[test]))
          colnames(d_test) <- paste0("v", seq_len(ncol(d_test)))
          l_train <- data_clean_roi_goodverts[train]
          for (i in seq_along(l_train)) {
            l_train[[i]] <- t(l_train[[i]])
            colnames(l_train[[i]]) <- paste0("v", seq_len(ncol(l_train[[i]])))
          }

          ## train

          fits <- lapply(
            seq_len(n_resamples),
            function(i) {
              idx1 <- resamples_train[[1]][i, ]
              idx2 <- resamples_train[[2]][i, ]
              X <- rbind(l_train[[1]][idx1, ], l_train[[2]][idx2, ])
              y <- rownames(X)
              fit_ridge <- mda::fda(y ~ X, method = gen.ridge, lambda = shrinkage_factor_ridge)
              fit_rda <- klaR::rda(x = X, grouping = y, gamma = shrinkage_factor_rda, lambda = 1)
              return(list(ridge = fit_ridge, rda = fit_rda))
            }
          )

          ## test
          projs_i <- lapply(
            fits,
            function(.x, .newdata) {
              fit_ridge <- .x$ridge
              fit_rda <- .x$rda
              tmp <- predict(fit_ridge, newdata = .newdata, type = "distances")
              res_ridge <- tmp[, "lo"] - tmp[, "hi"]
              res_rda <- ldf(fit_rda, newdata = .newdata, class_names = c("hi", "lo"))  ## positive, negative
              return(cbind(ridge = c(res_ridge), rda = c(res_rda)))
            },
            .newdata = d_test
          )
          projs_i <- abind(projs_i, rev.along = 0)
          #image(projs_i[, 1, ])
          #image(projs_i[, 2, ])
          proj_bar <- rowMeans(projs_i, dims = 2)

          ## bind into single dataframe:
          y_test <- rownames(d_test)
          projs[[roi_i]] <-
            data.table(
              value = proj_bar,
              variable = y_test,
              trial = seq_along(proj_bar),  ## trial_idxs[[test]]
              ## summarize performance with auc:
              auc_ridge = roc(y_test, proj_bar[, "ridge"], direction = "<", levels = classes)$auc,
              auc_rda = roc(y_test, proj_bar[, "rda"], direction = "<", levels = classes)$auc,
              uv = rowMeans(d_test)  ## univariate projections
              )

        }

        res <- rbindlist(projs, idcol = "roi")
        res[, ":=" (subj = subj_val, task = task_val, wave = wave_val)]  ## add subj/sess info

        projs_all[[test]] <- res


      }

      result <- rbindlist(projs_all, idcol = "test_session", fill = TRUE)
      list(result)

    },

    error = function(e) {
      msg <- paste("Error in task", tasks[task_i],
        "subject", subjs[subj_i], "wave", waves[wave_i], ":", e)
      print("")
      print(msg)
      print("")
      list(msg)
    }

  )
}
stopCluster(cl)

out <- rbindlist(allres)
fwrite(out, here("out", "spatial", file_name))