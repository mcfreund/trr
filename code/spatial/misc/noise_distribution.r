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
divnorm_vertex <- TRUE
divnorm_trial <- FALSE
demean_trial <- FALSE
classes <- c("lo", "hi")  ## -, +
tasks <- "Stroop"
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
n_cores <- 6

## out file name
file_name <- paste0(
  "noisepca__stroop",
  switch(divnorm_vertex + 1, "", "__divnorm_vertex"),
  switch(divnorm_trial + 1, "", "__divnorm_trial"),
  switch(demean_trial + 1, "", "__demean_trial"),
  ".csv"
  )
print(noquote(paste0("will save to file: ", file_name)))

## atlas info and other constants
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

## for dev/interactive use:
if (FALSE) {
  subj_i <- which(subjs == "DMCC1596165")
  task_i <- 1
  wave_i <- 2
  session_i <- 1
  roi_i <- 130
}


## execute ----

cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
allres <-
  foreach(task_i = seq_along(tasks), .inorder = FALSE, .combine = "c") %:%
  foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "c") %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %dopar% {

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
        sdev <- sqrt(Var(eps, 2))
        trials_session_c <- sweep(trials_session_c, 1, sdev, "/")
      }

      data_clean[[session_i]] <- trials_session_c

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


    projs <- enlist(rois)
    for (roi_i in seq_along(rois)) {

      ## extract roi_i for each session
      data_clean_roi <- lapply(data_clean_rois, "[[", rois[roi_i])  ## list of matrices of vertex by trial

      ## remove vertices with no BOLD variance
      is_good_vertex_session <- lapply(data_clean_roi, function(x) !is_equal(Var(x, 1, na.rm = TRUE), 0))
      is_good_vertex <- Reduce("&", is_good_vertex_session)  ## intersection
      stopifnot(mean(is_good_vertex) > 1/4)  ## stop if less than 1/4 vertices in ROI do not have signal
      data_clean_roi_goodverts <- lapply(data_clean_roi, function(x) x[is_good_vertex, ])

      l <- enlist(sessions)
      for (session in sessions) {
      ## within-condition variance
        x <- data_clean_roi_goodverts[[session]]
        levs <- unique(colnames(x))
        res <- enlist(levs)
        for (lev_i in seq_along(levs)) {
          x_k <- x[, colnames(x) %in% levs[lev_i]]
          if (any(is.na(x_k))) x_k <- x_k[, -which(is.na(colSums(x_k)))]
          x_k <- t(x_k)  ## rows are obs (trials), cols are spatial dims (vertices)
          n_vert <- ncol(x_k)
          pca <- prcomp(x_k, scale = FALSE)
          unif <- rep(1 / sqrt(n_vert), n_vert)
          cossim <- abs(c(unif %*% pca$rotation))
          res[[lev_i]] <- data.table(proj = cossim, eigval = pca$sdev^2)
        }
        l[[session]] <- rbindlist(res, idcol = "level")
      }

      projs[[roi_i]] <- rbindlist(l, idcol = "session")

    }

    d <- rbindlist(projs, idcol = "roi")
    d[, ":=" (subj = subj_val, task = task_val, wave = wave_val)]  ## add subj/sess info

    list(d)

}
stopCluster(cl)

out <- rbindlist(allres)
fwrite(out, here("out", "spatial", file_name))