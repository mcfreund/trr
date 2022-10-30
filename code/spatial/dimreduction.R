library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(ggplot2)
source(here("code", "_constants.R"))
source(here("code", "_funs.R"))

sessions <- "baseline"
variable <- "hilo_all"
classes <- c("lo", "hi")  ## -, +
tasks <- "Stroop"
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
glm_nm <- "null_2rpm"
resid_type <- "errts"
n_cores <- 18
demean_run <- TRUE
do_waves <- c(1, 2)
subjs <- switch(toString(do_waves),
  "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all
)
input_fname <- here("in", "behav",
  paste0("behavior-and-events_wave", do_waves[1], do_waves[2], "_alltasks.csv")
)
path_base <- "/data/nil-external/ccp/chenr/trr"
file_name_resamples <- file.path(path_base, "out", "spatial",
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
  "signalnoise__stroop",
  switch(demean_run + 1, "", "__demean_run"),
  "__baseline__wave", do_waves[1], do_waves[2], ".csv"
)

## read trial-wise coefficients:
alltrials <- read_results(
  waves = waves, tasks = tasks, sessions = sessions, subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) "errts_trials_target_epoch.RDS",
  read_fun = readRDS,
  n_cores = n_cores,
  path_base = file.path(path_base, "out", "timeseries")
)


## read trial data:
## behav$hilo indicates high-demand (incongruent) vs low-demand (congruent) for only bias trialtypes
## behav$hilo_all indicates high-demand (incongruent) vs low-demand (congruent) for both bias and pc50 trialtypes
behav <- fread(input_fname, na.strings = c("", "NA"))
cols <- c("subj", "wave", "task", "session", "trialtype", variable, "trialnum")
behav <- behav[task %in% tasks & session %in% sessions & wave %in% waves, ..cols]


## utilities ----


## for dev/interactive use:
if (FALSE) {
  subj_i <- 1
  task_i <- 1
  wave_i <- 1
  session_i <- 1
  roi_i <- 127
}



cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
allres <-
  foreach(task_i = seq_along(tasks), .inorder = FALSE, .combine = "c") %:%
  foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "c") %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %:%
  foreach(session_i = seq_along(sessions), .inorder = FALSE, .combine = "c") %dopar% {

    task_val <- tasks[task_i]
    subj_val <- subjs[subj_i]
    wave_val <- waves[wave_i]
    session_val <- sessions[session_i]

    trials <- alltrials[[paste0(wave_val, "_", task_val, "_", session_val, "_", subj_val)]]
    meta <- behav[subj == subj_val & wave == wave_val & session == session_val & task == task_val]
    y <- meta[[variable]]  ## classes / outcome vector
    rownames(trials) <- y

    ## identify run1 vs run2 trials
    is_run1 <- seq_len(nrow(trials)) < (n_trialspr[paste0(task_val, "_", session_val)] + 1)

    ## identify bad/uninteresting trials
    ok_trials_idx <- which(rowSums(is.na(trials)) == 0)  ## identify non-missing/non-censored trials
    ttype_of_interest_idx <- which(!is.na(meta[[variable]]))
    idx <- intersect(ok_trials_idx, ttype_of_interest_idx)

    ## subset
    trials_good <- trials[idx, ]
    y_good <- y[idx]
    is_run1_good <- is_run1[idx]
    run_labels <- ifelse(is_run1_good, "run1", "run2")

    ## regress run variance
    if (demean_run) {
      X <- cbind(indicator(y_good) * is_run1_good, indicator(y_good) * !is_run1_good)
      mu <- coef(.lm.fit(x = X, y = trials_good))  ## yeilds class means per run
      ## mean of class means per run (vertex by run):
      mu_bar <- average(t(mu), rep(c("run1", "run2"), each = n_classes))
      mu_bar <- tcrossprod(indicator(run_labels), mu_bar)  ## expand to match dims of (non-subsetted) data
      trials_good <- trials_good - mu_bar  ## center
    }

    ## segment images:
    parcs <- parcellate(t(trials_good), atlas, col_roi = roi_col)

    ## loop over ROIs and train/test models:
    res <-  mfutils::enlist(rois)
    for (roi_i in seq_along(rois)) {
      ## extract roi_i for each session
      parc <- parcs[[roi_i]]  ## vertex by trial
      ## remove vertices with no BOLD variance
      is_good_vertex <- !is_equal(Var(parc, 1, na.rm = TRUE), 0)
      b <- parc[is_good_vertex, ]
      n_vert <- nrow(b)

      unif <- cbind(rep(1, n_vert) / sqrt(n_vert))  ## uniform dimension, scaled to unit variance

      ## alignment between signal and unif
      signal_vec <- average(b, y_good) %*% cbind(c(1, -1)) / 2
      ssq_signal <- sum(signal_vec^2)
      signal_vec_scaled <- signal_vec / sqrt(ssq_signal)  ## for projectng noise vectors

      ## get noise distribution using pooled within-class covariance
      s <- vector("list", n_classes)
      for (i in seq_along(classes)) s[[i]] <- cov(t(b)[y_good == classes[i], ])
      sbar <- Reduce("+", s) / 2
      pca <- eigen(sbar)
      eigvecs <- pca$vectors %*% diag(sign(colMeans(pca$vectors)))  ## reflect axes to have positive mean
      eigvals <- pca$values
      eigvecs_unscaled <- scale(eigvecs, center = FALSE, scale = 1/eigvals)
      d <- data.frame(
        dimension = seq_len(n_vert),
        eigvals = eigvals,  ## variance (length^2) of each noise dim
        cossim_noise_unif = crossprod(eigvecs, unif),  ## cosine similarity between uniform and eigenvector
        proj_noise_unif = crossprod(eigvecs_unscaled, unif),  ## variance along uniform dim per noise dimension
        proj_noise_unif_scaled = crossprod(eigvecs_unscaled, unif) / sum(eigvals),  ## as proportion of total noise
        cossim_signal_noise = abs(crossprod(eigvecs, signal_vec_scaled)),
        proj_signal_noise = abs(crossprod(eigvecs, signal_vec)),  ## length of signal on each noise dim
        proj_signal_noise_scaled = crossprod(eigvecs, signal_vec)^2 / ssq_signal,  ## proportion of var on noise dim
        snr = proj_signal_noise^2 / eigvals,
        ssq_signal = ssq_signal,  ## these are sacalars so will be duplicated in output
        cossim_signal_unif = crossprod(signal_vec / sqrt(ssq_signal), unif)
      )

      res[[roi_i]] <- d

    }
    out <- rbindlist(res, idcol = "roi")
    out$subj <- subj_val
    out$wave <- wave_val
    out$session <- session_val

    return(list(out))

}
stopCluster(cl)

out <- rbindlist(allres)
fwrite(out, here("out", "spatial", file_name))


sum(d$proj_signal_noise^2) - ssq_signal
