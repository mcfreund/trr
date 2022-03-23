## setup ----

library(here)
library(mfutils)
library(dplyr)
library(tibble)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(ggplot2)
library(purrr)
library(gifti)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))

subjs <- subjs_wave12_all
do_waves <- c(1, 2)
resid_type <- "errts"

## iterators for dev:
# wave_i <- 1
# task_i <- 4
# session_i <- 1
# subj_i <- 1


## run ----

## 1. build averaging matrix A, of dimension trial by TR. ----
## when applied to BOLD time series E (TR by vertex), A aggregates across target TRs within each trial.
## that is, Z = A %*% E  is a trial-by-vertex matrix.
## To do this, rows (TRs) of A must be >0 when the TR belongs to the target window of a given trial, 0 otherwise.
## Further, rows must sum to one.

waves <- waves[do_waves]

A_list <- enlist(combo_paste(waves, tasks, sessions, subjs))

for (wave_i in seq_along(waves)) {
  for (task_i in seq_along(tasks)) {
    for (session_i in seq_along(sessions)) {
  
      name_wave_i <- waves[wave_i]
      name_task_i <- tasks[task_i]
      name_session_i <- sessions[session_i]
      dir_image <- file.path("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir_image[wave_i], "fMRIPrep_AFNI_ANALYSIS")
      dir_evts <- file.path("/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS", wavedir_evts[wave_i])

      n_tr <- n_trs[paste0(name_task_i, "_", name_session_i)]
      n_trial <- n_trialspr[paste0(name_task_i, "_", name_session_i)]*2
      
      ## read in stim times, censor list (all subjs):
      
      fname.stimtimes <-
        file.path(
          dir_evts, subjs, "evts", 
          paste0(subjs, "_", name_task_i, "_", name_session_i, "_allTrials.txt")
          )
      stimtimes <- lapply(fname.stimtimes, function(x) as.matrix(fread(x, header = FALSE, sep = " ")))
      
      fname.censor <-
        file.path(dir_image, subjs,  "INPUT_DATA", name_task_i, name_session_i, "movregs_FD_mask.txt")
      censor_list <- lapply(fname.censor, fread, header = FALSE, colClasses = "integer")
      
      
      ## build averaging matrix A:
      
      for (subj_i in seq_along(subjs)) {
        
        onsets <- c(
          stimtimes[[subj_i]][1, ], 
          stimtimes[[subj_i]][2, ] + sec_tr*n_tr/2    ## shift run 2 onsets by duration (in seconds) of run 1
          )
        
        if (sum(duplicated(onsets)) > 0) {
          print(
            paste0(
              "identical onsets: ", 
              subjs[subj_i], " ", name_task_i, " ", name_session_i, ", ", sum(duplicated(onsets))
              )
            )
        }
        overlapping_trials <- which(duplicated(onsets) | duplicated(onsets, fromLast = TRUE))
        intervals <- floor(onsets / 1.2) + 1  ## gives image (interval) to which onset belongs
  
        
        ## build indicator matrix (dummy):
        
        dummy <- matrix(0, nrow = n_tr, ncol = n_trial, dimnames = list(tr = NULL, trial = NULL))
        for (target_tr_i in seq_along(target_trs[[name_task_i]])) {
          
          ## assign unique value to each onset (i.e., each trial)
          
          target_tr_val <- target_trs[[name_task_i]][target_tr_i] + 1 ## plus 1 because of TENTzero/afni 0 based index
          trial_nums <- seq_along(intervals)

          tr_grid <- numeric(n_tr)
          tr_grid[intervals + target_tr_val] <- trial_nums
          
          ## set previous, so built in model.mat; for censoring overlapping trials
          tr_grid[which(tr_grid == overlapping_trials[2]) - 1] <- overlapping_trials[1]
  
          dummy <- dummy + model.matrix(~ as.factor(tr_grid))[, -1]  ## remove intercept
          
        }
        
        
        ## censor:
        
        dummy[, overlapping_trials] <- 0  ## censor overlapping trials
        dummy[censor_list[[subj_i]] < 1, ] <- 0  ## censor TRs
        
        ## scale:
        
        divisor <- colSums(dummy)
        dummy[, divisor < 1] <- NA ## when this is zero, all TRs of the trial were censored due to motion
        if (any(divisor < 1)) print(paste0("trial censored: ", subjs[subj_i], " ", name_task_i, " ", name_session_i))  ## just to look
        A <- sweep(dummy, 2, divisor, "/")  ## normalize elements of each column (trial) to sum to one
        
        ## save:
        
        nm <- paste0(waves[wave_i], "_", name_task_i, "_", name_session_i, "_", subjs[subj_i])
        A_list[[nm]] <- A
        # image(A)  ## view
        
      }
    }
  }
}




# name_subj_i <- subjs[subj_i]
# name_wave_i <- waves[wave_i]
# name_task_i <- tasks[task_i]

cl <- makeCluster(20, type = "FORK")
registerDoParallel(cl)
res <- 
  foreach(subj_i = seq_along(subjs), .inorder = FALSE) %:%
  foreach(task_i = seq_along(tasks), .inorder = FALSE) %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE) %:%
  foreach(session_i = seq_along(sessions), .inorder = FALSE) %dopar% {

    ## get constants, averaging matrix:

    name_subj_i <- subjs[subj_i]
    name_wave_i <- waves[wave_i]
    name_task_i <- tasks[task_i]    
    name_session_i <- sessions[session_i]
    name_glm <- paste0(name_session_i, "_", "null_2rpm_", waves[wave_i])
    dir_glm <- here("out", "timeseries", name_subj_i, "RESULTS", name_task_i, name_glm)
    n_tr <- n_trs[paste0(name_task_i, "_", name_session_i)]
    n_trial <- n_trialspr[paste0(name_task_i, "_", name_session_i)]*2
    nm <- paste0(waves[wave_i], "_", name_task_i, "_", name_session_i, "_", subjs[subj_i])
    A <- A_list[[nm]]
    
    ## load data:
    
    eps_name <- here(dir_glm, paste0(resid_type, "_", c("L", "R"), ".gii"))  ## L, R
    if (any(!file.exists(eps_name))) stop("no file!")
    E <- cbind(read_gifti2matrix(eps_name[1]), read_gifti2matrix(eps_name[2]))  ## L, R
    
    ## check for unexpected dimensions:
    dims_bad <- any(dim(E) != c(n_tr, n_vert))
    if (dims_bad) stop ("bad dims: error time-series")
    
    ## average and save:
    
    B <- crossprod(A, E)  ## trial by vertex matrix B
    # image(B)  ## view
    saveRDS(B, here(dir_glm, paste0(resid_type, "_trials_target_epoch.RDS")))
    saveRDS(A, here(dir_glm, paste0(resid_type, "_averaging_matrix.RDS")))
}
stopCluster(cl)








## misc spot checks ----



## <<<<<< BEGIN CHECKING DUMMY MATRIX >>>>>>

# read in design matrices (all subjs):
# fname.xmat <-
#   file.path(
#     dir_image, subjs,
#     "1TRpK_SURFACE_RESULTS", name_task_i, paste0(name_session_i, "_ON_MIXED_censored"), "X.xmat_L.1D"
#     )
# 
# ## NEED TO FORCE PYTHON2 USE
# # https://stackoverflow.com/questions/7237415/python-2-instead-of-python-3-as-the-temporary-default-python
# # mkdir ~/bin
# # PATH=~/bin:$PATH
# # ln -s /usr/bin/python2 ~/bin/python
# # rm ~/bin/python
# X <- read_xmat(fname.xmat[subj_i])
# X <- X > 0.9  ## threshold (b/c some non-'active' TRs will nevertheless have >0 value, though quite small)
# class(X) <- "numeric"
# is_censored <- is_equal(rowSums(X), 0)
# sum(is_censored)
# xmat_intervals <- unlist(apply(X[, c("ON_TRIALS#8", "ON_TRIALS#9", "ON_TRIALS#10")] > 0, 2, which), use.names = FALSE)
# onset_intervals <- unlist(apply(dummy > 0, 2, which))
# sum(!is_equal(sort(xmat_intervals), sort(onset_intervals)))
# sort(setdiff(xmat_intervals, onset_intervals))
# sort(setdiff(onset_intervals, xmat_intervals))
# which(censor_list[[subj_i]] < 1)
# length(xmat_intervals)
# length(onset_intervals)
# which(X[, "ON_TRIALS#0"] > 0) == (intervals + 1)
# (X[, "ON_TRIALS#10"] > 0) == (tr_grid > 0)
# m <- cbind(
#   # xmat = which(X[, "ON_TRIALS#0"] > 0),
#   # intervals + 1  ## plus 1 because of TENTzero
#   xmat = X[, "ON_TRIALS#1"] > 0,
#   own = tr_grid > 0
#   )
# image(m)
# plot(m)
# abline(0, 1)


## <<<<<< END CHECKING DUMMY MATRIX >>>>>>




## <<<<<< BEGIN CHECKING WITH PREVIOUS IMPLEMENTATION >>>>>>
# 
# A_list_ub55 <- readRDS(here("..", "ub55", "A_list.RDS"))
# subjs_common <- intersect(A_list_ub55$subj, subjs)
# A_list_ub55 <- A_list_ub55 %>% filter(subj %in% subjs_common)
# 
# subj_val = "448347"
# task_val = "Stroop"
# 
# named_vector <- function(type, nms) setNames(vector(type, length(nms)), nms)
# 
# r <- named_vector("logical", combo_paste(tasks, subjs_common))
# for (task_val in tasks) {
# 
#   for (subj_val in subjs_common) {
# 
#     a1 <- A_list_ub55 %>% filter(subj == subj_val, task == task_val) %>% select(mat)
#     a1_run1 <- a1[1, "mat"][[1]][[1]]
#     a1_run2 <- a1[2, "mat"][[1]][[1]]
# 
#     a1 <- rbind(
#       cbind(a1_run1, array(0, dim = dim(a1_run2))),
#       cbind(array(0, dim = dim(a1_run1)), a1_run2)
#     )
# 
#     nm <- paste0("wave1_", task_val, "_baseline_", subj_val)
#     a2 <- t(A_list[[nm]])
# 
#     r[[paste0(task_val, "_", subj_val)]] <- isTRUE(all.equal(c(a1), c(a2)))
#     # rowSums(is.na(a1)) > 0
#     # a2[rowSums(is.na(a1)) > 0, ]
# 
# 
#   }
# 
# }

## only diffs seem to be final trials at end of run, which sometimes are NA in original implementation.
## checking the censor lists, these trials/trs are not censored, so they should be in there.



## <<<<<< END CHECKING WITH PREVIOUS IMPLEMENTATION >>>>>>

