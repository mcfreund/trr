library(here)
library(data.table)
library(doParallel)
library(foreach)
library(mfutils)

source(here("code", "_constants.R"))
source(here("code", "_paths.R"))
source(here("code", "_subjects.R"))
source(here("code", "timeseries", "_utils_fmri.R"))


## input vars ----

atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
subjs <- subjs_wave12_all
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
do_tasks <- c("Stroop")
n_cores <- 20

if (FALSE) {  ## for dev
  subj_i <- 1
  wave_i <- 1
  task_i <- 1
  session_i <- 1
}

## execute ----

atlas <- get(atlas_nm)
waves <- waves[do_waves]
tasks <- do_tasks

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <-
  foreach(subj_i = seq_along(subjs), .inorder = FALSE) %:%
  foreach(task_i = seq_along(tasks), .inorder = FALSE) %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE) %:%
  foreach(session_i = seq_along(sessions), .inorder = FALSE) %dopar% {

    name_task_i <- tasks[task_i]
    name_wave_i <- waves[wave_i]
    name_session_i <- sessions[session_i]
    name_subj_i <- subjs[subj_i]

    dir_glm <-
      here("out", "timeseries", name_subj_i, "RESULTS", name_task_i,
        paste0(name_session_i, "_", glm_nm, "_", name_wave_i)
        )
    fname_trials <- here(dir_glm, paste0(resid_type, "_trials_target_epoch.RDS"))
    trials <- t(readRDS(fname_trials))  ## vertex by trial
    ## create list of vertex by trial matrices (one matrix per roi):
    parcels <- parcellate(trials, atlas, col_roi = roi_col)
    means <- do.call(rbind, lapply(parcels, colMeans))  ## average across vertices per parcel, bind into single matrix
    colnames(means) <- sprintf("trial_%03d", seq_len(ncol(means)))
    means <- as.data.table(means)
    filename <- here(dir_glm, paste0("trial-means", "_", atlas_nm, "-", roi_col, "_resid-", resid_type, ".csv"))
    fwrite(means, filename)

}
stopCluster(cl)
