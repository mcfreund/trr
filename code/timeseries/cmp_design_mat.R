# Compare design matrices from condition-level and selective averging model

library(here)

source(here("code", "_constants.R"))
source(here("code", "_funs.R"))

# Constants
debugging <- T
name_task <- "Stroop"
name_session <- "reactive"
conditions <- c("PC50Con#0", "PC50InCon#0", "buffCon#0",
    "biasCon#0", "biasInCon#0")  # the first TR in every trial type
all_glm_dir <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/"

subjs <- subjs_wave12_all
do_waves <- c(1, 2)
waves <- waves[do_waves]
resid_type <- "errts"
name_task_session <- paste0(name_task, "_", name_session)
n_targets <- length(target_trs[[name_task]])

for (subj_i in seq_along(subjs)) {

    for (wave_i in seq_along(waves)) {

        # Load X
        name_glm <- paste0(name_session, "_", name_glms_dmcc[name_task])
        dir_glm <- file.path(all_glm_dir, wavedir_image[waves[wave_i]],
            "fMRIPrep_AFNI_ANALYSIS", subjs[subj_i], "1TRpK_RESULTS", name_task, name_glm)
        x <- mikeutils::read_xmat(file.path(dir_glm, "X.xmat.1D"))

        # Get the indices of trial onsets (use > 0.9 due to rounding problems)
        x_trial_start_idx <- sort(which(x[, conditions] > 0.9, arr.ind = T)[, 1])

        # Target TR indices
        if (debugging) debug_offset <- 0 else debug_offset <- 1  # debugging
        x_target_idx <- rep(x_trial_start_idx, each = n_targets) +
            target_trs[[name_task]] - debug_offset
        dim(x_target_idx) <- c(n_targets, length(x_trial_start_idx))

        # Load A
        name_glm <- paste0(name_session, "_", "null_2rpm_", waves[wave_i])
        dir_glm <- here("out", "timeseries", subjs[subj_i], "RESULTS", name_task, name_glm)
        a <- readRDS(here(dir_glm, paste0(resid_type, "_averaging_matrix.RDS")))

        for (i in seq_along(x_trial_start_idx)) {
            if (sum(a[x_target_idx[, i], ], na.rm = T) < 0.99) {
                print(paste0("Error found for subject: ", subj_i, ", wave: ", wave_i,
                    ",  trial (in X): ", i))
                print(paste0("Sum of A[target_TRs, ]: ",
                    sum(a[x_target_idx[, i], ], na.rm = T)))
                stop()
            }
        }

        # if (debugging) {
        #     n_trial_x <- length(x_trial_start_idx)
        #     n_trial_a <- sum(colSums(a) != 0)
        #     if (n_trial_x != n_trial_a) {
        #         print(paste0("Subject: ", subj_i, ", wave: ", wave_i,
        #             ", number of trials found in X: ", n_trial_x,
        #             ", number of trials found in A: ", n_trial_a))
        #     }
        # }

        # if (debugging) break

    }
    # if (debugging) break
}