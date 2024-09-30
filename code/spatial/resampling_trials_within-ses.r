# Resamples trials stratified by confounding factors.
#
# Author: Michael Freund
#
# 07/21/2022: Fix a bug that fails to exclude trials with NAs
#
# Updated by Ruiqi Chen on 07/20/2022 for compatibility with wave13 & wave23

library(here)
library(tidyr)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(ggplot2)

source(here("code", "_constants.R"))
source(here("code", "_paths.R"))
source(here("code", "_subjects.R"))
source(here("code", "timeseries", "_utils_fmri.R"))

theme_set(theme_minimal(base_size = 14))

## input vars ----

sessions <- "baseline"
path_base <- "/data/nil-external/ccp/chenr/trr/"
variable <- "hilo_all"
classes <- c("lo", "hi")  ## -, +
task <- "Stroop"
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
subjs <- switch(toString(do_waves),
  "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all
)
input_fname <- file.path(path_base, "in", "behav",
  paste0("behavior-and-events_wave", do_waves[1], do_waves[2], "_", task, ".csv")
)
n_cores <- 4
n_resamples <- 100

file_name <- here("out", "spatial",
  paste0("trialidx_stroop_congruency_winses_wave", do_waves[1], do_waves[2], ".RDS")
)

## execute ----

print_tables <- TRUE
waves <- waves[do_waves]
tasks <- task  # In case of error
n_classes <- length(classes)

## read trial-wise coefficients:
alltrials <- read_results(
  waves = waves, tasks = task, sessions = sessions, subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) "errts_trials_target_epoch.RDS",
  read_fun = readRDS,
  n_cores = n_cores,
  path_base = file.path(path_base, "out/timeseries")
)
bads_idx <- lapply(alltrials, function(x) which(is.na(rowSums(x))))  ## get list of bad trials
bads <- data.table(
  trial = stack(bads_idx)$value,
  id = stack(bads_idx)$ind
)  # unlist() will rename duplicate ids (when there're multiple bad trials in a session)
#bads <- bads %>% tidyr::separate(id, c("wave", "task", "session", "subj"), sep = "_")

if (task != "Stroop" | variable != "hilo_all") stop("This script is only for Stroop!")

behav <- fread(input_fname, na.strings = c("", "NA"))
cols <- c("subj", "wave", "session", "run", "trial.num", "trialtype", "pc", "item", "color", "word")
behav <- behav[wave %in% do_waves, ..cols]
behav[, task := "Stroop"]
behav[, wave := paste0("wave", wave)]
behav[session == "bas", session := "baseline"]
behav[session == "pro", session := "proactive"]
behav[session == "rea", session := "reactive"]
for (ses in sessions) {
  behav[session == ses, trial := trial.num + n_trialspr[paste0("Stroop_", ses)] * (run - 1)]
}

## merge info on bad trials:
behav$is_good <- match(
  behav[, paste0(wave, "_", task, "_", session, "_", subj, "_", trial)],
  paste0(bads$id, "_", bads$trial)
) %>% is.na

## resample:

resample_to <- c(  ## per run*stimulus
  baseline = 3
  #proactive = 3,
  #reactive = 6
)
items_to_drop <- list(
  baseline_run1 = c("pinkGREEN", "blackYELLOW", "redBLUE", "whitePURPLE"),
  baseline_run2 = c("greenPINK", "yellowBLACK", "blueRED", "purpleWHITE")
)

resample_stroop <- function(nms, trials, n_resamples, resample_to) {
  trials_idx <- resample_idx(nms, n_resamples = n_resamples, resample_to = resample_to)
  ## convert index for trials into actual trial numbers.
  apply(trials_idx, 2, function(x) trials[x])
}

# Debugging
if (FALSE) {
  #subj_nm <- "448347"
  subj_nm <- subjs[1]
  wave_nm <- "wave1"
  session_nm <- "baseline"
}

res <- enlist(combo_paste(subjs, "__", waves, "__", sessions, "__", 1:2))
for (subj_nm in subjs) {
  for (wave_nm in waves) {
    for (session_nm in sessions) {
      for (run_nm in 1:2) {

        ## get relevant subset of trials, drop bad ones
        behav_i <- behav[subj == subj_nm & wave == wave_nm & session == session_nm & run == run_nm & is_good]

        ## drop certain items to uncorrelate confounding factors
        if (session_nm == "reactive") {
          behav_i <- behav_i[!grepl("RED|BLUE", item)]
        } else if (session_nm %in% c("baseline")) {
          drop_these <- items_to_drop[[paste0(session_nm, "_run", run_nm)]]
          behav_i <- behav_i[!item %in% drop_these]
        }
        conditions <- behav_i[, paste0(item, "_", trialtype)]
        global_min <- min(table(conditions))
        if (global_min < resample_to[session_nm]) {
          resamp2 <- global_min
        } else {
          resamp2 <- resample_to[session_nm]
        }

        trials <-
          resample_stroop(
            nms = conditions,
            trials = behav_i[, trial],
            n_resamples = n_resamples,
            resample_to = resamp2
          )
        #colnames(trials_resamp) <- gsub("(.*__)(.*)", "\\2", colnames(trials_resamp))

        nm <- paste0(subj_nm, "__", wave_nm, "__", session_nm, "__", run_nm)
        res[[nm]] <- trials

      }
    }
  }
}

## check that all counts are balanced

all_balanced <- vapply(
  res,
  function(m) {
    words <- gsub("[a-z]|_PC50|_bias|Con|InCon", "", colnames(m))
    colors <- gsub("[A-Z]|_PC50|_bias|Con|InCon", "", colnames(m))
    congruencies <- gsub("(^.*_bias|^.*_PC50)(.*)", "\\2", colnames(m))
    words_ok <- all(as.matrix(table(words, congruencies)) %*% cbind(c(-1, 1)) == 0)
    colors_ok <- all(as.matrix(table(colors, congruencies)) %*% cbind(c(-1, 1)) == 0)
    words_ok & colors_ok
  },
  logical(1)
) %>% all
stopifnot(all_balanced)

## save

saveRDS(res, file_name)


if (print_tables) {

  # lapply(
  #   res,
  #   function(m) {
  #     words <- gsub("[a-z]|_PC50|_bias|Con|InCon", "", colnames(m))
  #     congruencies <- gsub("(^.*_bias|^.*_PC50)(.*)", "\\2", colnames(m))
  #     table(words, congruencies)
  #   }
  # ) %>% print
  # lapply(
  #   res,
  #   function(m) {
  #     colors <- gsub("[A-Z]|_PC50|_bias|Con|InCon", "", colnames(m))
  #     congruencies <- gsub("(^.*_bias|^.*_PC50)(.*)", "\\2", colnames(m))
  #     table(colors, congruencies)
  #   }
  # ) %>% print
  counts <- sapply(res, ncol)
  print(counts)
  print(which(counts < 48))

}
