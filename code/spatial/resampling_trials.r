# Resamples trials stratified by confounding factors.
#
# Author: Michael Freund

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
source(here("code", "_funs.R"))

theme_set(theme_minimal(base_size = 14))

## input vars ----

variable <- "hilo_all"
classes <- c("lo", "hi")  ## -, +
task <- "Stroop"
train <- c("proactive", "reactive")
test <- c("baseline")
subjs <- subjs_wave12_complete
glm_nm <- "null_2rpm"
resid_type <- "errts"
do_waves <- c(1, 2)
n_cores <- 12
n_resamples <- 100

file_name <- here("out", "spatial", "trialidx_stroop_congruency.RDS")

## execute ----

draw_plots <- TRUE
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
bads_idx <- lapply(alltrials, function(x) which(is.na(rowSums(x))))  ## get list of bad trials
bads <- data.table(
  trial = unlist(bads_idx),
  id = names(unlist(bads_idx))
)
#bads <- bads %>% tidyr::separate(id, c("wave", "task", "session", "subj"), sep = "_")

if (task == "Stroop" & variable == "hilo_all") {

behav <- fread(here::here("in", "behav", "behavior-and-events_wave12_Stroop.csv"), na.strings = c("", "NA"))
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

ttypes_to_resample <- c(  ## per run
  baseline = "biasCon",
  proactive = "biasInCon",
  reactive = "biasInCon"
)

resample_to <- c(  ## per run*stimulus
  baseline = 4, ## actually 4.5
  proactive = 1, ## actually 1.5
  reactive = 6
)

resample_stroop <- function(nms, trials, n_resamples, resample_to) {
  trials_idx <- resample_idx(nms, n_resamples = n_resamples, resample_to = resample_to)
  ## convert index for trials into actual trial numbers.
  apply(trials_idx, 2, function(x) trials[x])
}

res <- enlist(combo_paste(subjs, "__", waves, "__", sessions))
for (subj_nm in subjs) {
  for (wave_nm in waves) {
    for (session_nm in sessions) {

      ## get relevant subset of trials, drop bad ones
      behav_i <- behav[subj == subj_nm & wave == wave_nm & session == session_nm & is_good]
      #behav_i <- behav[subj == subj_nm & wave == wave_nm & session == session_nm]

      ## drop "buffcon" in reactive
      if (session_nm == "reactive") behav_i <- behav_i[!grepl("RED|BLUE", item)]

      ## trials to resample (bias)
      behav_i_bias_resamp <- behav_i[trialtype == ttypes_to_resample[session_nm], ]
      global_min <- min(table(behav_i_bias_resamp[, paste0(item, "_", run, "__", trialtype)]))
      if (global_min < resample_to[session_nm]) {
        resamp2 <- global_min
      } else {
        resamp2 <- resample_to[session_nm]
      }
      trials_bias_resamp <-
        resample_stroop(
          nms = behav_i_bias_resamp[, paste0(item, "_", run, "__", trialtype)],
          trials = behav_i_bias_resamp[, trial],
          n_resamples = n_resamples,
          resample_to = resamp2
        )

      if (session_nm %in% c("baseline", "proactive")) {

        ## because baseline session requires resampling the trials (over both runs) of each stimulus to an odd nubmer
        ## (9), fpr each stimulus, randomly resample one run to 4 and other to 5.
        n_extra_cols <- switch(session_nm, baseline = 4, proactive = 12)
        names_extra_cols <- rep(ttypes_to_resample[session_nm], n_extra_cols)
        stimuli_extra <- switch(
          session_nm,
          baseline = c("blueBLUE", "redRED", "purplePURPLE", "whiteWHITE"),
          proactive =
            c("redBLUE", "purpleBLUE", "whiteBLUE", "blueRED", "purpleRED", "whiteRED", "bluePURPLE", "redPURPLE",
              "whitePURPLE", "blueWHITE", "redWHITE", "purpleWHITE"
            )
          )

        extra_trials <- matrix(NA, nrow = n_resamples, ncol = n_extra_cols)
        colnames(extra_trials) <- rep(ttypes_to_resample[session_nm], n_extra_cols)
        for (i in seq_len(n_resamples)) {

          extra_in_run1 <- sample(stimuli_extra, length(stimuli_extra) / 2)
          extra_in_run2 <- setdiff(stimuli_extra, extra_in_run1)
          trials_run1 <- resample_stroop(
            behav_i[item %in% extra_in_run1 & run == 1 & !trial %in% trials_bias_resamp[i, ], item],
            behav_i[item %in% extra_in_run1 & run == 1 & !trial %in% trials_bias_resamp[i, ], trial],
            1,
            1
          )
          trials_run2 <- resample_stroop(
            behav_i[item %in% extra_in_run2 & run == 2 & !trial %in% trials_bias_resamp[i, ], item],
            behav_i[item %in% extra_in_run2 & run == 2 & !trial %in% trials_bias_resamp[i, ], trial],
            1,
            1
          )
          extra_trials[i, ] <- c(trials_run1, trials_run2)
        }

        trials_bias_resamp <- cbind(trials_bias_resamp, extra_trials)

      }
      colnames(trials_bias_resamp) <- gsub("(.*__)(.*)", "\\2", colnames(trials_bias_resamp))

      ## PC50 (don't resample)
      behav_i_pc50 <- behav_i[pc == "PC50", ]
      trials_pc50 <-
        resample_stroop(
          nms = behav_i_pc50[, trialtype],
          trials = behav_i_pc50[, trial],
          n_resamples = n_resamples,
          resample_to = NULL
        )

      ## bias trials to NOT resample
      behav_i_bias <- behav_i[trialtype != ttypes_to_resample[session_nm] & pc == "bias", ]
      trials_bias <-
        resample_stroop(
          nms = behav_i_bias[, trialtype],
          trials = behav_i_bias[, trial],
          n_resamples = n_resamples,
          resample_to = NULL
        )

      dim(trials_bias_resamp)
      dim(trials_bias)
      dim(trials_pc50)
      #stopifnot(length(unique(rowSums(trials_pc50))) == 1)  ## check that all resamples have same set of trials
      #stopifnot(length(unique(rowSums(trials_bias))) == 1)  ## check that all resamples have same set of trials

      ## bind all together
      trials <- cbind(trials_bias_resamp, trials_pc50, trials_bias)
      colnames(trials) <- gsub("bias|PC50", "", colnames(trials))

      nm <- paste0(subj_nm, "__", wave_nm, "__", session_nm)
      res[[nm]] <- trials

    }
  }
}

## save

saveRDS(res, file_name)


if (draw_plots) {

  ## check that classes are balanced

  counts <- enlist(names(res))
  for (ii in seq_along(res)) {

    trials <- res[[ii]]
    nms <- strsplit(names(res)[ii], "__", c("subj", "wave", "session"))[[1]]
    behav_i <- behav[subj == nms[1] & wave == nms[2] & session == nms[3], c("color", "word", "trial", "run", "pc")]
    behav_i[, congruency := ifelse(color == tolower(word), "congr", "incon")]

    counts[[ii]] <- lapply(
      seq_len(nrow(trials)),
      function(resample_i) {
        dt <- behav_i[trial %in% trials[resample_i, ], c("color", "word", "run", "pc", "congruency")]
        as.data.table(table(dt))
      }
    ) %>%
    rbindlist(idcol = "iter")

  }
  allcounts <- rbindlist(counts, idcol = "id")
  allcounts <- separate(allcounts, id, c("subj", "wave", "session"), "__")
  allcounts <- allcounts[N > 0]


  allcounts %>%
    ggplot(aes(color, N, fill = congruency)) +
    stat_summary(fun = sum, geom = "bar", position = position_dodge(width = 1/3), width = 0.25) +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_color.pdf"), height = 6, width = 6, dev = "pdf")

  allcounts %>%
    ggplot(aes(word, N, fill = congruency)) +
    stat_summary(fun = sum, geom = "bar", position = position_dodge(width = 1/3), width = 0.25) +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_word.pdf"), height = 6, width = 6, dev = "pdf")

  allcounts %>%
    ggplot(aes(run, N, fill = congruency)) +
    stat_summary(fun = sum, geom = "bar", position = position_dodge(width = 1/3), width = 0.25) +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_run.pdf"), height = 3, width = 6, dev = "pdf")


  ## boxplots:

  allcounts[, .(N = sum(N)), by = c("subj", "wave", "session", "iter", "congruency", "color")] %>%
    ggplot(aes(color, N, fill = congruency)) +
    geom_boxplot() +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_color_box.pdf"), height = 6, width = 6, dev = "pdf")

  allcounts[, .(N = sum(N)), by = c("subj", "wave", "session", "iter", "congruency", "word")] %>%
    ggplot(aes(word, N, fill = congruency)) +
    geom_boxplot() +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_word_box.pdf"), height = 6, width = 6, dev = "pdf")

  allcounts[, .(N = sum(N)), by = c("subj", "wave", "session", "iter", "congruency", "run")] %>%
    ggplot(aes(run, N, fill = congruency)) +
    geom_boxplot() +
    facet_grid(cols = vars(wave, session)) +
    scale_fill_manual(values = c(congr = "steelblue", incon = "firebrick")) +
    coord_flip() +
    theme(legend.position = "none") +
    labs(caption = "congr = blue, incon = red")
  ggsave(here("out", "spatial", "class-balancing_run_box.pdf"), height = 3, width = 6, dev = "pdf")

}
