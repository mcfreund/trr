library(here)
library(tidyr)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(progress)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(viridis)
library(purrr)
source(here("code", "_constants.R"))
source(here("code", "_funs.R"))

theme_set(theme_half_open())

## input vars ----

atlas_nm <- "schaefer2018_7_400_fsaverage5"
do_waves <- c(1, 2)
subjs <- subjs_wave12_all
hi <- c(Axcpt = "BX", Cuedts = "InConInc", Stern = "LL5RN", Stroop = "biasInCon")
lo <- c(Axcpt = "BY", Cuedts = "ConInc", Stern = "LL5NN", Stroop = "biasCon")
atlas <- get(atlas_nm)
is_core32 <- atlas$data %in% core32  ## indicates vertices that belong to "core32", a set of 32 ROIs (parcels)
waves <- waves[do_waves]

behav_wave12 <- fread(here("in", "behav", "behavior-and-events_wave12_alltasks.csv"))
behav_wave12 <- behav_wave12[subj %in% subjs]



## run ----

fn <- here("out", "icc", "trial-level-recovery_hilo_core32_wave12.csv")


if (file.exists(fn)) {
  
  d <- fread(fn)
  
} else {
  
  ## iterators for dev
  # wave_i = 2
  # task_i = 4
  # session_i = 3
  # task_val <- tasks[task_i]
  # wave_val <- waves[wave_i]
  # session_val <- sessions[session_i]

  
  res <- enlist(combo_paste(waves, tasks, sessions))
  for (wave_i in seq_along(waves)) {  
    for (task_i in seq_along(tasks)) {
      for (session_i in seq_along(sessions)) {
       
        task_val <- tasks[task_i]
        wave_val <- waves[wave_i]
        session_val <- sessions[session_i]
        trialtype_val <- c(hi = hi[[task_val]], lo = lo[[task_val]])
        
        
        ## read trial-level coefficients (from selective averaging) and condition-level coefficients (from GLMs)
        
        ## trial-level coefs from selective averaging:
        resid <- read_results(
          wave_val, task_val, session_val, subjs,
          glmname = "null_2rpm", 
          filename_fun = function(...) "errts_trials_target_epoch.RDS",
          read_fun = readRDS
        )
        names(resid) <- gsub(paste0(wave_val, "_", task_val, "_", session_val, "_"), "", names(resid))

        betas <- read_betas_dmcc(
          .subjs = subjs, 
          .task = task_val, 
          .glm = paste0(session_val, "_", name_glms_dmcc[task_val]), 
          .dir = file.path(
            "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS", wavedir_image[wave_val], "fMRIPrep_AFNI_ANALYSIS"
            )
        )  ## condition-level coefs from GLM
        
        
        ## extract core32 etc...
        resid32 <- lapply(resid, function(x) x[, is_core32])
        betas32 <- betas[is_core32, trialtype_val, target_trs[[task_val]], ]
        rm(resid, betas)  ## remove from global env
        gc()  ## free memory
        
        
        ## get high-demand minus low-demand contrast

        ## first for selective-averaged coefs:
        
        ## get vector A that averages over trials (by trialtype) and contrasts averages across trialtypes:
        A <- 
          behav_wave12[wave == wave_val & task == task_val & session == session_val] %>% 
          split(.$subj) %>%
          lapply(function(x) averaging_matrix(x$trialtype)[, trialtype_val] %*% rbind(1, -1))
        resid32 <- resid32[names(A)]  ## make sure have same order

        resid_sum <- 
          map2(
            A, resid32, 
            function(.x, .y) {
                is_censored_trial <- is.na(rowSums(.y))
                crossprod(.x[!is_censored_trial, ], .y[!is_censored_trial, ])
                }
          )  ## apply
        resid_sum <- abind(resid_sum, along = 1)
        ## resid_sum is a subject by vertex matrix. 
        ## each vertex indicates the average difference (over trials) in BOLD signal between high-demand and low-demand 
        ## conditions, within our frontoparietal regions of interest (core32)        
        
        ## now, for condition-level GLM coefficients:
        betas32 <- aperm(betas32, c(1, 2, 4, 3))  ## put TR on 'outside'        
        b <- rowMeans(betas32, dims = 3)  ## average over TRs
        b <- b[, trialtype_val["hi"], ] - b[, trialtype_val["lo"], ]
        
        ## reduce to common set of subjects
        resid_sum <- t(resid_sum)
        subjs_intersect <- intersect(colnames(resid_sum), colnames(b))
        resid_sum <- resid_sum[, subjs_intersect]
        b <- b[, subjs_intersect]
        
        ## correlate:
        
        r <- colSums(scale2unit(center(resid_sum)) * scale2unit(center(b)))  ## linear corr
        neg_rmse <- -sqrt(colMeans((resid_sum - b) * (resid_sum - b)))  ## negative root mean square
        
        nm <- paste0(wave_val, "_", task_val, "_", session_val)
        res[[nm]] <- data.table(subj = colnames(b), r = r, neg_rmse = neg_rmse)
        
        print(nm)
  
        
      } 
    }  
  }
  
  
  
  d <- rbindlist(res, idcol = "id")
  d <- separate(d, id, c("wave", "task", "session"))
  
  fwrite(d, here("out", "timeseries", "trial-level-recovery_hilo_core32_wave12.csv"))
  
}


## plot ----


p_cor_box <- d %>%
  
  mutate(wave = as.factor(ifelse(wave == "wave1", 1, 2))) %>%
  
  ggplot(aes(wave, r)) +
  geom_boxplot(width = 0.25, fill = "grey40") +
  facet_wrap(
    vars(task, session), nrow = 1, 
    labeller = labeller(
      session = c(baseline = "bas", proactive = "pro", reactive = "rea")
      )
    ) +
  labs(y = "linear correlation") +
  theme(strip.background = element_blank())

p_cor_box
ggsave(here("out", "figs", "timeseries-comparison_glm-vs-selav_hilo-correlation_core32.pdf"), dev = "pdf", width = 7, height = 4)



d %>%
  ggplot(aes(r, subj)) +
  geom_boxplot(width = 0.2, fill = "grey40")

ggsave(
  here("out", "figs", "timeseries-comparison_glm-vs-selav_hilo-correlation_core32_subjs.pdf"), 
  dev = "pdf", width = 4.5, height = 5
  )
