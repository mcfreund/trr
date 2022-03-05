library(colorout)
library(httpgd)
library(here)
library(reticulate)
library(tidyr)
library(dplyr)
library(magrittr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mikeutils)
library(progress)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(viridis)
library(purrr)

theme_set(theme_half_open())
hgd()


source(here("code", "_constants.R"))
source(here("code", "_atlases.R"))
source(here("..", "ub55", "code", "_funs.R"))  ## for read_betas_dmcc()
source(here("code", "_funs.R"))
source(here("code", "_read_behav_wave12.R"))



## input vars ----

atlas <- "schaefer"
subjs <- subjs_wave12
do_waves <- c(1, 2)
hi <- c(Axcpt = "BX", Cuedts = "InConInc", Stern = "LL5RN", Stroop = "biasInCon")
lo <- c(Axcpt = "BY", Cuedts = "ConInc", Stern = "LL5NN", Stroop = "biasCon")
waves <- waves[do_waves]


## wrangle behavioral data:

behav_wave12$Axcpt <- rename(behav_wave12$Axcpt, rt = target.rt)
cols <- c("wave", "run", "session", "subj", "trialtype", "rt", "acc", "trial.num")
behav_wave12 <- lapply(behav_wave12, function(x) x[, ..cols])
behav_wave12 <- rbindlist(behav_wave12, idcol = "task")
behav_wave12 <- behav_wave12[order(wave, task, session, subj, run, trial.num), ]
behav_wave12[, trialnum := 1:.N, by = c("wave", "task", "subj", "session")]
behav_wave12[, c("run", "trial.num") := NULL]
behav_wave12[, wave := paste0("wave", wave)]

behav_wave12[session == "bas"]$session <- "baseline"
behav_wave12[session == "pro"]$session <- "proactive"
behav_wave12[session == "rea"]$session <- "reactive"



## execute ----


fn <- here("out", "icc", "trial-level-recovery_hilo_core32_wave12.csv")


if (file.exists(fn)) {
  
  d <- fread(fn)
  
} else {
  
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
        
        
        ## read trial-level residuals, condition-level betas
        
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
        )
        
        
        ## extract core32 etc...
        
        is_core32 <- schaefer10k %in% core32
        resid32 <- lapply(resid, function(x) x[, is_core32])
        betas32 <- betas[is_core32, trialtype_val, target_trs[[task_val]], ]
        rm(resid, betas)
        gc()
        
        
        ## residuals: get hilo contrast
        
        l <- 
          behav_wave12[wave == wave_val & task == task_val & session == session_val] %>% 
          split(.$subj)
        
        A <- 
          lapply(
            l, 
            function(x) {
              A <- model.matrix(~ 0 + trialtype, x)
              A <- sweep(A, 2, colSums(A), "/")
              colnames(A) <- gsub("trialtype", "", colnames(A))
              A[, trialtype_val] %*% rbind(1, -1)  ## gives hilo contrast
            }
            )
        resid32 <- resid32[names(A)]  ## make sure have same order
        
        resid_sum <- map2(A, resid32, ~crossprod(.x, .y))
        resid_sum <- abind(resid_sum, along = 1)
        
        
        ## betas: get hilo contrast
        
        betas32 <- aperm(betas32, c(1, 2, 4, 3))  ## put TR on 'outside'
        
        b <- rowMeans(betas32, dims = 3)
        b <- b[, trialtype_val["hi"], ] - b[, trialtype_val["lo"], ]
        
        
        ## correlate:
        
        resid_sum <- t(resid_sum)
        
        r <- colSums(scale2unit(meancenter(resid_sum)) * scale2unit(meancenter(b)))  ## linear corr
        neg_rmse <- -sqrt(colMeans((resid_sum - b) * (resid_sum - b)))  ## negative root mean square
        
        nm <- paste0(wave_val, "_", task_val, "_", session_val)
        res[[nm]] <- data.table(subj = colnames(b), r = r, neg_rmse = neg_rmse)
        
        print(nm)
  
        
      }
      
    }
    
  }
  
  
  
  d <- rbindlist(res, idcol = "id")
  d <- separate(d, id, c("wave", "task", "session"))
  
  fwrite(d, here("out", "icc", "trial-level-recovery_hilo_core32_wave12.csv"))
  
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
ggsave(here("out", "test_retest", "figs", "recovery_box_cor_hilo_core32.pdf"), p_cor_box, dev = "pdf", width = 7, height = 4)



p_cor_box_subj <- d %>%
  ggplot(aes(r, subj)) +
  geom_boxplot(width = 0.2, fill = "grey40")

p_cor_box_subj

ggsave(
  here("out", "test_retest", "figs", "recovery_box_subj_cor_hilo_core32.pdf"), 
  p_cor_box_subj, 
  dev = "pdf", width = 4.5, height = 5
  )
