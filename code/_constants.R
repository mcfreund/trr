## image and design information

sec_tr <- 1.2  ## seconds per tr
n_vert <- 20484  ## surface hcp mesh
n_trs <- c(
  Axcpt_baseline   = 1220,
  Axcpt_proactive  = 1220,
  Axcpt_reactive   = 1220,
  Cuedts_baseline  = 1300,
  Cuedts_proactive = 1300,
  Cuedts_reactive  = 1300,
  Stern_baseline   = 1200,
  Stern_proactive  = 1200,
  Stern_reactive   = 1200,
  Stroop_baseline  = 1080,
  Stroop_proactive = 1080,
  Stroop_reactive  = 1180
)  ## BOTH RUNS INCLUDED
n_trialspr <- c(
  Axcpt_baseline = 72,
  Axcpt_proactive = 72,
  Axcpt_reactive = 72,
  Cuedts_baseline = 54,
  Cuedts_proactive = 54,
  Cuedts_reactive = 54,
  Stern_baseline = 45,
  Stern_proactive = 45,
  Stern_reactive = 45,
  Stroop_baseline = 108,
  Stroop_proactive = 108,
  Stroop_reactive = 120
)  ## number of trials (events) per subj*run


tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
taskruns <- sort(mfutils::combo_paste(tasks, c("run1", "run2")))
sessions <- c("baseline", "proactive", "reactive")
waves <- c("wave1", "wave2", "wave3")
wavedir_image <- c(wave1 = "HCP_SUBJECTS_BACKUPS", wave2 = "DMCC_Phase3", wave3 = "DMCC_Phase4")
wavedir_evts <- c(wave1 = "DMCC2", wave2 = "DMCC3", wave3 = "DMCC4")

name_glms_dmcc <- c(
  Axcpt = "Cues_EVENTS_censored",
  Cuedts = "CongruencyIncentive_EVENTS_censored",
  Stern = "ListLength_EVENTS_censored",
  Stroop = "Congruency_EVENTS_censored"
)

hi <- c(Axcpt = "BX", Cuedts = "InConInc", Stern = "LL5RN", Stroop = "biasInCon")
lo <- c(Axcpt = "BY", Cuedts = "ConInc", Stern = "LL5NN", Stroop = "biasCon")

## TRs of interest

## from jo:
# Axcpt: Cues, BX high, BY low. 17, 8:10
# Cuedts: CongruencyIncentive, InConNoInc high, ConNoInc low. 19, 9:11
# Stern: ListLength, LL5RN high, LL5NN low. 21, 12:14.
# Stroop: Congruency, biasInCon high, biasCon low. 13, 3:5

target_trs <- list(
  Axcpt = 8:10,
  Cuedts = 9:11,
  Stern = 12:14,
  Stroop = 3:5
)
