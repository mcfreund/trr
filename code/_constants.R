## system info

n_core <- parallel::detectCores()

## image, design, analysis info

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



subjs_wave12_all <- data.table::fread(here::here("in", "subjects_wave12_all.txt"), header = FALSE)$V1
subjs_wave13_all <- data.table::fread(here::here("in", "subjects_wave13_all.txt"), header = FALSE)$V1
subjs_wave23_all <- data.table::fread(here::here("in", "subjects_wave23_all.txt"), header = FALSE)$V1
subjs_wave12_good <- data.table::fread(here::here("in", "subjects_wave12_good.txt"), header = FALSE)$V1
subjs_old_complete <- data.table::fread(here::here("in", "subjects_old_complete.txt"),
  header = FALSE)$V1  # The `subjs_wave12_complete` in previous codes. 18 subjects.

## DMCC5009144: DMCC2 (wave1) Stroop baseline missing
## 197449: Jo recommends not using Stroop baseline
## 178243: Stroop wave 1 baseline run 2 ended a bit early
## DMCC6418065 & DMCC6671683: missing Stroop behavioral data

subjs_wave12_all <- setdiff(subjs_wave12_all, c("DMCC5009144", "197449", "178243"))
subjs_wave12_good <- setdiff(subjs_wave12_good, c("DMCC5009144", "197449", "178243"))
subjs_wave12_complete <- setdiff(subjs_wave12_all, c("DMCC6418065", "DMCC6671683"))

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

## regions of interest

## from schaefer-400 atlas
core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345, 
  349, 350, 351, 352, 354, 361, 365, 387
)
dmcc35 <- c(
  99, 127, 130, 132, 140, 141, 148, 337, 340, 349, 350, 361, 77, 78, 86, 87, 88, 90, 91, 101, 103, 105, 110, 139, 172, 
  175, 185, 189, 290, 306, 311, 314, 346, 347, 353
)


## strings

models <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma", "summarystat")
models_hbm <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma")
responses <- c("rda", "uv")
models_hbm <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma")
models_hbm <- setNames(models_hbm, c("Full", "ILS", "ILS Sym", "Homog"))
models <- c(models_hbm, "Sum Stat" = "summarystat")
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
rois <- mfutils::schaefer2018_17_400_fsaverage5$key[[roi_col]]
core32_nms <- rois[core32]  ## replace indices with names
dmcc35_nms <- rois[dmcc35]
path_out <- switch(
    Sys.info()["nodename"],
    puter2 = "/mnt/d/trr_data/",
    CCPLINUX1 = "/data/nil-external/ccp/chenr/trr/out"
)


## plotting
legend_network8 <- data.frame(
  network = c("Cont", "Default", "DorsAttn", "SalVentAttn", "SomMot", "Temp", "Vis", "Limbic"),
  color = c(qualitative_hcl(7, palette = "Dark 3"), "grey50"),
  label = c(
    "Cont", "\nDefault", "\n\nDorsAttn", "\n\n\nSalVentAttn", "\n\n\n\nSomMot", "\n\n\n\n\nTemp", "\n\n\n\n\n\nVis",
    "\n\n\n\n\n\n\nLimbic"
    ),
    stringsAsFactors = FALSE
)
colors_network8 <- setNames(legend_network8$color, legend_network8$network)
colors_response <- setNames(diverging_hcl(10, palette = "Purple-Brown")[c(2, 9)], c("rda", "uv"))
#colors_response <- setNames(diverging_hcl(2, palette = "Blue-Red2"), responses)
colors_models2 <- c(summarystat = "grey40", no_lscov_symm = "firebrick")
colors_models <- setNames(sequential_hcl(5, palette = "Purple-Blue"), models)
colors_nois <- c(
  ContA = "firebrick", ContB = "firebrick", DorsAttnA = "firebrick",
  DorsAttnB = "firebrick", other = "grey60"
  )
example_rois <- c(
    "17Networks_LH_ContA_PFCl_2",
    "17Networks_LH_ContA_IPS_4"
  )
example_roi_names <- c(
  `17Networks_LH_ContA_PFCl_2` = "left PFCl_2",
  `17Networks_LH_ContA_IPS_4` = "left IPS_4"
  )

axis_text_x_angle <- 30