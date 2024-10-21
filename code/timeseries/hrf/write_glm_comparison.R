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
library(rhdf5)
library(gifti)
source(here("code", "_constants.R"))
source(here("code", "_subjects.R"))
source(here("code", "_atlas.R"))
source(here("code", "_paths.R"))
source(here("code", "timeseries", "_utils_fmri.R"))

theme_set(theme_half_open())

## functions ----

construct_filename_gifti <- function(
  subject, wave, session, run, glmname, hemi,
  task = "Stroop",
  base_dir = file.path("/data/nil-external/ccp/freund/stroop-rsa-pc/", "out", "glms"),
  prefix = "STATS",
  suffix = "_REML.func.gii"
){
  
  arg <- as.list(environment())
  if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

  file.path(
    base_dir, subject, wave, "RESULTS", task, paste0(session, "_", glmname),
    paste0(prefix, "_", subject, "_", run, "_", hemi, suffix)
  )
}

construct_filenames_gifti <- function(
    subjects, waves, sessions, runs, glmnames,
    hemis = c("L", "R"), 
    task = "Stroop", 
    base_dir = file.path("/data/nil-external/ccp/freund/stroop-rsa-pc/", "out", "glms"),
    prefix = "STATS", 
    suffix = "_REML.func.gii",
    returnDT = TRUE
    ){
    
    d <- expand.grid(
        subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmnames,
        hemi = hemis,
        stringsAsFactors = FALSE
    )

    missing_col <- names(d)[!names(d) %in% c("subject", "wave", "session", "run", "glmname", "hemi")]
    if (length(missing_col) > 0) stop("missing column:", missing_col)
    
    filename <- vector("character", nrow(d))
    for (row_i in seq_len(nrow(d))) {
        x <- d[row_i, ]
        filename[row_i] <- construct_filename_gifti(
            subject = x$subject, wave = x$wave, 
            session = x$session, run = x$run, 
            glmname = x$glmname, hemi = x$hemi,
            task = task, base_dir = base_dir, prefix = prefix, suffix = suffix
        )
    }

    if (returnDT) {
        d$filename <- filename
        data.table::setDT(d, key = c("subject", "wave", "session", "run", "hemi"))
        d
    } else {
        filename
    }


}


cor_diag <- function(x, y) {
    
  x <- scale(x, scale = FALSE)
  x <- sweep(x, 2, sqrt(colSums(x^2)), "/")

  y <- scale(y, scale = FALSE)
  y <- sweep(y, 2, sqrt(colSums(y^2)), "/")

  colSums(x * y)

}


## input vars ----

tasks <- "Stroop"
n_cores <- 18
atlas_nm <- "schaefer2018_7_400_fsaverage5"
do_waves <- c(1, 2)
waves <- waves[do_waves]
glmname <- "lsall_1rpm"
#glmname <- "lssep_1rpm"

if (glmname == "lssep_1rpm") {
    suffix <- ".gii"
    pattern <- NULL
} else if (glmname %in% c("lsall_1rpm", "condition_1rpm", "lsall_acompcor06_1rpm")) {
    suffix <- "_REML.func.gii"
    pattern <- "_Coef"
} else {
    stop("not configured for provided glmname")
}
subjs <- switch(toString(do_waves), "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all)
hi <- c(Axcpt = "BX", Cuedts = "InConInc", Stern = "LL5RN", Stroop = "biasInCon")
lo <- c(Axcpt = "BY", Cuedts = "ConInc", Stern = "LL5NN", Stroop = "biasCon")
atlas <- get(atlas_nm)
is_core32 <- atlas$data %in% core32  ## indicates vertices that belong to "core32", a set of 32 ROIs (parcels)
b <- data.table::fread(file.path(
  "/data/nil-external/ccp/freund/stroop-rsa-pc/in",
  "behavior-and-events_stroop_2021-10-20_nice.csv")
)
b <- dplyr::arrange(b, "subj", "wave", "session", "run", "trial_num")


## read ----

## read selective-averaged coefficients:

alltrials <- read_results(
  waves = waves,
  tasks = tasks,
  sessions = sessions,
  subjs = subjs,
  glmname = "null_2rpm",
  filename_fun = function(...) "errts_trials_target_epoch.RDS",
  read_fun = readRDS,
  n_cores = n_cores,
  path_base = file.path(path_out, "timeseries")
)

## get info to read GLM-all coefficients:

input <- construct_filenames_gifti(
    subject = subjs, wave = waves, session = sessions, run = c("run1", "run2"),
    glmname = glmname,
    suffix = suffix
    )

## compute:

cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
allres <-
  foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "c") %:%
  foreach(sess_i = seq_along(sessions), .inorder = FALSE, .combine = "c") %:%
  foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %dopar% {

    ## get info for glm-all coefs:
    input_i <- input[subject == subjs[subj_i] & session == sessions[sess_i] & wave == waves[wave_i], ]
    if (any(!file.exists(input_i$filename))) return(NULL)
    
    ## get info for selective-averaged coefs:
    selav_id <- paste(waves[wave_i], "Stroop", sessions[sess_i], subjs[subj_i], sep = "_")
    
    ## get trial IDs and contrast matrix:
    congruency <- 
      b[subj == subjs[subj_i] & task == "Stroop" & wave == waves[wave_i] & session == sessions[sess_i],
      trial_type
      ]
    w <- mfutils::averaging_matrix(congruency)[, c("incon", "congr")] %*% rbind(1, -1)

    ## glm-all:
    ## read and bind both runs:
    input1 <- input_i[run == "run1", ]
    input2 <- input_i[run == "run2", ]
    giftis1 <- lapply(input1[, filename], read_gifti) %>% concat_hemis(input1$hemi, pattern = pattern)
    giftis2 <- lapply(input2[, filename], read_gifti) %>% concat_hemis(input2$hemi, pattern = pattern)
    giftis <- cbind(giftis1, giftis2)
    ## get Stroop contrast:
    giftis_stroop <- giftis %*% w
    ## parcellate:
    parcs_stroop_glmall <- parcellate(giftis_stroop, atlas, col_roi = "parcel")

    ## selective-averaged:
    selav_all <- t(alltrials[[selav_id]])
    ## get Stroop contrast:
    selav_stroop <- selav_all %*% w
    ## parcellate:
    parcs_stroop_selav <- parcellate(selav_stroop, atlas, col_roi = "parcel")
    
    ## compare and store:
    r_contr <- Map(function(x, y) cor(x, y), parcs_stroop_glmall, parcs_stroop_selav)
    x_mean <- Map(function(x, y) colMeans(cbind(glmall = x, selav = y)), parcs_stroop_glmall, parcs_stroop_selav)
    r_trial <- sapply(parcellate(as.matrix(cor_diag(t(selav_all), t(giftis))), atlas, col_roi = "parcel"), median, na.rm = TRUE)
    out <- data.table(
      r_spatial = c(do.call(rbind, r_contr)),
      r_temporal = r_trial,
      do.call(rbind, x_mean),
      atlas$key,
      subject = subjs[subj_i],
      wave = waves[wave_i],
      session = sessions[sess_i]
    )

    list(out)
 
  }
res <- rbindlist(allres)
fwrite(res, here("out", "timeseries", paste0("glm_comparison_", glmname, "-vs-selav.csv")))

