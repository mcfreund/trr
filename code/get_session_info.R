library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(data.table)
library(brms)
library(posterior)
library(mfutils)

source(here("code", "inferential", "_parameters_reliability.R"))
source(here("code", "inferential", "_utils_reliability.R"))
source(here("code", "_atlas.R"))
source(here("code", "_constants.R"))
source(here("code", "_paths.R"))
source(here::here("code", "_subjects.R"))

## read subject list from brms models ----

responses <- c("uv", "rda")
model_info <- rbind(
  ## all models for core32 parcels:
  CJ(model_nm = models_hbm, response = responses, roi_nm = core32_nms, session = "baseline"),
  ## only winning model (no_lscov_symm) for rest of parcels:
  CJ(model_nm = "no_lscov_symm", response = responses, roi_nm = rois[!rois %in% core32_nms], session = "baseline")
)

plan(sequential, split = TRUE)
responses <- c("uv", "rda")
model_info <- rbind(
  ## all models for core32 parcels:
  CJ(model_nm = models_hbm, response = responses, roi_nm = core32_nms, session = "baseline"),
  ## only winning model (no_lscov_symm) for rest of parcels:
  CJ(model_nm = "no_lscov_symm", response = responses, roi_nm = rois[!rois %in% core32_nms], session = "baseline")
)

modlist <- future_pmap(
    model_info[1],
    read_summarize_hbm,
    base_path = file.path(path_out, "inferential"),
    atlas_nm = atlas_nm
)
subjs_models <- rownames(ranef(modlist[[1]])$subj)

## compare to subj lists
subj_list <- list(
    wave12_all = subjs_wave12_all,
    wave23_all = subjs_wave23_all,
    wave13_all = subjs_wave13_all,
    wave12_complete = subjs_wave12_complete,
    wave12_good = subjs_wave12_good
)
subjs_all <- Reduce(union, subj_list)

## from spatial modeling scripts:
do_waves <- c(1, 2)
subjs_spatial <- switch(toString(do_waves),
  "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all
)
length(subjs_all)
length(subjs_models)
all(subjs_all %in% subjs_models)
all(subj_list$wave12_all %in% subjs_models)
all(subj_list$wave23_all %in% subjs_models)
all(subj_list$wave13_all %in% subjs_models)
all(subj_list$wave12_complete %in% subjs_models)
all(subj_list$wave12_good %in% subjs_models)
all(subjs_spatial %in% subjs_models)

## read subject info from ndar files ---

filename_sessinfo <- "/mnt/c/Users/mfreu/Downloads/ndar_subject01_DMCC.csv"
filename_bdays <- "/mnt/c/Users/mfreu/Downloads/NDARlookup.csv"

sessinfo <- fread(filename_sessinfo, skip = 1)
bdays <- fread(filename_bdays)

bdays <- rename(bdays[sub.id %in% subjs_models][order(birthday), c("sub.id", "birthday")], subj = sub.id)
sessinfo <- sessinfo[src_subject_id %in% subjs_models, c("src_subject_id", "visit", "interview_date", "sex")]

bdays[, birthday := as.Date(birthday, format = "%m/%d/%Y")]
sessinfo <- rename(sessinfo, subj = src_subject_id, date_baseline = interview_date)
sessinfo <- sessinfo[bdays, on = "subj"]

sessinfo[, date_baseline := as.Date(date_baseline, format = "%m/%d/%Y")]
sessinfo <- dcast(sessinfo, ... ~ visit, value.var = "date_baseline")

stopifnot(nrow(sessinfo) == length(subjs_models))

sessinfo[, age_baseline := as.numeric(difftime(wave1, birthday, units = "days") / 365.25)]
sessinfo[, time_12 := as.numeric(difftime(wave2, wave1, units = "days"))]
sessinfo[, time_23 := as.numeric(difftime(wave3, wave2, units = "days"))]

range(sessinfo$time_12)
median(sessinfo$time_12)
round(c(quantile(sessinfo$time_12, 0.25), quantile(sessinfo$time_12, 0.75)))
sort(sessinfo$time_23[complete.cases(sessinfo$time_23)])

table(sessinfo$sex)

round(range(sessinfo$age_baseline))
median(sessinfo$age_baseline)
round(c(quantile(sessinfo$age_baseline, 0.25), quantile(sessinfo$age_baseline, 0.75)))
sort(sessinfo$time_23[complete.cases(sessinfo$time_23)])
