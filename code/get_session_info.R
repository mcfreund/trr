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
table(sessinfo$sex)

round(range(sessinfo$age_baseline))
median(sessinfo$age_baseline)
round(c(quantile(sessinfo$age_baseline, 0.25), quantile(sessinfo$age_baseline, 0.75)))


## get dates ----

wave1 <- fread(here("in", "dates_DMCC2.txt"))
wave2 <- fread(here("in", "dates_DMCC3.txt"))
wave3 <- fread(here("in", "dates_DMCC4.txt"))

rep1 <- rbind(
    wave1[sub.id %in% c(subjs_wave12_complete, subjs_wave13_all)],
    wave2[sub.id %in% subjs_wave23_all]
    ) %>%
    select(-V1) %>% 
    melt(id.vars = c("sub.id", "task.id"), variable.name = "run", value.name = "date") %>%
    mutate(date = as.Date(date, format = "%m-%d-%Y")) %>%
    group_by(sub.id) %>%
    summarize(
        start_rep1 = min(date, na.rm = TRUE),
        end_rep1 = max(date, na.rm = TRUE)
    )
rep2 <- rbind(
    wave2[sub.id %in% subjs_wave12_complete],
    wave3[sub.id %in% c(subjs_wave13_all, subjs_wave23_all)]
    ) %>% 
    select(-V1) %>% 
    melt(id.vars = c("sub.id", "task.id"), variable.name = "run", value.name = "date") %>%
    mutate(date = as.Date(date, format = "%m-%d-%Y")) %>%
    group_by(sub.id) %>%
    summarize(
        start_rep2 = min(date, na.rm = TRUE),
        end_rep2 = max(date, na.rm = TRUE)
    )

data <- full_join(rep1, rep2, by = "sub.id")

time_btw <- sort(data$start_rep2 - data$end_rep1)
time_win1 <- sort(data$end_rep1 - data$start_rep1)
time_win2 <- sort(data$end_rep2 - data$start_rep2)

time_btw
time_win1
time_win2
f <- c(range, median, IQR)
lapply(f, function(x) x(time_btw))
lapply(f, function(x) x(time_win1))
lapply(f, function(x) x(time_win2))
lapply(f, function(x) x(c(time_win1, time_win2)))

