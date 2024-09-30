## subject IDs
subjs_wave12_all <- data.table::fread(here::here("in", "subjects_wave12_all.txt"), header = FALSE)$V1
subjs_wave13_all <- data.table::fread(here::here("in", "subjects_wave13_all.txt"), header = FALSE)$V1
subjs_wave23_all <- data.table::fread(here::here("in", "subjects_wave23_all.txt"), header = FALSE)$V1
subjs_wave12_good <- data.table::fread(here::here("in", "subjects_wave12_good.txt"), header = FALSE)$V1
## The `subjs_wave12_complete` in previous codes. 18 subjects.
#subjs_old_complete <- data.table::fread(here::here("in", "subjects_old_complete.txt"), header = FALSE)$V1
## DMCC5009144: DMCC2 (wave1) Stroop baseline missing
## 197449: Jo recommends not using Stroop baseline
## 178243: Stroop wave 1 baseline run 2 ended a bit early
## DMCC6418065 & DMCC6671683: missing Stroop behavioral data
subjs_wave12_all <- setdiff(subjs_wave12_all, c("DMCC5009144", "197449", "178243"))
subjs_wave12_good <- setdiff(subjs_wave12_good, c("DMCC5009144", "197449", "178243"))
subjs_wave12_complete <- setdiff(subjs_wave12_all, c("DMCC6418065", "DMCC6671683"))
