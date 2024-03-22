library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(data.table)
library(colorspace)
library(patchwork)


## read data

d <- readRDS(here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12.RDS"))


## pick first resample and first ROI in core32 to characterize
## plot: cosine similarity per dimension
## plot: trial-level variance per dimension
## plot: total variance (sum over all dimensions)

## example roi

x <- d[test_session == "baseline" & task == "Stroop"]
x[, proj_uv_scaled := abs(proj_uv_scaled)]
x[, proj_rda_scaled := abs(proj_rda_scaled)]

cols_values <- c("proj_rda_scaled", "proj_mv_relvar", "weights_cossim")
cols_by <- c("proj_uv_scaled", "proj_uv_relvar", "dimension", "var_total", "subj", "wave", "roi")
x_mean <- x[, lapply(.SD, mean), .SDcols = cols_values, by = cols_by]
x_resamplei <- x[resample_idx == 1]

saveRDS(x_mean, here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12_means_baseline.RDS"))
saveRDS(x_resamplei, here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12_resample1_baseline.RDS"))

#rsync -avzP freundm@${ccplinux_ip}:/data/nil-external/ccp/freund/trr/out/spatial/noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12_resample1_baseline.RDS /media/mcf/WD_BLACK/wustl_ccplinux1_backup/trr/out/spatial/
#rsync -avzP freundm@${ccplinux_ip}:/data/nil-external/ccp/freund/trr/out/spatial/noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12_means_baseline.RDS /media/mcf/WD_BLACK/wustl_ccplinux1_backup/trr/out/spatial/