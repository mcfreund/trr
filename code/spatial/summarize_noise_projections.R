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

d <- readRDS(here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess.RDS"))

x <- d[test_session == "baseline" & task == "Stroop"]
x[, proj_uv_scaled := abs(proj_uv_scaled)]
x[, proj_rda_scaled := abs(proj_rda_scaled)]

cols_values <- c("proj_rda_scaled", "proj_mv_relvar", "weights_cossim")
cols_by <- c("proj_uv_scaled", "proj_uv_relvar", "dimension", "var_total", "subj", "wave", "roi")
x_mean <- x[, lapply(.SD, mean), .SDcols = cols_values, by = cols_by]
x_resamplei <- x[resample_idx == 1]

saveRDS(x_mean, here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_means_baseline.RDS"))
saveRDS(x_resamplei, here("out", "spatial", "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_resample1_baseline.RDS"))


## ----

x[, proj_mv_relvar := abs(proj_mv_relvar)]
x[, proj_uv_relvar := abs(proj_uv_relvar)]
x[, proj_uv_absvar := abs(proj_uv_relvar)*var_total]
x[, proj_mv_absvar := abs(proj_mv_relvar)*var_total]
x[, var_dim_uv := proj_uv_absvar / proj_uv_scaled]  ## variance per dimension
x[, var_dim_mv := proj_mv_absvar / proj_rda_scaled]
stopifnot(all.equal(x$var_dim_uv, x$var_dim_mv))  ## check for equality
x <- x %>%
  rename(var_dim = var_dim_mv, region = roi) %>%
  mutate(var_dim_uv = NULL)


## get noise alignment of each weight vector per dimension per region
x_l <- melt(x,
    measure.vars = c("proj_uv_scaled", "proj_rda_scaled"),
    id.vars = c("subj", "region", "dimension", "resample_idx", "var_dim"))
x_sum <- x_l[, .(
    cosine_sim = tanh(mean(atanh(value))),
    sd_dim = mean(sqrt(var_dim))),
    by = .(subj, variable, region, dimension)]

## get total noise (average SD per region*subject)
x_total <- x[, .(
    proj_mv_absvar = mean(sqrt(proj_mv_absvar)),
    proj_uv_absvar = mean(sqrt(proj_uv_absvar))
    ),
    by = c("subj", "region", "wave", "resample_idx")]
x_total[, delta_absvar := log(proj_uv_absvar / proj_mv_absvar)]
x_total <- x_total[, .(
    delta_absvar = mean(delta_absvar),
    proj_mv_absvar = mean(proj_mv_absvar),
    proj_uv_absvar = mean(proj_uv_absvar)
    ),
    by = c("region", "subj")]

fwrite(x_total, here("out", "spatial", "noise_projs_mean_total.csv"))
fwrite(x_sum, here("out", "spatial", "noise_projs_mean_spectrum.csv"))
