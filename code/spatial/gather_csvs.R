# gather_csvs.R - Combine data from different waves

# Author: Ruiqi Chen
#
# Recode all wave IDs to "wave1" and "wave2" for following analysis
# and gather all CSVs together.

library(here)
library(data.table)

in_path <- here("out", "spatial")
wave12_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave12.csv"
wave13_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave13.csv"
wave23_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave23.csv"
output_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess.csv"

wave12 <- fread(here(in_path, wave12_file))
wave13 <- fread(here(in_path, wave13_file))[
  wave == "wave3", wave := "wave2"]
wave23 <- fread(here(in_path, wave23_file))[
  wave == "wave2", wave := "wave1"][
    wave == "wave3", wave := "wave2"]  ## Order matters!

waves <- rbind(wave12, wave13, wave23)

fwrite(waves, here(in_path, output_file))
