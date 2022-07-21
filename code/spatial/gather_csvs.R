# gather_csvs.R - Combine data from different waves

# Author: Ruiqi Chen
#
# Recode all wave IDs to "wave1" and "wave2" for following analysis
# and gather all CSVs together.
#
# Note: it takes several minutes for IO. Should use `data.table` instead (20x faster).

library(tidyverse)
library(here)

in_path <- here("out", "spatial")
wave12_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave12.csv"
wave13_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave13.csv"
wave23_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess_wave23.csv"
output_file <- "projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess.csv"

wave12 <- read_csv(here(in_path, wave12_file))
wave13 <- read_csv(here(in_path, wave13_file)) %>%
  mutate(wave = recode(wave, "wave3" = "wave2"))
wave23 <- read_csv(here(in_path, wave23_file)) %>%
  mutate(wave = recode(wave, "wave2" = "wave1", "wave3" = "wave2"))

wave <- bind_rows(wave12, wave13, wave23)

unique(wave$wave)  # Should be "wave1" "wave2"

write.csv(wave, here(in_path, output_file))
