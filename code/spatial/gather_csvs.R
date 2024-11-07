# gather_csvs.R - Combine data from different waves

# Author: Ruiqi Chen
#
# Recode all wave IDs to "wave1" and "wave2" for following analysis
# and gather all CSVs together.

library(here)
library(data.table)
library(purrr)

gather_csvs <- function(
  output_type,
  filename_components = c("stroop__rda__n_resamples100__demean_run", "__cv_allsess"),
  in_path = here("out", "spatial"),
  do_network = FALSE,
  do_write = FALSE
  ) {
  
  stopifnot(output_type %in% c("projections", "noise_projs", "weights") && length(output_type) == 1)
  stopifnot(length(filename_components) == 2)
  network_str <- switch(do_network + 1, "", "__network")
  is_noiseprojs <- output_type == "noise_projs"
  file_type <- switch(is_noiseprojs + 1, "csv", "RDS")
  
  f_prefix <- paste0(output_type, "__", filename_components[1], network_str, filename_components[2])
  wave12_file <- paste0(f_prefix, "_wave12.", file_type)
  wave13_file <- paste0(f_prefix, "_wave13.", file_type)
  wave23_file <- paste0(f_prefix, "_wave23.", file_type)
  output_file <- paste0(f_prefix, ".", file_type)
  funcs <- switch(
    is_noiseprojs + 1,
    list(read_fun = fread, write_fun = fwrite),
    list(read_fun = readRDS, write_fun = saveRDS)
  )
  wave12 <- funcs$read_fun(here(in_path, wave12_file))
  wave13 <- funcs$read_fun(here(in_path, wave13_file))[wave == "wave3", wave := "wave2"]
  ## Order matters:
  wave23 <- funcs$read_fun(here(in_path, wave23_file))[wave == "wave2", wave := "wave1"]
  wave23 <- wave23[wave == "wave3", wave := "wave2"]
  waves <- rbind(wave12, wave13, wave23)
  if (do_write) funcs$write_fun(waves, here(in_path, output_file))
 
  waves

}

g <- expand.grid(
  output_type = c("projections", "noise_projs", "weights"),
  do_network = c(FALSE, TRUE)
)

pwalk(g, gather_csvs, do_write = TRUE)