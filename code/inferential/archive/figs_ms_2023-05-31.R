## NB: don't source this script, but run in interactive! weirdness wiht data.table referencing/envs...

library(here)
library(rmarkdown)
library(mfutils)

roi_obj <- c("top_t")
#roi_obj <- c("core32", "dmcc35", "top_t", "top_z")
pars <- data.frame(roi_obj = roi_obj, stringsAsFactors = FALSE)

for (i in seq_len(nrow(pars))) {

  p <- c(pars[i, , drop = FALSE])

  render_report(
    name = "figs_ms_2023-05-31",
    src_dir = "",
    params = p,
    base_dir = here("code", "inferential", "viz"),
    output_dir = here("reports"),
    envir = parent.frame()   ## avoids data.table weirdness with accessing cols
    )

}

