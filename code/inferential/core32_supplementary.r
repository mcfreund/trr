## NB: don't source this script, but run in interactive! weirdness wiht data.table referencing/envs...

library(here)
library(rmarkdown)
library(mfutils)

render_report(
  name = "core32_supplementary",
  src_dir = "",
  base_dir = here("code", "inferential"),
  output_dir = here("reports"),
  envir = parent.frame()   ## avoids data.table weirdness with accessing cols
)
