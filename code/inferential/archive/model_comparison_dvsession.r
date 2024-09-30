## NB: don't source this script, but run in interactive! weirdness wiht data.table referencing/envs...

library(here)
library(rmarkdown)
library(mfutils)

responses <- c("uv", "rda")
pars <- data.frame(response = responses, stringsAsFactors = FALSE)

for (i in seq_len(nrow(pars))) {

  p <- c(pars[i, , drop = FALSE])

  render_report(
    name = "model_comparison_dvsession",
    src_dir = "",
    params = p,
    base_dir = here("code", "inferential"),
    output_dir = here("reports"),
    envir = parent.frame()   ## avoids data.table weirdness with accessing cols
    )

}
