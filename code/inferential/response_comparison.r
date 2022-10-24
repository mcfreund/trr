## NB: don't source this script, but run in interactive! weirdness wiht data.table referencing/envs...

library(here)
library(rmarkdown)
library(mfutils)

mv <- "rda"
model <- "no_lscov_symm"
session <- c("reactive", "proactive", "baseline")
pars <- data.frame(mv = mv, model = model, session = session, stringsAsFactors = FALSE)

for (i in seq_len(nrow(pars))) {

  p <- c(pars[i, , drop = FALSE])

  render_report(
    name = "response_comparison_400",
    src_dir = "",
    params = p,
    base_dir = here("code", "inferential"),
    output_dir = here("reports"),
    envir = parent.frame()   ## avoids data.table weirdness with accessing cols
    )

}
