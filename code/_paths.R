## paths
source(here::here("code", "_atlas.R"))

if (Sys.info()["nodename"] == "puter2") {

  path_out <- "/mnt/d/trr_data/"

} else if (tolower(Sys.info()["nodename"]) == "ccplinux1") {

  path_out <- "/data/nil-external/ccp/chenr/trr/out"

}

path_reliab <- here::here("out", "inferential", atlas_nm)

path_figs <- here::here("figs")
path_figs_results <- file.path(path_figs, "main_results")
path_figs_modcomp <- file.path(path_figs, "model_comparison")
path_figs_results_supp <- file.path(path_figs_results, "supp")
path_figs_modcomp_supp <- file.path(path_figs_modcomp, "supp")