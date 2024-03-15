## paths
source(here::here("code", "_atlas.R"))

if (Sys.info()["nodename"] == "puter2") {

  path_out <- "/mnt/d/trr_data/"

} else if (Sys.info()["nodename"] == "CCPLINUX1") {

  path_out <- "/data/nil-external/ccp/chenr/trr/out"

}

path_reliab <- here("out", "inferential", atlas_nm)
path_figs <- here::here("figs")