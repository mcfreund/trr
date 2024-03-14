## regions of interest

## from schaefer-400 atlas
core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345,
  349, 350, 351, 352, 354, 361, 365, 387
)
dmcc35 <- c(
  99, 127, 130, 132, 140, 141, 148, 337, 340, 349, 350, 361, 77, 78, 86, 87, 88, 90, 91, 101, 103, 105, 110, 139, 172,
  175, 185, 189, 290, 306, 311, 314, 346, 347, 353
)
atlas_nm <- "schaefer2018_17_400_fsaverage5"
rois <- mfutils::schaefer2018_17_400_fsaverage5$key[["parcel"]]
core32_nms <- rois[core32]  ## replace indices with names
dmcc35_nms <- rois[dmcc35]
#roi_col <- "parcel"  ## "parcel" or "network"

# if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
#   rois <- get(atlas_nm)$key[[roi_col]]
#   atlas <- schaefer17_400
#   atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
# } else {
#   stop("not configured for atlas")
# }