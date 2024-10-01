## regions of interest

## from schaefer-400 atlas
core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345,
  349, 350, 351, 352, 354, 361, 365, 387
)
## for schaefer-400 atlas, 17 network:
dmcc35 <- c(126, 69, 82, 83, 84, 130,  85, 105,  91, 106,  94,  99, 122, 124, 125, 129, 131, 132, 143, 186,
                  187, 182, 180, 282, 310, 312, 299, 341, 309, 330, 331, 332, 334, 344, 350);
atlas_nm <- "schaefer2018_17_400_fsaverage5"
rois <- mfutils::schaefer2018_17_400_fsaverage5$key[["parcel"]]
core32_nms <- rois[core32]  ## replace indices with names
dmcc35_nms <- rois[dmcc35]
