## plotting

library(ggplot2)
library(colorspace)
library(viridis)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(mfutils)


## ggplot themes

# theme_set(theme_minimal(base_size = 12))
# theme_update(
#   strip.background = element_rect(fill = "transparent", color = "transparent"),
#   axis.line.y.left = element_line(),
#   axis.line.x.bottom = element_line(),
#   axis.ticks = element_line(),
#   panel.grid = element_blank()
# )

theme_default <- function(
  base_size = 12, base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22,
  ...
) {
  theme_minimal(base_size = base_size, base_family, base_line_size, base_rect_size) %+replace%
    theme(
      strip.background = element_rect(fill = "transparent", color = "transparent"),
      axis.line.y.left = element_line(),
      axis.line.x.bottom = element_line(),
      axis.ticks = element_line(),
      panel.grid = element_blank()
    ) +
    theme(...)
}


theme_surface <- function(
  base_size = 12, base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22,
  ...
) {
  theme_minimal(base_size = base_size, base_family, base_line_size, base_rect_size) %+replace%
    theme(
      axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
      axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_text(size = 7),
      legend.background = element_blank(), legend.text = element_text(size = 7), legend.direction = "horizontal",
      legend.key.height = unit(1 / 8, "cm"),
      legend.key.width = unit(1 / 3, "cm")
    ) +
    theme(...)
}


theme_set(theme_default())

## labels

legend_network8 <- data.frame(
  network = c("Cont", "Default", "DorsAttn", "SalVentAttn", "SomMot", "Temp", "Vis", "Limbic"),
  color = c(qualitative_hcl(7, palette = "Dark 3"), "grey50"),
  label = c(
    "Cont", "\nDefault", "\n\nDorsAttn", "\n\n\nSalVentAttn", "\n\n\n\nSomMot", "\n\n\n\n\nTemp", "\n\n\n\n\n\nVis",
    "\n\n\n\n\n\n\nLimbic"
  ),
  stringsAsFactors = FALSE
)


## color scales

colors_network8 <- setNames(legend_network8$color, legend_network8$network)
colors_response <- setNames(diverging_hcl(10, palette = "Purple-Brown")[c(2, 9)], c("rda", "uv"))
#colors_response <- setNames(diverging_hcl(2, palette = "Blue-Red2"), responses)
colors_models2 <- c(summarystat = "grey40", no_lscov_symm = "firebrick")
colors_models <- setNames(
  sequential_hcl(5, palette = "Purple-Blue"),
  c("full", "no_lscov", "no_lscov_symm", "fixed_sigma", "summarystat"))
colors_roi <- c("TRUE" = "firebrick", "FALSE" = "grey60")
colors_nois <- c(
  ContA = "firebrick", ContB = "firebrick", DorsAttnA = "firebrick",
  DorsAttnB = "firebrick", other = "grey60"
)


## sizes
axis_text_x_angle <- 30
size_line_roi <- c("TRUE" = 0.275, "FALSE" = 0.2)
color_line_roi <- c("TRUE" = "white", "FALSE" = "black")

## examples

example_rois <- c(
  "17Networks_LH_ContA_PFCl_2",
  "17Networks_LH_ContA_IPS_4"
)
example_roi_names <- c(
  `17Networks_LH_ContA_PFCl_2` = "left PFCl_2",
  `17Networks_LH_ContA_IPS_4` = "left IPS_4"
)

## atlas
rois <- mfutils::schaefer2018_17_400_fsaverage5$key[["parcel"]]
atlas <- ggsegSchaefer::schaefer17_400
atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)


## posterior stats

stat_names <- c(
  Mean = "mean", Median = "median", MAP = "map",
  `Lower 95% CI` = "hdi95_lower",
  `Upper CI 95` = "hdi95_upper",
  `Lower CI 89` = "hdi89_lower",
  `Upper CI 89` = "hdi89_upper",
  SD = "sd",
  IQR = "iqr",
  MAD = "mad",
  `95%ile` = "q95",
  `10%ile` = "q10",
  `5%ile` = "q05",
  rhat = "rhat",
  `Bulk ESS` = "ess_bulk",
  `Tail ESS` = "ess_tail"
)
