## plotting

library(ggplot2)
library(colorspace)
library(viridis)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(mfutils)

source(here::here("code", "_paths.R"))

## atlas

rois <- mfutils::schaefer2018_17_400_fsaverage5$key[["parcel"]]
atlas <- ggsegSchaefer::schaefer17_400
atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
n_bins_rois <- round(1 + log2(length(rois)))


## ggplot themes

theme_default <- function(
  base_size = 10, base_family = "",
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
      panel.grid = element_blank(),
      plot.subtitle = element_text(size = rel(1), hjust = 0),
      plot.title = element_text(size = rel(1), hjust = 0)
    ) +
    theme(...)
}


theme_surface <- function(
  base_size = 8, base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22,
  ...
) {
  theme_minimal(base_size = base_size, base_family, base_line_size, base_rect_size) %+replace%
    theme(
      axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
      axis.ticks = element_blank(),
      legend.position = c(0.5, 0.5),
      legend.title = element_text(size = 7),
      legend.background = element_blank(),
      legend.text = element_text(size = 7),
      legend.direction = "horizontal",
      legend.key.height = unit(1 / 8, "cm"),
      legend.key.width = unit(1 / 3, "cm"),
      plot.margin = unit(c(0, 0, 0, 0), "mm")
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
## examples for posterior densities:
example_rois <- c(
  "17Networks_LH_ContA_PFCl_2",
  "17Networks_LH_ContA_IPS_4"
)
example_roi_names <- c(
  `17Networks_LH_ContA_PFCl_2` = "left PFCl_2",
  `17Networks_LH_ContA_IPS_4` = "left IPS_4"
)
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


## color scales

colors_network8 <- setNames(legend_network8$color, legend_network8$network)
colors_response <- setNames(diverging_hcl(10, palette = "Purple-Brown")[c(2, 9)], c("rda", "uv"))
colors_models2 <- c(summarystat = "grey40", no_lscov_symm = "firebrick")
colors_models <- setNames(
  sequential_hcl(5, palette = "Purple-Blue"),
  c("full", "no_lscov", "no_lscov_symm", "fixed_sigma", "summarystat"))
colors_models_comparison <- c(no_lscov_symm = "#D306CD", summarystat = "#21D1FE")
colors_roi <- c("TRUE" = "firebrick", "FALSE" = "grey60")
colors_nois <- c(
  ContA = "firebrick", ContB = "firebrick", DorsAttnA = "firebrick",
  DorsAttnB = "firebrick", other = "grey60"
)
color_line_roi <- c("TRUE" = "white", "FALSE" = "black")


## sizes
axis_text_x_angle <- 30
size_line_roi <- c("TRUE" = 0.3, "FALSE" = 0.2)
two_column_width <- 183  ## mm
one_column_width <- 89  ## mm
oneandhalf_column_width <- 120 ## mm


## components of main figures

margins <- c(6, 6, 6, 6)
univ_stat_lab <- bquote(t^"+")
responses <- c("uv", "rda")

titles <- c(
  icc_hbm_uv_trr = bquote(
    bold("test-retest reliability") ~ "estimates in" ~ bold("univariate") ~ "Stroop contrasts"
  ),
  icc_hbm_mv_trr = bquote(
    bold("test-retest reliability") ~ "estimates in" ~ bold("multivariate") ~ "Stroop contrasts"
  ),
  uv_mv_trrprecision = bquote(
    bold("precision of") ~ "test-retest reliability estimates:"~ bold("univariate vs. multivariate")~"Stroop contrasts"
  )
)

comparison_factor_labs <- list(
  icc_hbm_trr = c(
    "summarystat" = "Summary statistic\n(intra-class corr.)",
    "no_lscov_symm" = "Hierarchical Bayes\n(MAP estimate)"
  ),
  icc_hbm_trr_scatter = c(
    "summarystat" = "Summary statistic",
    "no_lscov_symm" = "Hierarchical Bayes"
  ),
  icc_hbm_trr_short = c(
    "summarystat" = "sum. stat.",
    "no_lscov_symm" = "HBM"
  ),
  uv_mv = c("uv" = "Univariate\nStroop contrast", "rda" = "Multivariate\nStroop contrast"),
  uv_mv_short = c("uv" = "univar.", "rda" = "multivar.")
)


params_comparison_plots <- list(

  icc_hbm_uv = list(
    i_expr = quote(response == "uv" & sum_fun %in% c("pointest", "map")),
    contrast_expr = "atanh(no_lscov_symm) - atanh(summarystat)",
    comparison_factor = "model_nm",
    comparison_factor_order = comparison_factor_labs$icc_hbm_trr,
    title = titles$icc_hbm_uv,
    path = path_figs_results,
    filename = "icc_hbm_uv_trr",
    params_plot_surface = list(
      alpha_col = NULL,
      underlay = NULL,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      fill_label = "TRR (r)"
    ),
    params_plot_hist_means = list(
      text_label = comparison_factor_labs$icc_hbm_trr_short,
      colors = colors_models_comparison,
      x_lab = "TRR (r)"
    ),
    params_plot_hist_diff = list(
      colors = colors_roi,
      x_lab = "\U0394TRR: Bayes \U2212 Sum. Stat.\n(z difference)",
      text_y = c(50, 80),
      text_x = c(-0.85, -0.85),
      text_label = c("'ROIs'" = TRUE, "all parcels" = FALSE)
    ),
    params_plot_scatter = list(
      color_col = "tplus",
      axis_labs = comparison_factor_labs$icc_hbm_trr_scatter,
      legend_position = c(0.75, 0.25),
      limits_x = c(-0.5, 1)
    )
  ),

  icc_hbm_mv = list(
    i_expr = quote(response == "rda" & sum_fun %in% c("pointest", "map")),
    contrast_expr = "atanh(no_lscov_symm) - atanh(summarystat)",
    comparison_factor = "model_nm",
    comparison_factor_order = comparison_factor_labs$icc_hbm_trr,
    title = titles$icc_hbm_mv,
    path = path_figs_results,
    filename = "icc_hbm_mv_trr",
    params_plot_surface = list(
      alpha_col = NULL,
      underlay = NULL,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      fill_label = "TRR (r)"
    ),
    params_plot_hist_means = list(
      text_label = comparison_factor_labs$icc_hbm_trr_short,
      colors = colors_models_comparison,
      x_lab = "TRR (r)"
    ),
    params_plot_hist_diff = list(
      colors = colors_roi,
      x_lab = "\U0394TRR: Bayes \U2212 Sum. Stat.\n(z difference)",
      text_y = c(50, 95),
      text_x = c(-0.85, -0.85),
      text_label = c("'ROIs'" = TRUE, "all parcels" = FALSE)
    ),
    params_plot_scatter = list(
      color_col = "tplus",
      axis_labs = comparison_factor_labs$icc_hbm_trr_scatter,
      legend_position = c(0.75, 0.25),
      limits_x = c(-0.5, 1)
    )
  ),

  uv_mv_hbm_trrprecision = list(
    i_expr = quote(sum_fun == "sd"),
    j_expr = quote(value := 1 / log(value)),
    contrast_expr = "rda - uv",
    comparison_factor = "response",
    comparison_factor_order = comparison_factor_labs$uv_mv,
    title = titles$uv_mv_trrprecision,
    path = path_figs_results,
    filename = "uv_mv_hbm_trrprecision",
    params_plot_surface = list(
      alpha_col = NULL,
      underlay = NULL,
      #limits = c(-2, 1),
      #breaks = c(-1, 0, 1),
      fill_label = "Precision(TRR) =\n1/(log SD(TRR))"
    ),
    params_plot_hist_means = list(
      text_label = comparison_factor_labs$uv_mv_short,
      colors = colors_response,
      text_x = c(-1, -1),
      text_y = c(110, 75),
      x_lab = "Precision(TRR)"
    ),
    params_plot_hist_diff = list(
      colors = colors_roi,
      x_lab = "\U0394Precision(TRR):\nmultivar. \U2212 univar.",
      text_y = c(50, 95),
      text_x = c(-0.75, -0.75),
      text_label = c("'ROIs'" = TRUE, "all\nparcels" = FALSE)
    ),
    params_plot_scatter = list(
      color_col = "tplus", axis_labs = comparison_factor_labs$uv_mv,
      limits_x = NULL,
      limits_y = NULL,
      breaks_x = c(-1.3, -1, -0.7),
      breaks_y = c(-1.3, -1, -0.7),
      add_zero_lines = FALSE,
      legend_position = c(0.3, 0.95)
    )
  )

)
