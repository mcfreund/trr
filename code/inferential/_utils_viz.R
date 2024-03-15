library(ggplot2)

arsinh <- scales::trans_new("arsinh", transform = function(x) asinh(x), inverse = function(x) sinh(x))


highlight_overlay <- function(
  x, upper,
  lower = min(x),
  upper_alpha = 1,
  lower_alpha = 0.1,
  scaling_fun = \(x) x^2
) {
  alphas <- numeric(length(x))
  is_upper <- x > upper
  is_lower <- x < lower
  is_mid <- !is_upper & !is_lower

  alphas[is_upper] <- upper_alpha
  alphas[is_lower] <- lower_alpha

  alphas[is_mid] <- lower_alpha + (upper_alpha - lower_alpha) * scaling_fun((x[is_mid] - lower) / (upper - lower))

  alphas

}



plot_surface <- function(
  data,
  statistic_col,
  border_color_col = "is_roi",
  border_color_values = color_line_roi,
  border_size_col = "is_roi",
  border_size_values = size_line_roi,
  alpha_col = "alpha",
  position = ggseg::position_brain(. ~ hemi + side),
  underlay = c(color = "grey", fill = "grey"),
  atlas_data = atlas,
  scale_fill = \(...) scale_fill_viridis_c(option = "magma", na.value = "white", ...),
  scale_fill_breaks = NULL,
  theme_ = theme_surface(),
  guides_ = guides(
    fill = guide_colorbar(title.position = "left", title.vjust = 0.9),
    alpha = "none", color = "none", size = "none"
  ),
  ...
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  # Global aesthetics

  p <- ggplot(data, aes(fill = .data[[statistic_col]]))
  
  # Check if border_color is provided

  if (!is.null(border_color_col)) {
    p <- p + aes(color = .data[[border_color_col]])
    if (!any(is.null(border_color_values))) {
      p <- p + scale_color_manual(values = border_color_values)
    }
  }
  if (!is.null(border_size_col)) {
    p <- p + aes(size = .data[[border_size_col]])
    if (!any(is.null(border_size_values))) {
      p <- p + scale_size_manual(values = border_size_values)
    }
  }

  ## underlay

  if (any(!is.null(underlay))) {
    p <- p +
      ggseg::geom_brain(color = underlay[["color"]], fill = underlay[["fill"]], atlas = atlas_data, position = position)
  }

  ## overlay

  if (!is.null(alpha_col)) {
    p <- p + ggseg::geom_brain(aes(alpha = .data[[alpha_col]]), atlas = atlas_data, position = position)
  } else {
    p <- p + ggseg::geom_brain(atlas = atlas_data, position = position)
  }

  ## scales, guides, themes
  if (is.null(scale_fill_breaks)) {
    p <- p + scale_fill()
  } else {
    p <- p + scale_fill(breaks = scale_fill_breaks)
  }
  p <- p + guides_ + theme_

  p

}


# plot_parcel_stats <- function(
#   data, statistic_col, color_col, color_values, id_vars,
#   n_bins = 10,
#   hist_means_axis_titles = c(x = "ICC(Stroop) (r)", y = "Number of parcels")
#   hist_diff_axis_titles = c(x = "ICC(Stroop) (r)", y = "Number of parcels")
# ) {

  
#   p_hist_means <- data %>%
#     ggplot(aes(.data[[statistic_col]], fill = .data[[color_col]])) +
#     geom_histogram(position = "identity", alpha = 0.5, color = "black", bins = n_bins) +
#     scale_fill_manual(values = color_values) +
#     labs(x = hist_means_axis_titles[["x"]], y = hist_means_axis_titles[["y"]]) +
#     theme(legend.position = "none") +
#     ## TODO: convert to geomtext:
#     annotate("text", y = 75, x = -0.3, label = "univar.", color = color_values[["uv"]]) +
#     annotate("text", y = 75, x = 0.75, label = "multiv.", color = color_values[["rda"]])

#   # p_icc_hist_diff <- data %>%
#   #   pivot_wider(
#   #     id_cols = c("region", "is_roi", "alpha", "dprime", "p_plus"),
#   #     names_from = "response", values_from = "r"
#   #   ) %>%
#   #   mutate(difference = atanh(rda - uv)) %>%
#   #   ggplot(aes(difference)) +
#   #   geom_histogram(position = "identity", fill = colors_roi[["FALSE"]], color = "black", bins = n_bins) +
#   #   geom_histogram(
#   #     data = . %>% filter(is_roi),
#   #     position = "identity", fill = colors_roi[["TRUE"]], color = "black", bins = n_bins) +
#   #   scale_x_continuous(breaks = c(-1, 0, 1)) +
#   #   annotate("text", y = 100, x = -0.8, label = "all", color = colors_roi[["FALSE"]]) +
#   #   annotate("text", y = 50, x = -0.8, label = "'ROI'", color = colors_roi[["TRUE"]]) +
#   #   labs(
#   #     x = "\u0394ICC(Stroop):\nmultiv. \u2212 univar. (z)",
#   #     y = "Number of parcels"
#   #   )



# }


create_icc_hist_means <- function(data, response_col, r_col, n_bins, colors_response, x_label, y_label) {
  p_icc_hist_means <- data %>%
    ggplot(aes(.data[[r_col]], fill = .data[[response_col]])) +
    geom_histogram(position = "identity", alpha = 0.5, color = "black", bins = n_bins) +
    scale_fill_manual(values = colors_response) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme(legend.position = "none") +
    annotate("text", y = 75, x = -0.3, label = "univar.", color = colors_response[["uv"]]) +
    annotate("text", y = 75, x = 0.75, label = "multiv.", color = colors_response[["rda"]])

  return(p_icc_hist_means)
}

create_icc_hist_diff <- function(data, difference_col, is_roi_col, n_bins, colors_roi, x_label, y_label) {
  p_icc_hist_diff <- data %>%
    ggplot(aes(.data[[difference_col]])) +
    geom_histogram(position = "identity", fill = colors_roi[["FALSE"]], color = "black", bins = n_bins) +
    geom_histogram(
      data = . %>% filter(.data[[is_roi_col]]),
      position = "identity", fill = colors_roi[["TRUE"]], color = "black", bins = n_bins
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    annotate("text", y = 100, x = -0.8, label = "all", color = colors_roi[["FALSE"]]) +
    annotate("text", y = 50, x = -0.8, label = "'ROI'", color = colors_roi[["TRUE"]]) +
    labs(
      x = x_label,
      y = y_label
    )

  return(p_icc_hist_diff)
}

create_icc_scatter <- function(data, uv_col, rda_col, dprime_col, highlight_overlay, x_label, y_label, color_label) {
  p_icc_scatter <- data %>%
    arrange(.data[[dprime_col]]) %>%
    ggplot(aes(.data[[uv_col]], .data[[rda_col]])) +
    geom_abline() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_point(
      aes(
        color = .data[[dprime_col]],
        alpha = highlight_overlay(.data[[dprime_col]], upper = quantile(.data[[dprime_col]], 0.9), lower = min(.data[[dprime_col]]), lower_alpha = 0.5)
      ),
      stroke = 0, size = 1
    ) +
    scale_alpha_identity(guide = "none") +
    scale_color_continuous_diverging("Blue-Red 3", breaks = c(-7, 0, 7)) +
    theme(
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      legend.position = c(1, 0.75),
      legend.key.width = unit(1 / 16, "cm"),
      legend.key.height = unit(1 / 6, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
    ) +
    labs(
      x = x_label,
      y = y_label,
      color = color_label,
      fill = color_label
    ) +
    coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1))

  return(p_icc_scatter)
}



label_regions <- function(labels, width = 18) {
  labels <- gsub("17Networks_", "", labels) %>% gsub("_", " ", .)
  labels <- stringr::str_wrap(labels, width = width)
  return(labels)
}
