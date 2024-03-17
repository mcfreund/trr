library(ggplot2)
library(dplyr)
library(rlang)

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


label_regions <- function(labels, width = 18) {
  labels <- gsub("17Networks_", "", labels) %>% gsub("_", " ", .)
  labels <- stringr::str_wrap(labels, width = width)
  return(labels)
}
label_regions_eg <- \(x) paste0("q", 1:4, "\n", label_regions(eg))



plot_surface <- function(
  data,
  statistic_col,
  border_color_col = "is_roi",
  border_size_col = border_color_col,
  border_color_values = color_line_roi,
  border_size_values = size_line_roi,
  alpha_col = "alpha_level",
  position = ggseg::position_brain(. ~ hemi + side),
  underlay = c(color = "grey", fill = "grey"),
  atlas_data = atlas,
  scale_fill = \(...) scale_fill_viridis_c(option = "magma", na.value = "white", ...),
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
  p <- p + scale_fill(...) + guides_ + theme_

  p

}



plot_hist_means <- function(
  data, value_col, fill_col, text_label, colors,
  alpha_level = 0.5,
  n_bins = 10,
  text_y = c(125, 75),
  text_x = c(-2/3, -2/3),
  text_size = 4,
  x_lab = "statistic",
  y_lab = "Number of parcels"
) {
  
  ## create object to hold colored text labels
  text_data <- data.frame(label = names(text_label), y = text_y)
  text_data[[value_col]] <- text_x
  text_data[[fill_col]]  <- text_label

  data %>%
    ggplot(aes(.data[[value_col]], fill = .data[[fill_col]])) +
    geom_histogram(position = "identity", alpha = alpha_level, color = "black", bins = n_bins) +
    geom_text(data = text_data, aes(y = y, label = label, color = .data[[fill_col]]), size = text_size) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(x = x_lab, y = y_lab) +
    theme(legend.position = "none")

}



plot_hist_diff <- function(
  data, id_cols, names_from_col, values_from_col, contrast,
  fill_col = "is_roi",
  text_label = c("'ROIs'" = TRUE, "all parcels" = FALSE),
  colors = colors_roi,
  text_x = c(-1, -1),
  text_y = c(50, 70),
  text_size = 4,
  n_bins = 10,
  x_breaks = c(-1, 0, 1),
  x_lab = "statistic",
  y_lab = "Number of parcels"
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  ## compute contrast:
  data <- data %>%
    as_tibble() %>%
    pivot_wider(id_cols = id_cols, names_from = names_from_col, values_from = values_from_col) %>%
    mutate(difference := {{ contrast }})

  ## create object to hold colored text labels
  text_data <- data.frame(label = names(text_label), y = text_y)
  text_data$difference <- text_x
  text_data[[fill_col]]  <- text_label

  data %>%
    ggplot(aes(difference, fill = .data[[fill_col]], color = .data[[fill_col]])) +
    geom_histogram(position = "identity", fill = colors[["FALSE"]], color = "black", bins = n_bins) +
    geom_histogram(
      data = . %>% filter(.data[[fill_col]]),
      fill = colors[["TRUE"]],
      color = "black",
      position = "identity",
      bins = n_bins
    ) +
    geom_text(data = text_data, aes(y = y, label = label), size = text_size) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = x_breaks) +
    labs(x = x_lab, y = y_lab) +
    theme(legend.position = "none")

}




plot_scatter <- function(
  data, x_col, y_col,
  id_cols, names_from, values_from,
  color_col = "q05",
  linetype = "dashed",
  alpha_level = 0.75,
  point_size = 1,
  breaks = c(-1, 0, 1),
  x_lab = "statistic",
  y_lab = "Number of parcels",
  color_lab = "d'",
  scale_color = colorspace::scale_color_continuous_diverging("Blue-Red 3", breaks = c(-7, 0, 7)),
  pivot = TRUE
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  if (pivot) {
    data <- pivot_wider(data, id_cols = all_of(id_cols), names_from = all_of(names_from), values_from = all_of(values_from))
  }

  data %>%
    ## sort by color col to control plotting order of points
    arrange({{ color_col }}) %>%
    ggplot(aes(.data[[x_col]], .data[[y_col]])) +
    geom_abline(linetype = linetype) +
    geom_vline(xintercept = 0, linetype = linetype) +
    geom_hline(yintercept = 0, linetype = linetype) +
    geom_point(aes(color = .data[[color_col]]), alpha = alpha_level, stroke = 0, size = point_size) +
    scale_color +
    theme(
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      legend.position = c(1, 0.75),
      legend.key.width = unit(1 / 16, "cm"),
      legend.key.height = unit(1 / 6, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
    ) +
    labs(x = x_lab, y = y_lab, color = color_lab) +
    scale_x_continuous(breaks = breaks) +
    scale_y_continuous(breaks = breaks)

}


arrange_plots <- function(
  brains,
  lower_panels,
  plot_layout = c(
    patchwork::area(t = 0, b = 65, l = 0, r = 220),
    patchwork::area(t = 66, b = 100, l = 20, r = 200)
  ),
  filename = NULL,
  path = NULL,
  width = 8,
  height = 4.5,
  dev = cairo_pdf,
  plot_annotations = NULL,
  ...
) {

  p <- (patchwork::free(brains) / lower_panels) + patchwork::plot_layout(design = plot_layout)

  if (!is.null(plot_annotations)) {
    p <- p + plot_annotations
  }

  if (!is.null(filename) && !is.null(path)) {
    ggsave(file.path(path, paste0(filename, ".pdf")), p, dev = cairo_pdf, height = height, width = width, ...)
  }

  p

}


arrange_figure_comparison <- function(
  data, comparison_factor, comparison_factor_order, colors_comparison,
  margins, title, path_figs_results, id_cols, color_col, color_lab, contrast_expr, filename,
  alpha_level = 0.5,
  value_col = "value",
  breaks = c(-1, 0, 1),
  limits = c(-1, 1),
  fill_labs = "TRR (r)"
) {
  

  # Set fill_col and names_from to be the character version of comparison_factor
  fill_col <- as_string(ensym(comparison_factor))
  names_from <- as_string(ensym(comparison_factor))
  
  # Set x_col and y_col based on comparison_factor_order
  x_col <- comparison_factor_order[1]
  y_col <- comparison_factor_order[2]

  # Generate the surface plot
  p_brains <- data %>%
    mutate({{comparison_factor}} := factor({{comparison_factor}}, levels = comparison_factor_order)) %>%
    group_by({{comparison_factor}}) %>%
    plot_surface(
      statistic = value_col,
      limits = limits,
      breaks = breaks,
      alpha_col = NULL,
      underlay = NULL
    ) +
    labs(fill = fill_labs) +
    facet_grid(
      vars({{comparison_factor}}),
      labeller = labeller(
        {{comparison_factor}} := setNames(names(comparison_factor_order), comparison_factor_order)
        ),
      switch = "y"
    )
  
  # Generate the histogram of means
  p_means <- plot_hist_means(
    data,
    value_col = value_col,
    fill_col = fill_col,
    alpha_level = alpha_level,
    colors = colors_comparison,
    text_label = comparison_factor_order
  ) +
  theme(plot.margin = unit(margins, "points"))
  
  # Generate the histogram of differences
  p_diff <- plot_hist_diff(
    data,
    id_cols = id_cols,
    names_from = names_from,
    values_from = value_col,
    contrast = eval(parse(text = contrast_expr))
  ) +
  theme(plot.margin = unit(margins, "points"))
  
  # Generate the comparison scatterplot
  p_scatter <- plot_scatter(
    data,
    x_col = x_col,
    y_col = y_col,
    id_cols = id_cols,
    names_from = names_from,
    values_from = value_col,
    color_col = color_col,
    color_lab = color_lab
  ) +
  theme(plot.margin = unit(margins, "points"))
  
  # Arrange and save the plots
  p <- arrange_plots(
    p_brains,
    p_means + p_diff + p_scatter,
    plot_annotations = patchwork::plot_annotation(title = title),
    filename = filename,
    path = path_figs_results
  )

  list(
    brains = p_brains, means = p_means, diff = p_diff, scatter = p_scatter,
    arranged = p
    )

}



pivot_and_contrast <- function(
  data, contrast, id_cols, names_from, values_from,
  new_col = difference,
  ...
  ) {

  ## compute contrast:
  data <- data %>%
    pivot_wider(
      id_cols = {{ id_cols }},
      names_from = {{ names_from }},
      values_from = {{ values_from }},
      ...
      ) %>%
    mutate({{ new_col }} := {{ contrast }})
  
  data

}
