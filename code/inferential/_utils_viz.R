library(ggplot2)
library(dplyr)
library(rlang)

arsinh <- scales::trans_new("arsinh", transform = function(x) asinh(x), inverse = function(x) sinh(x))


highlight_overlay <- function(
  x, upper,
  lower = min(x),
  upper_alpha = 1,
  lower_alpha = 0.1,
  scaling_fun = function(x) x^2
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


create_colorgrid <- function(
  colorscale, limits_color_stat, limits_alpha_stat,
  n = 50,
  ...
) {
  color_stat <- seq(limits_color_stat[1], limits_color_stat[2], length.out = n)
  alpha_stat <- seq(limits_alpha_stat[1], limits_alpha_stat[2], length.out = n)
  clrs <- colorscale$map(color_stat)
  colorgrid <- expand.grid(color = clrs, alpha_stat = alpha_stat)
  color_stat <- setNames(color_stat, clrs)
  colorgrid$color_stat <- color_stat[colorgrid$color]
  colorgrid$alpha <- highlight_overlay(colorgrid$alpha_stat, ...)
  as_tibble(colorgrid)
}



label_regions <- function(labels, width = 18) {
  labels <- gsub("17Networks_", "", labels) %>% gsub("_", " ", .)
  labels <- stringr::str_wrap(labels, width = width)
  return(labels)
}
label_regions_eg <- function(x) paste0("q", 1:4, "\n", label_regions(eg))



pivot_and_contrast <- function(
  data, contrast_expr, id_cols, names_from, values_from,
  new_col = "difference",
  ...
) {
  ## compute contrast:

  contrast <- rlang::parse_expr(contrast_expr)
  
  data <- data %>%
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = all_of(names_from),
      values_from = all_of(values_from),
      ...
    ) %>%
    mutate(!!new_col := !!contrast)
  
  data
}


plot_surface <- function(
  data,
  statistic_col,
  border_color_col = "is_roi",
  border_size_col = "is_roi",
  border_color_values = color_line_roi,
  border_size_values = size_line_roi,
  alpha_col = "alpha_level",
  position = ggseg::position_brain(. ~ hemi + side),
  underlay = c(color = "grey", fill = "grey"),
  atlas_data = atlas,
  scale_fill = function(...) scale_fill_viridis_c(option = "magma", na.value = "white", ...),
  fill_label = "statistic",
  theme_ = theme_surface(strip.text = element_text(size = rel(1))),
  guides_ = guides(
    fill = guide_colorbar(title.position = "left", title.vjust = 0.9),
    alpha = "none", color = "none", size = "none"
  ),
  facet_factor = NULL,
  facet_factor_order,
  ...
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  ## group_by if facetting:
  if (!is_null(facet_factor)) {
    facet_factor_sym <- rlang::parse_expr(facet_factor)
    data <- data %>%
      group_by(.data[[facet_factor]]) %>%
      mutate({{facet_factor_sym}} := factor({{facet_factor_sym}}, levels = names(facet_factor_order)))
  }
  
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
    p <- p + 
      ggseg::geom_brain(aes(alpha = .data[[alpha_col]]), atlas = atlas_data, position = position) +
      scale_alpha_identity()
  } else {
    p <- p + ggseg::geom_brain(atlas = atlas_data, position = position)
  }

  ## facet
  
  if (!is_null(facet_factor)) {
    facet_factor_sym <- rlang::parse_expr(facet_factor)
    p <- p + 
      facet_grid(
        vars({{facet_factor_sym}}),
        labeller = labeller(
          {{facet_factor_sym}} := setNames(facet_factor_order, names(facet_factor_order))
        ),
        switch = "y")
  }

  ## scales, guides, themes
  p <- p + scale_fill(...) + labs(fill = fill_label) + guides_ + theme_

  p

}



plot_hist_means <- function(
  data, value_col, fill_col, text_label, colors,
  alpha_level = 0.5,
  n_bins = 10,
  text_y = c(140, 75),
  text_x = c(-1/2, -1/2),
  text_size = 4,
  x_lab = "statistic",
  y_lab = "Number of parcels"
) {

  ## create object to hold colored text labels
  text_data <- data.frame(label = text_label, y = text_y)
  text_data[[value_col]] <- text_x
  text_data[[fill_col]] <- names(text_label)

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
  data,
  value_col = "difference",
  fill_col = "is_roi",
  text_label = c("'ROIs'" = TRUE, "all" = FALSE),
  colors = colors_roi,
  text_x = c(-0.9, -0.9),
  text_y = c(50, 70),
  text_size = 4,
  n_bins = 10,
  x_breaks = c(-1, 0, 1),
  x_lab = "statistic",
  y_lab = "Number of parcels"
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  ## create object to hold colored text labels
  text_data <- data.frame(label = names(text_label), y = text_y)
  text_data$difference <- text_x
  text_data[[fill_col]]  <- text_label

  data %>%
    ggplot(aes(.data[[value_col]], fill = .data[[fill_col]], color = .data[[fill_col]])) +
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
  data, x_col, y_col, color_col, axis_labs,
  linetype = "dotted",
  alpha_level = 0.75,
  point_size = 1,
  breaks_x = c(-1, 0, 1),
  breaks_y = c(-1, 0, 1),
  limits_x = c(-1, 1),
  limits_y = c(-1, 1),
  legend_position = c(0.5, -0.5),
  legend_direction = "horizontal",
  color_lab = univ_stat_lab,
  scale_color = colorspace::scale_color_continuous_diverging("Blue-Red 3", breaks = c(-6, 0, 6)),
  add_zero_lines = TRUE
) {

  source(here::here("code", "inferential", "_parameters_viz.R"))
  color_sym <- rlang::parse_expr(color_col)

  p <- data %>%
    ## sort by color col to control plotting order of points
    as_tibble %>%
    ungroup %>%
    arrange({{ color_sym }}) %>%
    ggplot(aes(.data[[x_col]], .data[[y_col]])) +
    geom_abline(linetype = linetype, color = "grey50")
  
  if (add_zero_lines) {
    p <- p +
      geom_vline(xintercept = 0, linetype = linetype, color = "grey50") +
      geom_hline(yintercept = 0, linetype = linetype, color = "grey50")
  }
  
  p <- p +
    geom_point(aes(color = .data[[color_col]]), alpha = alpha_level, stroke = 0, size = point_size) +
    scale_color +
    theme(
      legend.position = legend_position,
      legend.direction = legend_direction,
      legend.key.height = unit(1 / 8, "cm"),
      legend.key.width = unit(1 / 4, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
    ) +
    labs(
      x = axis_labs[names(axis_labs) == x_col],
      y = axis_labs[names(axis_labs) == y_col],
      color = color_lab
    ) +
    scale_x_continuous(breaks = breaks_x, limits = limits_x) +
    scale_y_continuous(breaks = breaks_y, limits = limits_y)

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
  width = 183,
  height = 110,
  dev = cairo_pdf,
  plot_annotations = NULL,
  units = "mm",
  ...
) {

  p <- (patchwork::free(brains) / lower_panels) + patchwork::plot_layout(design = plot_layout)

  if (!is.null(plot_annotations)) {
    p <- p + plot_annotations
  }

  if (!is.null(filename) && !is.null(path)) {
    ggsave(
      file.path(path, paste0(filename, ".pdf")),
      p, dev = cairo_pdf, height = height, width = width, units = units, ...)
  }

  p

}


create_args <- function(plot_func, ...) {
  # Get the default arguments of the plotting function
  default_args <- formals(plot_func)
  #browser()
  
  # Create an environment from the default arguments
  args_env <- list2env(as.list(default_args))
  
  # Update the environment with any additional arguments provided
  list2env(list(...), envir = args_env)
  
  # Return the updated arguments as a list
  as.list(args_env)
}



arrange_figure_comparison <- function(
  data, contrast_expr, comparison_factor, comparison_factor_order,
  title, path, filename,
  params_plot_surface = NULL,
  params_plot_hist_means = NULL,
  params_plot_hist_diff = NULL,
  params_plot_scatter = NULL,
  id_cols = c("region", "tplus", "q05", "is_roi"),
  value_col = "value"
) {
  source(here::here("code", "inferential", "_parameters_viz.R"))
  
  ## order factor levels (defines order of surface rows):
  comparison_factor_sym <- rlang::parse_expr(comparison_factor)
  data <- data %>% 
    mutate({{comparison_factor_sym}} := factor({{comparison_factor_sym}}, levels = names(comparison_factor_order)))

  ## create data for hist diff and scatter: pivot to wide and contrast
  data_wide <- pivot_and_contrast(
    data,
    contrast_expr = contrast_expr,  ## evaluate in context of data
    id_cols = id_cols,
    names_from = comparison_factor,
    values_from = value_col
  )

  ## create argument lists based on defaults, common args, and unique args (params_*)
  ## common:
  new_surface <- list(plot_func = plot_surface,
    data = data, statistic_col = value_col, facet_factor = comparison_factor,
    facet_factor_order = comparison_factor_order
  )
  new_hist_means <- list(plot_func = plot_hist_means,
    data = data, value_col = value_col, fill_col = comparison_factor
  )
  new_hist_diff <- list(plot_func = plot_hist_diff, data = data_wide)
  new_scatter <- list(
    plot_func = plot_scatter, data = data_wide,
    x_col = names(comparison_factor_order)[1],
    y_col = names(comparison_factor_order)[2]
  )
  
  ## add unique and defaults:
  args_plot_surface <- do.call(create_args, c(new_surface, params_plot_surface))
  args_plot_hist_means <- do.call(create_args, c(new_hist_means, params_plot_hist_means))
  args_plot_hist_diff <- do.call(create_args, c(new_hist_diff, params_plot_hist_diff))
  args_plot_scatter <- do.call(create_args, c(new_scatter, params_plot_scatter))
    
  # Generate plots:
  p_brains <- do.call(plot_surface, args_plot_surface)
  p_means <- do.call(plot_hist_means, args_plot_hist_means)
  p_diff <- do.call(plot_hist_diff, args_plot_hist_diff)
  p_scatter <- do.call(plot_scatter, args_plot_scatter)
  
  # Arrange and save the plots
  p <- arrange_plots(
    p_brains,
    p_means + p_diff + p_scatter,
    plot_annotations = patchwork::plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5))),
    filename = filename,
    path = path
  )

  list(
    brains = p_brains, means = p_means, diff = p_diff, scatter = p_scatter,
    arranged = p
  )

}


