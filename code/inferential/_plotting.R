# Plotting function & Generate figures
#
# Author: Ruiqi Chen
#
# 07/29/2022 update: remove all unused codes
#
# When being sourced, this script provides a function `brain_plot()` that can plot
# some statistics for each parcel over the brain. Please refer to the comments above
# the definition of `brain_plot()` for its usage.
#
# When executing directly, the script converts an `.rmd` file to an `.md` report under
# the `reports/` directory. Note that it will change the working directory to `/reports`
# during knitting.

library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(mfutils)

# ROIs
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"

# Atlas
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
}

# Theme
theme_set(theme_bw(base_size = 12))
theme_surface <- list(
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_text(size = 7),
    legend.background = element_blank(), legend.text = element_text(size = 7), legend.direction = "horizontal",
    legend.key.height = unit(1 / 4, "cm"), legend.key.width = unit(1 / 3, "cm")
  )
)

# **brain_plot(): Plot a column in a tibble onto the brain**
#
# Inputs:
# - `df`: a tibble, the names of the parcels should be in `df$region`
# - `stat_term`: a string, name of the column in `df` to plot
# - `eff_term` and `eff`: `eff_term` is a string, `eff` can be a string or a sequence of string,
#     if specified, `df` will be filtered by `df[[eff_term]] %in% eff`
# - `lim`: a sequence of length 2, the limit of colormap to use (by default automatically determined)
# - `direct`: can set as -1 to reverse the colormap
# - `fig_title`: a string, title of the plot
# - `savename`: NULL (don't save) or a string, full path to save the figure
#
# Output:
# - `fig`: A `ggplot2()` figure plotted with `geom_brain()` from package `ggseg`. You sometimes
#     need to `print()` it to display the figure.
#
brain_plot <- function(df, stat_term = "Estimate", eff_term = NULL, eff = NULL,
                    lim = NULL, direct = 1, fig_title = "", savename = NULL) {
  if (!is.null(eff_term)) {
    df <- df %>% filter(.data[[eff_term]] %in% .env$eff)
    # df <- df %>% group_by(.data[[eff_term]])  # Will fail due to error in brain_join()
  }
  fig <- df %>%
    ggplot() +
    geom_brain(aes(fill = .data[[stat_term]]),
      atlas = atlas, position = position_brain(side ~ hemi)) +
    scale_fill_viridis_c(
      limits = lim,
      direction = direct,
      option = "magma", na.value = "grey",
      breaks = scales::extended_breaks(4)
      ) +
    theme_surface +
    labs(title = fig_title, fill = NULL)

  # Saving or displaying
  if (!is.null(savename)) {
    ggsave(savename, plot = fig)
  }
  return(fig)
}

# **get_diff_data() - calculate the difference between terms**
#
# Inputs:
# - x: a tibble
# - name_term: "model" or "response"
# - base_level: base level for comparison between `name_term`s, e.g. "uv". By default
#   it's the first value in `name_term` if `name_term` is not factor, or it's the first
#   level that appears in `name_term` if `name_term` is a factor.
# - val_term: value term, by default "Estimate"
# - id_term: terms that uniquely define an observation, by default all other terms
# - pivoting: "long" or "wide", by default "long"
#
# Ouput:
# - res: a tibble with `name_term`s replaced by the difference between each level
#   (expcept base level) and base level.
get_diff_dat <- function(x, name_term, base_level = NA, val_term = "Estimate",
  id_term = NULL, pivoting = "long") {

  if (!is.null(id_term)) {
    res <- select(x, .env$name_term, .env$val_term, .env$id_term)
  } else res <- x
  if (is.factor(x[[name_term]])) {  # Sort the column according to levels
    res <- res %>% arrange(.env$name_term)
  }
  term_levels <- unique(res[[name_term]])
  if (is.na(base_level)) base_level <- term_levels[[1]]

  res <- pivot_wider(res, names_from = .env$name_term, values_from = .env$val_term)
  for (term_level in term_levels) {
    if (term_level != base_level) res[[paste(term_level, "-", base_level)]] <- (
      res[[term_level]] - res[[base_level]])
  }
  res <- select(res, !.env$term_levels)

  if (pivoting == "long") {
    new_names <- paste(term_levels[term_levels != base_level], "-", base_level)
    res <- pivot_longer(res, all_of(new_names),
      names_to = name_term, values_to = val_term)
    if (is.factor(x[[name_term]])) res[[name_term]] <- factor(res[[name_term]], new_names)
  }

  res

}

# Main function - knitting an rmarkdown report
#
# Note: It will change the working directory to the directory of the rmarkdown file.
# Make sure the rmarkdown file works in this way before knitting it!
if (sys.nframe() == 0) {
  f <- "model_comparison.rmd"
  cwd <- getwd()
  setwd(here::here("reports"))
  knitr::knit(f)
  setwd(cwd)
}