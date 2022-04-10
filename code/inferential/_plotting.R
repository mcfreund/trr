library(here)
library(readr)
library(tidyr)
library(dplyr)
library(mfutils)
library(ggplot2)
library(cowplot)
library(viridis)
library(patchwork)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)

# Constants
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"

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

# Atlas
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
  rois <- get(atlas_nm)$key[[roi_col]]
  atlas <- schaefer17_400
  atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
  stop("not configured for atlas")
}

# Plotting function
brain_plot <- function(df, eff_term = NULL, eff = NULL, stat_term = "tstat",
                    lim = NULL, direct = 1, fig_title = "Example figure", savename = NULL) {
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
  if (is.null(savename)) {
    message()
    message("Note: you need to use 'print(brain_plot(...))'",
      " instead of 'brain_plot(...)' to display the result.")
    message()
  } else {
    ggsave(savename, plot = fig)
  }
  return(fig)
}


# Main function when running this script directly
if (sys.nframe() == 0) {
  fname1 <- "multivariate_linear_model.csv"
  b1 <- read_csv(here("out", "spatial", fname1)) %>%
    mutate(term = ifelse(term == "hilo_alllo", "hilo_allhi", term)) %>%
    mutate(b = ifelse(term == "hilo_allhi", -b, b),
      tstat = ifelse(term == "hilo_allhi", -tstat, tstat))
  f1 <- brain_plot(b1, eff_term = "term", eff = "hilo_allhi", lim = c(-6, 12),
    fig_title = "multivariate")

  fname2 <- "univariate_linear_model.csv"
  b2 <- read_csv(here("out", "spatial", fname2))
  f2 <- brain_plot(b2, eff_term = "term", eff = "hilo_allhi", lim = c(-6, 12),
    fig_title = "univariate")

  f3 <- brain_plot(b1, eff_term = "term", eff = "wavewave2", lim = c(-2, 2),
    fig_title = "multivariate")
  f4 <- brain_plot(b2, eff_term = "term", eff = "wavewave2", lim = c(-2, 2),
    fig_title = "univariate")
  f3 + f4 + plot_annotation(title = "t-statistics for wave effect in Stroop baseline")
  ggsave(here("out", "spatial", "wave.png"))
}