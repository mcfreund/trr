library(here)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)

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
atlas_nm <- "schaefer2018_17_400_fsaverage5"
if (atlas_nm == "schaefer2018_17_400_fsaverage5") {
    atlas <- schaefer17_400
    atlas$data$region <- gsub("^lh_|^rh_", "", atlas$data$label)
} else {
    stop("not configured for atlas")
}

# Plotting function
brain_plot <- function(df, eff = "hilo_alllo", lim = c(-12, 3), direct = -1,
                        fig_title = "t-statistics for hi-lo, multivariate method, stroop, baseline",
                        savename = NULL, eff_term = "term", stat_term = "tstat") {

    fig <- df %>%
        filter(.data[[eff_term]] %in% .env$eff) %>%
        # group_by(.data[[eff_term]]) %>%  # Doesn't work due to error in brain_join()
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
        # facet_wrap(as.formula(paste0("~", eff_term))) +  # Will generate NA for non-regions in the atlas
        labs(title = fig_title, fill = NULL)

    # Saving or displaying
    if (is.null(savename)) {
        print(fig)
        readline("The ggplot will disappear after function call. Press enter to continue.")
    } else {
        ggsave(savename, plot = fig, dpi = "screen")
    }
}


# Main function when running this script directly
if (sys.nframe() == 0) {
    fname <- "multivariate_linear_model.csv"
    b <- read_csv(here("out", "spatial", fname))
    brain_plot(b)
}