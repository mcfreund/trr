library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(data.table)
library(brms)
library(posterior)
library(colorspace)
library(mfutils)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(patchwork)

source(here("code", "inferential", "_parameters_viz.R"))
source(here("code", "inferential", "_parameters_reliability.R"))
source(here("code", "inferential", "_utils_viz.R"))
source(here("code", "inferential", "_utils_reliability.R"))
source(here("code", "_atlas.R"))
source(here("code", "_constants.R"))
source(here("code", "_paths.R"))

## other parameters and constants

selected_model <- "no_lscov_symm"
selected_session <- "baseline"
stopifnot(length(selected_model) == 1 && length(selected_session) == 1)

dir.create(path_figs_results_supp, recursive = TRUE, showWarnings = FALSE)


## read data ----

posterior_samples <- load_summaries(
  data_type = "posterior_samples",
  prefixes = c("popef", "ratio", "trr"),
  session = selected_session,
  models_order = selected_model,
  sum_fun = subset_and_order
)

posterior_summaries <- load_summaries(
  data_type = "posterior_summaries",
  prefixes = c("popef", "ratio", "trr"),
  sess = selected_session,
  models_order = selected_model,
  sum_fun = subset_and_order
)

summarystats <- load_summaries(
  data_type = "summarystat",
  prefixes = c("popef", "trr"),
  session = selected_session,
  sum_fun = subset_and_order,
  models_order = NULL
)


## minor computations ----

## get dprime and ROIs column
popef <-
  posterior_summaries$popef %>%
  filter(response == "uv", sum_fun %in% c("mean", "sd")) %>%
  dcast(... ~ sum_fun, value.var = "value") %>%
  mutate(
    dprime = mean / sd,
    is_roi = quantile(dprime, 0.9) < dprime,
    alpha_roi = highlight_overlay(dprime, upper = quantile(dprime, 0.9), lower = min(dprime), lower_alpha = 0),
    network = get_network(region)
  )

popef_network <- popef %>%
  group_by(network) %>%
  summarize(
    dprime = mean(dprime),
    n_parcel = n(),
    n_roi = sum(is_roi),
    prop_roi = mean(is_roi)
  )


## format for plotting ----

## reorder popef:
popef <- popef %>% mutate(network = factor(network, levels = popef_network$network[order(popef_network$dprime)]))
popef_network <- popef_network %>%
  mutate(network = factor(network, levels = popef_network$network[order(popef_network$dprime)]))

## create dataframes for each plot:
uv_mv_icc_hbm <-
  bind_rows(
    summarystats$trr %>% rename(value = r) %>% mutate(sum_fun = "pointest", statistic = "trr"),
    posterior_summaries$trr
  ) %>%
  full_join(popef %>% select(region, dprime, is_roi, alpha_roi), by = "region")


## population-level contrast on univariate response ----

## plot surfaces

p_pop_brain_thresh <-
  popef %>%
  plot_surface(statistic = "dprime") +
  labs(fill = "d'") +
  theme(legend.position = c(0.5, 0))
p_pop_brain_unthresh <-
  popef %>%
  plot_surface(statistic = "dprime", underlay = NULL, alpha_col = NULL) +
  labs(fill = "d'") +
  theme(legend.position = c(0.5, 0))

ggsave(file.path(path_figs_results_supp, "pop_dprime_thresh.pdf"), p_pop_brain_thresh, height = 2, width = 4.5)
ggsave(file.path(path_figs_results_supp, "pop_dprime_unthresh.pdf"), p_pop_brain_unthresh, height = 2, width = 4.5)

## check stats within network assignments

p_parcel_dprime <-
  popef %>%
  ggplot(aes(dprime)) +
  geom_histogram(bins = n_bins_rois) +
  geom_vline(xintercept = quantile(popef$dprime, 0.9)) +
  labs(x = "Population-Level\nStroop Effect (d')", y = "Number of parcels")

p_network_dprime <-
  popef %>%
  ggplot(aes(dprime, network)) +
  geom_vline(xintercept = quantile(popef$dprime, 0.9)) +
  geom_boxplot(fill = "grey", color = "black", width = 0.2) +
  theme(
    axis.line.y.left = element_blank(),
    panel.grid.major.y = element_line(linewidth = rel(0.25), linetype = "dashed")
  ) +
  labs(y = "Schaefer 400-17 Network", x = "Population-Level\nStroop Effect (d')")

p_network_props <-
  popef_network %>%
  ggplot(aes(prop_roi, network)) +
  geom_bar(fill = "grey", color = "black", width = 0.5, stat = "identity") +
  geom_text(aes(label = paste0(n_roi, "/", n_parcel), x = 0, y = network), hjust = 0, size = 2) +
  theme(axis.line.y.left = element_blank()) +
  labs(y = "Schaefer 400-17 Networks\nwith 'ROI' parcels", x = "Proportion of\n'ROI' parcels")

ggsave(
  file.path(path_figs_results_supp, "pop_network_dprime.pdf"),
  p_parcel_dprime + p_network_dprime + p_network_props,
  height = 3, width = 9
)


## Impact of HBM on TRR | Univariate, Multivariate ----

## Univariate: ICC versus HBM

## surface plot
p_trr_uv_brains <-
  uv_mv_icc_hbm[response == "uv" & sum_fun %in% c("pointest", "map")] %>%
  mutate(model_nm = factor(model_nm, levels = c("summarystat", selected_model))) %>%
  group_by(model_nm) %>%
  plot_surface(
    statistic = "value",
    limits = c(-0.4, 0.9),
    breaks = c(0, 0.9),
    alpha_col = NULL
  ) +
  labs(fill = "TRR (r)") +
  facet_grid(
    vars(model_nm),
    labeller = labeller(model_nm = setNames(c("ICC", "HBM\n(MAP)"), c("summarystat", selected_model))), switch = "y"
  )

## histogram of means
p_trr_uv_means <- plot_hist_means(
  uv_mv_icc_hbm[response == "uv" & sum_fun %in% c("pointest", "map")],
  value_col = "value",
  fill_col = "model_nm",
  alpha_level = 0.5,
  colors = colors_models_comparison,
  text_label = models_comparison,
)

## histogram of differences
p_trr_uv_diff <- plot_hist_diff(
  uv_mv_icc_hbm[response == "uv" & sum_fun %in% c("pointest", "map")],
  id_cols = c("region", "is_roi", "alpha_roi", "dprime"),
  names_from = "model_nm",
  values_from = "value",
  contrast = atanh(no_lscov_symm) - atanh(summarystat),
  x_lab = "\u0394ICC(Stroop):\nmultiv. \u2212 univar. (z)"
)

## comparison scatterplot
p_trr_uv_scatter <-
  uv_mv_icc_hbm[response == "uv" & sum_fun %in% c("pointest", "map")] %>%
  pivot_wider(id_cols = c(region, is_roi, dprime, alpha_roi), names_from = model_nm, values_from = value) %>%
  plot_scatter(
    x_col = "no_lscov_symm",
    y_col = "summarystat",
    color_col = "dprime",
    x = "HBM MAP(TRR): (r)",
    y = "ICC TRR: Univariate (r)",
    color_lab = "d'"
  )

## arrange and save
p_trr_uv <- arrange_plots(
  p_trr_uv_brains,
  p_trr_uv_means,
  p_trr_uv_diff,
  p_trr_uv_scatter,
  filename = "trr_uv",
  path = path_figs_results
)

p_trr_uv

