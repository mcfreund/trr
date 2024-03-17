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

## get population stats and ROIs columns

popef <-
  posterior_summaries$popef %>%
  filter(response == "uv", sum_fun %in% c("mean", "sd", "map", "q05")) %>%
  dcast(... ~ sum_fun, value.var = "value") %>%
  mutate(
    is_roi = quantile(q05, 0.9) < q05,
    alpha_level = highlight_overlay(mean, upper = quantile(mean, 0.9), lower = min(mean), lower_alpha = 0),
    tplus = mean / sd,
    is_roi_mean = quantile(mean, 0.9) < mean,
    is_roi_map = quantile(map, 0.9) < map,
    is_roi_tplus = quantile(tplus, 0.9) < tplus,
    network = get_network(region)
  )

## examine similarities:
if (FALSE) {
  sum(popef$is_roi & popef$is_roi_mean)  ## agreement (out of 40)
  sum(popef$is_roi & popef$is_roi_tplus)
  pairs(popef[, c("map", "mean", "q05", 'sd', "tplus")])
}

popef_network <- popef %>%
  group_by(network) %>%
  summarize(
    q05 = mean(q05),
    n_parcel = n(),
    n_roi = sum(is_roi),
    prop_roi = mean(is_roi)
  )



## format for plotting ----

## reorder popef:
popef <- popef %>% mutate(network = factor(network, levels = popef_network$network[order(popef_network$q05)]))
popef_network <- popef_network %>%
  mutate(network = factor(network, levels = popef_network$network[order(popef_network$q05)]))

## create dataframes for each plot:
uv_mv_icc_hbm <-
  bind_rows(
    summarystats$trr %>% rename(value = r) %>% mutate(sum_fun = "pointest", statistic = "trr"),
    posterior_summaries$trr
  ) %>%
  full_join(popef %>% select(region, q05, is_roi, alpha_level), by = "region")


## population-level contrast on univariate response ----

## plot surfaces

p_pop_brain_thresh <-
  popef %>%
  plot_surface(
    statistic = "mean",
    #limits = c(-0.1, 0.13),
    breaks = c(-0.1, 0, 0.1)
  ) +
  labs(fill = "d") +
  theme(legend.position = c(0.5, 0))
p_pop_brain_unthresh <-
  popef %>%
  plot_surface(
    statistic = "mean",
    #limits = c(. %>% min(mean), . %>% max(mean)),
    breaks = c(-0.1, 0, 0.1),
    underlay = NULL,
    alpha_col = NULL
  ) +
  labs(fill = "d") +
  theme(legend.position = c(0.5, 0))

ggsave(file.path(path_figs_results_supp, "pop_d_thresh.pdf"), p_pop_brain_thresh, height = 2, width = 4.5)
ggsave(file.path(path_figs_results_supp, "pop_d_unthresh.pdf"), p_pop_brain_unthresh, height = 2, width = 4.5)

## check stats within network assignments

thresh_roi <- unname(quantile(popef$q05, 0.9))

p_parcel_q05 <-
  popef %>%
  ggplot(aes(q05, fill = is_roi)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = quantile(popef$q05, 0.9)) +
  labs(x = "Population-Level\nUnivar. Stroop Effect (d)", y = "Number of parcels") +
  scale_fill_manual(values = colors_roi) +
  scale_x_continuous(breaks = c(-0.1, 0.1, round(thresh_roi, 2))) +
  theme(legend.position = "none")

p_network_q05 <-
  popef %>%
  ggplot(aes(q05, network)) +
  geom_vline(xintercept = thresh_roi) +
  geom_boxplot(fill = "grey", color = "black", width = 0.4) +
  theme(
    axis.line.y.left = element_blank(),
    panel.grid.major.y = element_line(linewidth = rel(0.25), linetype = "dashed"),
    axis.text.y = element_text(size = rel(0.75))
  ) +
  scale_x_continuous(breaks = c(-0.1, 0.1, round(thresh_roi, 2))) +
  labs(y = "Cortical Network\n(Schaefer 400-17)", x = "Population-Level\nStroop Effect (d')")

p_network_props <-
  popef_network %>%
  mutate(is_roi = TRUE) %>%
  ggplot(aes(prop_roi, network, fill = is_roi)) +
  geom_bar(color = "black", width = 0.6, stat = "identity") +
  geom_text(
    aes(label = paste0(n_roi, "/", n_parcel), y = network, x = prop_roi + 0.002),
    hjust = 0,
    color = "black",
    size = 1.75
  ) +
  scale_fill_manual(values = colors_roi) +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)) +
  #scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  theme(
    axis.line.y.left = element_blank(),
    axis.text.y = element_text(size = rel(0.75)),
    legend.position = "none"
  ) +
  labs(y = "Cortical Network\n(Schaefer 400-17)", x = "Proportion of\n'ROI' parcels")

ggsave(
  file.path(path_figs_results_supp, "pop_network_q05.pdf"),
  p_parcel_q05 + p_network_q05 + p_network_props,
  height = 3, width = 9
)

arrange_plots(
  p_pop_brain_unthresh,
  p_parcel_q05 + p_network_q05 + p_network_props,
  filename = "pop_uv",
  path = path_figs_results,
  plot_layout = c(
    area(t = 0, b = 62, l = 0, r = 220),
    area(t = 55, b = 100, l = 20, r = 200)
  ),
  height = 5
)


## location(TRR): ICC versus HBM, separately within univariate and multivariate ----

plots_icc_vs_hbm <- enlist(responses)
for (response_i in seq_along(responses)) {
  plots_icc_vs_hbm[[response_i]] <- arrange_figure_comparison(
    data = uv_mv_icc_hbm[response == responses[[response_i]] & sum_fun %in% c("pointest", "map")],
    comparison_factor = model_nm,
    comparison_factor_order = models_comparison,
    colors_comparison = colors_models_comparison,
    margins = margins,
    title = titles[[response_i]],
    path_figs_results = path_figs_results,
    id_cols = c("region", "is_roi", "alpha_level", "q05"),
    color_col = "q05",
    color_lab = "d'",
    contrast_expr = "atanh(no_lscov_symm) - atanh(summarystat)",
    filename = paste0("trr_", responses[[response_i]])
  )
}

## precision(TRR): univariate versus multivariate ----

plot_trr_precision <- arrange_figure_comparison(
  data = uv_mv_icc_hbm[sum_fun == "sd"] %>% mutate(value = 1/log(value)),
  limits = NULL,
  comparison_factor = response,
  comparison_factor_order = responses,
  colors_comparison = colors_response,
  margins = margins,
  title = titles[[response_i]],
  path_figs_results = path_figs_results,
  id_cols = c("region", "is_roi", "alpha_level", "q05"),
  color_col = "q05",
  color_lab = "d'",
  contrast_expr = "rda - uv",
  filename = "trrprecision"
)

plot_trr_precision


## TRR thresholded: univariate versus multivariate ----

data <- full_join(
  uv_mv_icc_hbm %>% filter(sum_fun == "mean") %>% select(-sum_fun),
  uv_mv_icc_hbm %>%
    filter(sum_fun == "q05") %>%
    select(-sum_fun) %>%
    group_by(response) %>%
    mutate(
      lb = value,
      alpha_trr = highlight_overlay(value, upper = 0, lower = -1, lower_alpha = 0),
      value = NULL
    )
)

## check highlight_alpha range:
if (FALSE) {
  data %>%
    arrange(response, alpha_trr) %>%
    ggplot(aes(lb, alpha_trr)) +
    geom_point() +
    facet_wrap(~response)
}

p_trr_thresh_brains <-
  data %>%
  mutate(response = factor(response, levels = responses)) %>%
  group_by(response) %>%
  plot_surface(
    statistic_col = "value",
    #underlay = NULL,
    underlay = c(color = "grey", fill = "white"),
    border_size_values = c("TRUE" = 0.3, "FALSE" = 0.2),
    limits = c(-0.5, 1),
    breaks = c(0, 0.5, 1),
    alpha_col = "alpha_trr"
  ) +
  scale_alpha_identity(guide = TRUE) +
  labs(fill = "mean(TRR)") +
  facet_grid(
    vars(response),
    labeller = labeller(response = setNames(names(responses), responses)),
    switch = "y"
  )
p_trr_thresh_brains
ggsave(file.path(path_figs_results, "trr_thresh.pdf"), p_trr_thresh_brains, height = 2, width = 4.5)

# uv_mv_icc_hbm %>%
#   filter(sum_fun == "q05") %>%
#   pivot_wider(id_cols = c("region", "is_roi"), names_from = "response", values_from = "value") %>%
#   ggplot(aes(uv, rda)) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_point(aes(color = is_roi), size = 4) +
#   theme_default(base_size = 12)


## complementary benefits and example posteriors ----

## scatterplot panel

## change in SD(TRR) in multiv vs univ contrasts (x axis):
hbm_vs_icc <-
  uv_mv_icc_hbm[response == "rda" & sum_fun %in% c("map", "pointest")] %>%
  pivot_and_contrast(
    contrast = atanh(no_lscov_symm) - atanh(summarystat),
    id_cols = c("region", "is_roi"),
    names_from = "model_nm",
    values_from = "value",
    new_col = delta_loc
  )
## change in TRR in HBM versus SumStat estimates (y axis):
uv_vs_mv <-
  posterior_summaries$trr[sum_fun == "sd"] %>%
  pivot_and_contrast(
    contrast = log(uv) - log(rda),
    id_cols = "region",
    names_from = "response",
    values_from = "value",
    new_col = delta_sd
  )
d <- full_join(hbm_vs_icc, uv_vs_mv)

## contingincy table:
d_ctab <- d %>%
  mutate(sign_delta_loc = sign(delta_loc), sign_delta_sd = sign(delta_sd)) %>%
  group_by(sign_delta_loc, sign_delta_sd, is_roi) %>%
  summarise(n = n())
p_quadrants_scatter <-
  d %>%
  ggplot(aes(delta_sd, delta_loc)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)  +
  geom_point(shape = 21, size = 3) +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  scale_y_continuous(limits = c(-2, 2))

## posterior density panel

## plot example posteriors:
dir.create(file.path(path_figs_results_supp, "example_posteriors"), showWarnings = FALSE)
examples <- list(
  ul = d %>% filter(delta_loc > 1, delta_sd < 0.1) %>% pull(region),
  ur = d %>% filter(delta_loc > 1, delta_sd > 0.5) %>% pull(region),
  br = d %>% filter(delta_loc < -0.5, delta_sd > 0) %>% pull(region),
  bl = d %>% filter(delta_loc < -0.5, delta_sd < 0) %>% pull(region)
)
p_eg <- enlist(names(examples))
for (q in seq_along(examples)) {
  p_eg[[q]] <-
    posterior_samples$trr[region %in% examples[[q]]] %>%
    ggplot(aes(value, color = response, fill = response)) +
    geom_density(linewidth = 1.5, aes(fill = NULL, group = response)) +
    facet_wrap(vars(region), ncol = 6, scales = "free", labeller = as_labeller(label_regions)) +
    scale_color_manual(values = colors_response) +
    scale_fill_manual(values = colors_response) +
    theme(legend.position = "none") +
    labs(x = "TRR", y = "posterior density")
  ggsave(
    file.path(path_figs_results_supp, "example_posteriors", paste0("quadrant", q, ".pdf")),
    p_eg[[q]],
    dev = cairo_pdf,
    width = 8,
    height = ceiling(length(examples[[q]]) / 6) * 1.5
  )
}

eg <- c(
  "17Networks_LH_ContA_IPS_2",
  "17Networks_LH_ContA_PFCd_1",
  "17Networks_RH_LimbicB_OFC_3",
  "17Networks_RH_SomMotB_S2_8"
)

p_quadrants_density <-
  posterior_samples$trr[region %in% eg] %>%
  ggplot(aes(value, color = response, fill = response)) +
  geom_density(linewidth = 1.5, aes(group = response), alpha = 0.5) +
  geom_vline(
    data = summarystats$trr %>% filter(region %in% eg),
    aes(xintercept = r, color = response),
    linetype = "dashed", linewidth = 1
  ) +
  facet_wrap(vars(region), ncol = 2, labeller = as_labeller(label_regions_eg)) +
  scale_color_manual(values = colors_response) +
  scale_fill_manual(values = colors_response) +
  theme(
    legend.position = "none",
    #strip.text = element_text(size = 6),
    axis.line.y.left = element_blank()
  ) +
  labs(x = "TRR", y = "posterior density") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = 0)

p_quadrants <- p_quadrants_scatter + p_quadrants_density
ggsave(file.path(path_figs_results, "quadrants.pdf"), p_quadrants, dev = cairo_pdf, width = 130, height = 130/1.75, unit = "mm")


## association between SD(TRR) and reduction in trial-level noise ----


## xaxis: change in variability ratio, multivariate - univariate
response_ratio <-
  posterior_summaries$ratio[statistic == "ratio" & sum_fun == "map"] %>%
  pivot_and_contrast(
    contrast = ratio_uv - ratio_rda,
    id_cols = region,
    names_from = response,
    values_from = value,
    new_col = delta_ratio,
    names_prefix = "ratio_"
  )
## yaxis: change in sd(TRR), multivariate - univariate
response_sdtrr <-
  posterior_summaries$trr[sum_fun == "sd"] %>%
  pivot_and_contrast(
    contrast = log(sdtrr_uv) - log(sdtrr_rda),
    id_cols = region,
    names_from = response,
    values_from = value,
    new_col = delta_sd,
    names_prefix = "sdtrr_"
  )
d_ratio_vs_sdtrr <- full_join(response_ratio, response_sdtrr, by = c("region"))

p_ratio <-
  d_ratio_vs_sdtrr %>%
  ggplot(aes(delta_ratio, delta_sd)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")  +
  geom_abline(linetype = "dotted")  +
  stat_ellipse(type = "t", linetype = "dotted") +
  geom_point(shape = 21, size = 1.25) +
  labs(
    x = "\u0394Variab. Ratio, univ.\u2212multiv.\n(log ratio)",
    y = "\u0394SD(TRR), univ.\u2212multiv.\n(log ratio)"
  ) +
  theme(legend.position = "none") +
  geom_text(
    data = . %>%
      summarize(r = cor(delta_ratio, delta_sd), rho = cor(delta_ratio, delta_sd, method = "spearman")) %>%
      mutate(label = paste0("r = ", round(r, 2), "\n", "\n\u03C1 = ", round(rho, 2))),
    aes(x = -1, y = 1, label = label),
  ) +
  scale_y_continuous(limits = c(-1.4, 1.4)) +
  scale_x_continuous(limits = c(-1.3, 1.3))
## add axis text labels:
## relative certainty in TRR: increases with univariate -- increases with multivariate
## amount of trial vs. subject variability: less trial variability in mutlivariate -- less trial variability in univariate

ggsave(
  file.path(path_figs_results, "variability_ratio_vs_precision.pdf"),
  p_ratio,
  dev = cairo_pdf, width = 4.5, height = 4, scale = 0.8
  )



## weight vector analysis ----




