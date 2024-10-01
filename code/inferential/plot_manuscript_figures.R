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
library(knitr)
library(kableExtra)

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

noiseprojs <- readRDS(
  here(
    "out", "spatial",
    "noise_projs__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12_resample1_baseline.RDS"
  )
)


## minor computations and formatting for plotting ----

## get population stats and ROIs columns

popef <-
  posterior_summaries$popef %>%
  filter(response == "uv", sum_fun %in% c("mean", "sd", "map", "q05")) %>%
  dcast(... ~ sum_fun, value.var = "value") %>%
  mutate(
    is_roi = region %in% dmcc35_nms,
    tplus = mean / sd,
    is_roi_q05 = quantile(q05, 0.9) < q05,
    is_roi_mean = quantile(mean, 0.9) < mean,
    is_roi_map = quantile(map, 0.9) < map,
    is_roi_tplus = quantile(tplus, 0.9) < tplus,
    network = get_network(region)
  ) %>%
  rename(m = mean)  ## population-level univariate stroop contrast

## examine similarities:
if (FALSE) {
  sum(popef$is_roi & popef$is_roi_mean)  ## agreement (out of 40)
  sum(popef$is_roi & popef$is_roi_tplus)
  pairs(popef[, c("m", "map", "q05", "sd", "tplus")])
}

popef_network <- popef %>%
  group_by(network) %>%
  summarize(
    q05 = mean(q05),
    n_parcel = n(),
    n_roi = sum(is_roi),
    prop_roi = mean(is_roi)
  )

## format weight-vector data

noiseprojs[, proj_mv_relvar := abs(proj_mv_relvar)]
noiseprojs[, proj_uv_relvar := abs(proj_uv_relvar)]
noiseprojs[, proj_uv_absvar := abs(proj_uv_relvar)*var_total]
noiseprojs[, proj_mv_absvar := abs(proj_mv_relvar)*var_total]
noiseprojs[, var_dim_uv := proj_uv_absvar / proj_uv_scaled]  ## variance per dimension
noiseprojs[, var_dim_mv := proj_mv_absvar / proj_rda_scaled]
stopifnot(all.equal(noiseprojs$var_dim_uv, noiseprojs$var_dim_mv))  ## check for equality
noiseprojs <- noiseprojs %>%
  rename(var_dim = var_dim_mv, region = roi) %>%
  mutate(var_dim_uv = NULL)


## reorder popef

popef <- popef %>%
  mutate(network = factor(network, levels = popef_network$network[order(popef_network$q05)]))
popef_network <- popef_network %>%
  mutate(network = factor(network, levels = popef_network$network[order(popef_network$q05)]))

uv_mv_icc_hbm <-
  bind_rows(
    summarystats$trr %>% rename(value = r) %>% mutate(sum_fun = "pointest", statistic = "trr"),
    posterior_summaries$trr
  ) %>%
  full_join(popef %>% select(region, tplus, q05, is_roi), by = "region")

## change in SD(TRR) in multiv vs univ contrasts

hbm_vs_icc <-
  uv_mv_icc_hbm[response == "rda" & sum_fun %in% c("map", "pointest")] %>%
  pivot_and_contrast(
    contrast_expr = "atanh(no_lscov_symm) - atanh(summarystat)",
    id_cols = c("region", "is_roi", "tplus"),
    names_from = "model_nm",
    values_from = "value",
    new_col = "delta_loc"
  )

## change in TRR in HBM versus SumStat estimates

uv_vs_mv <-
  posterior_summaries$trr[sum_fun == "sd"] %>%
  pivot_and_contrast(
    contrast = "log(uv) - log(rda)",
    id_cols = "region",
    names_from = "response",
    values_from = "value",
    new_col = "delta_sd"
  )
d <- full_join(hbm_vs_icc, uv_vs_mv)

## change in variability ratio, multivariate - univariate

response_ratio <-
  posterior_summaries$ratio[statistic == "ratio" & sum_fun == "map"] %>%
  pivot_and_contrast(
    contrast = "ratio_uv - ratio_rda",
    id_cols = "region",
    names_from = "response",
    values_from = "value",
    new_col = "delta_ratio",
    names_prefix = "ratio_"
  )
d_ratio_vs_sdtrr <- full_join(response_ratio, d, by = c("region"))



## population-level contrast on univariate response ----

## plot surfaces

p_pop_brain_unthresh <-
  popef %>%
  plot_surface(
    statistic = "tplus",
    underlay = NULL,
    alpha_col = NULL
  ) +
  labs(fill = univ_stat_lab) +
  theme(legend.position = c(0.5, 0))
ggsave(file.path(path_figs_results, "pop_tplus_unthresh.pdf"), p_pop_brain_unthresh, height = 2, width = 4.5)
#ggsave(file.path(path_figs_results_supp, "pop_tplus_unthresh.pdf"), p_pop_brain_unthresh, height = 2, width = 4.5)

## check stats within network assignments

thresh_roi <- unname(quantile(popef$q05, 0.9))

p_parcel_q05 <-
  popef %>%
  ggplot(aes(q05, fill = is_roi)) +
  geom_histogram(bins = 30) +
  #geom_vline(xintercept = quantile(popef$q05, 0.9)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Population-level univariate\nStroop contrast (5%ile)", y = "Number of parcels") +
  annotate(
    geom = "text", x = thresh_roi + 0.005, y = 30, label = "'ROI'", color = colors_roi[["TRUE"]], size = 3,
    hjust = 0
  ) +
  scale_fill_manual(values = colors_roi) +
  scale_x_continuous(breaks = c(-0.1, 0.1, round(thresh_roi, 2))) +
  theme(legend.position = "none")

p_network_q05 <-
  popef %>%
  ggplot(aes(q05, network)) +
  #geom_vline(xintercept = thresh_roi) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "grey", color = "black", width = 0.5) +
  theme(
    axis.line.y.left = element_blank(),
    panel.grid.major.y = element_line(linewidth = rel(0.25), linetype = "dashed"),
    axis.text.y = element_text(size = rel(0.75))
  ) +
  scale_x_continuous(breaks = c(-0.1, 0.1, round(thresh_roi, 2))) +
  labs(
    y = "Cortical Network\n(Schaefer 400-17)",
    x = "Population-level univariate\nStroop contrast (5%ile)"
  )

p_network_props <-
  popef_network %>%
  mutate(is_roi = TRUE) %>%
  ggplot(aes(prop_roi, network, fill = is_roi)) +
  geom_bar(color = "black", width = 0.6, stat = "identity") +
  geom_text(
    aes(label = paste0(n_roi, "/", n_parcel), y = network, x = prop_roi + 0.004),
    hjust = 0,
    color = "black",
    size = 1.75
  ) +
  scale_fill_manual(values = colors_roi) +
  scale_x_continuous(limits = c(0, 0.55), breaks = c(0, 0.25, 0.5)) +
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

p_uv <- arrange_plots(
  p_pop_brain_unthresh,
  p_parcel_q05 + p_network_q05 + p_network_props,
  filename = "supp/pop_uv",
  path = path_figs_results,
  plot_layout = c(
    area(t = 1, b = 62, l = 1, r = 220),
    area(t = 63, b = 100, l = 21, r = 200)
  ),
  width = two_column_width,
  height = two_column_width/1.9
)


## comparison plots ----
## see params_comparison_plots list in _parameters_viz.R

plots_comparison <- enlist(names(params_comparison_plots))
for (plot_i in seq_along(params_comparison_plots)) {
  pars <- params_comparison_plots[[plot_i]]
  ## subset
  data <- uv_mv_icc_hbm[eval(pars$i_expr), ]
  ## modify
  if (!is.null(pars$j_expr)) data[, eval(pars$j_expr)]
  ## plot
  arrange_args <- c(list(data = data), pars[!names(pars) %in% c("i_expr", "j_expr")])
  plots_comparison[[plot_i]] <- do.call(arrange_figure_comparison, arrange_args)
}


## TRR thresholded: univariate versus multivariate ----

## set params:

## parcels with trr > upper: full opacity
## parcels with trr < upper but > lower: linear function btw 0.1 and 1
## parcels with trr < lower: 0.1
highlight_pars <- list(
  upper = 0, lower = -1, lower_alpha = 0.2,
  scaling_fun = function(x) x^2
)
color_scale <- scale_fill_viridis_c(
  option = "magma", na.value = "white", oob = scales::squish, limits = c(0, 1), breaks = c(0, 0.5, 1)
)

data <- full_join(
  uv_mv_icc_hbm %>% filter(sum_fun == "mean") %>% select(-sum_fun),
  uv_mv_icc_hbm %>%
    filter(sum_fun == "q05") %>%
    select(-sum_fun) %>%
    group_by(response) %>%
    mutate(
      trr_q05 = value,
      alpha_trr = highlight_overlay(
        value, upper = highlight_pars$upper, lower = highlight_pars$lower,
        lower_alpha = highlight_pars$lower_alpha,
        scaling_fun = highlight_pars$scaling_fun
        ),
      value = NULL
    )
)

## check highlight_alpha range:
if (FALSE) {
  data %>%
    arrange(response, alpha_trr) %>%
    ggplot(aes(trr_q05, alpha_trr)) +
    geom_point() +
    facet_wrap(~response)
}

wrapped_labs <- stringr::str_wrap(comparison_factor_labs$uv_mv, width = 7)
names(wrapped_labs) <- names(comparison_factor_labs$uv_mv)
p_trr_thresh_brains <- data %>%
  plot_surface(
    statistic_col = "value",
    underlay = c(color = "grey50", fill = "grey50"),
    #border_size_values = c(`TRUE` = 0.5, `FALSE` = 0.1),
    alpha_col = "alpha_trr",
    facet_factor = "response",
    facet_factor_order = wrapped_labs,
    fill_label = "Mean(TRR) (r)",
    guides_ = guides(
      fill = guide_colorbar(title.position = "top", title.vjust = 0.9),
      color = "none", size = "none"
    ),
    scale_fill = function(...) color_scale
  ) +
  theme(
    legend.position = c(0.25, -0.1),
    plot.caption.position = "plot"
  ) +
  labs(
    caption = "fully opaque parcels have 5%ile(TRR) > 0   "
  )

colorgrid <- create_colorgrid(
  colorscale = color_scale,
  limits_color_stat = color_scale$limits,
  limits_alpha_stat = c(-1, 1),
  upper = highlight_pars$upper,
  lower = highlight_pars$lower,
  lower_alpha = highlight_pars$lower_alpha,
  scaling_fun = highlight_pars$scaling_fun,
  n = 101
)
colorscale <- colorgrid %>%
  filter(alpha_stat <= 0.5) %>%
  ggplot(aes(as.factor(color_stat), as.factor(alpha_stat))) +
  geom_raster(fill = "grey50") +
  geom_raster(aes(fill = color, alpha = alpha)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  scale_y_discrete(breaks = c(-1, -0.5, 0, 0.5), labels = c("-1", "-0.5", "0", ">0")) +
  scale_x_discrete(breaks = c(0, 0.5, 1), labels = c("<=0", "0.5", "1")) +
  labs(x = "Mean TRR (r)", y = "5%ile TRR (r)") +
  theme(
    axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank(), axis.title = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.ticks = element_line(linewidth = 1/4, color = "grey40"),
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  )
ggsave(
  file.path(path_figs_results, "trr_thresh_new_colorscale.pdf"),
  colorscale,
  width = 20,
  height = 15,
  units = "mm"
)

ggsave(
  file.path(path_figs_results, "trr_thresh_new.pdf"),
  p_trr_thresh_brains,
  width = oneandhalf_column_width,
  height = 55,
  units = "mm"
)



## complementary benefits and example posteriors ----

## scatterplot panel

## contingincy table:
d_ctab <- d %>%
  mutate(sign_delta_loc = sign(delta_loc), sign_delta_sd = sign(delta_sd)) %>%
  group_by(sign_delta_loc, sign_delta_sd, is_roi) %>%
  summarise(n = n()) %>%
  mutate(
    quadrant = case_when(
      sign_delta_sd == -1 & sign_delta_loc == 1 ~ "q1",
      sign_delta_sd == 1  & sign_delta_loc == 1 ~ "q2",
      sign_delta_sd == -1 & sign_delta_loc == -1 ~ "q3",
      sign_delta_sd == 1 & sign_delta_loc == -1 ~ "q4"
    )
  )


## posterior density panel

## plot example posteriors:
dir.create(file.path(path_figs_results_supp, "example_posteriors"), showWarnings = FALSE)
examples <- list(
  q1 = d %>% filter(delta_sd < 0.1, delta_loc > 1) %>% pull(region),
  q2 = d %>% filter(delta_sd > 0.5, delta_loc > 1) %>% pull(region),
  q3 = d %>% filter(delta_sd < 0, delta_loc < -0.5) %>% pull(region),
  q4 = d %>% filter(delta_sd > 0, delta_loc < -0.5) %>% pull(region)
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
  "17Networks_RH_ContA_PFCl_3",
  "17Networks_RH_LimbicB_OFC_3",
  "17Networks_RH_SomMotB_S2_8"
)

if (FALSE) {
  ## have a look at anatomical locations:
  data_eg <- uv_mv_icc_hbm %>% filter(sum_fun == "mean", response == "rda") %>% mutate(value = ifelse(region %in% eg, 1, NA))
  plot_surface(
    data = data_eg,
    statistic = "value",
    underlay = NULL,
    alpha_col = NULL
  )
}


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
  labs(x = "Test\U2012retest reliability (r)", y = "posterior density") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = 0) +
  geom_text(
    data = data.frame(
      value = -0.6,
      response = c("uv", "rda"),
      label = c("univar.", "\n\nmultivar."),
      region = eg[1]
    ),
    aes(label = label, y = 2),
    size = 3
  )


## scatter plot:

p_quadrants_scatter <-
  d %>% 
  mutate(is_example_region = region %in% eg) %>%
  arrange(is_example_region, tplus) %>%
  ggplot(aes(delta_sd, delta_loc)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")  +
  geom_point(
    aes(color = tplus, shape = is_example_region, alpha = is_example_region, stroke = is_example_region),
    size = 1.5
  ) +
  geom_text(
    data = d_ctab %>% 
      group_by(quadrant) %>%
      summarize(n = sum(n)) %>%
      mutate(
        delta_sd = case_when(
          quadrant == "q1" ~ -0.6,
          quadrant == "q2" ~ 0.55,
          quadrant == "q3" ~ -0.6,
          quadrant == "q4" ~ 1
        ),
        delta_loc = case_when(
          quadrant == "q1" ~ 2,
          quadrant == "q2" ~ 2,
          quadrant == "q3" ~ -0.5,
          quadrant == "q4" ~ -0.5
        ),
    ),
    aes(label = paste0(quadrant, ", N = ", n)),
    size = 3,
    color = "grey50"
  ) +
  scale_color_continuous_diverging("Blue-Red 3", breaks = c(-6, 0, 6)) +
  scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16)) +
  scale_alpha_manual(values = c("TRUE" = 4/5, "FALSE" = 3/4)) +
  labs(
    x = "\U0394Precision(TRR): multivar. \U2212 univar.\n(log)",
    y = "\U0394Multivar. TRR: MAP \U2212 ICC\n(z difference)",
    color = bquote("t"^"+")
  ) +
  guides(shape = "none", alpha = "none", color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme(
    legend.position = c(0.85, 0.1),
    legend.direction = "horizontal",
    legend.key.height = unit(1 / 8, "cm"),
    legend.key.width = unit(1 / 4, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )


p_quadrants <- p_quadrants_scatter + p_quadrants_density
ggsave(
  file.path(path_figs_results, "quadrants.pdf"),
  p_quadrants,
  dev = cairo_pdf, width = 130, height = 130 / 1.75, unit = "mm"
)


## weight vector analysis ----

## example ROI

region_noiseproj <- "17Networks_RH_ContA_PFCl_3"
noiseprojs_eg <- noiseprojs[region == region_noiseproj]

p_alignment_eg <-
  noiseprojs_eg %>%
  melt(measure.vars = c("proj_uv_scaled", "proj_rda_scaled")) %>%
  group_by(subj, region, dimension, variable) %>%
  summarize(
    cosine_sim = tanh(mean(atanh(value))),
    sd_dim = mean(sqrt(var_dim))
  ) %>% 
  ggplot(aes(x = dimension)) +
  stat_summary(aes(y = sd_dim), geom = "point", fun = "mean", size = 1/3) +
  stat_summary(aes(y = sd_dim), geom = "line", fun = "mean", linewidth = 1/3) +
  stat_summary(aes(y = sd_dim), geom = "errorbar", fun.data = "mean_cl_boot", width = 0, linewidth = 1/3) +
  stat_summary(
    aes(y = cosine_sim * weight_vector_coeff, group = variable, color = variable),
    geom = "point", fun = "mean", size = 1/3
  ) +
  stat_summary(
    aes(y = cosine_sim * weight_vector_coeff, group = variable, color = variable),
    geom = "line", fun = "mean", linewidth = 1/3
  ) +
  stat_summary(
    aes(y = cosine_sim * weight_vector_coeff, group = variable, color = variable),
    geom = "errorbar", fun.data = "mean_cl_boot", width = 0, linewidth = 1/3
  ) +
  scale_y_continuous(
    name = "SD(trial-level noise)",
    sec.axis = sec_axis(~ . / weight_vector_coeff, name = "Noise alignment\n(cosine similarity)")
  ) +
  scale_x_continuous(trans = "log", breaks = c(1, 2, 4, 8, 16, 32, 58)) +
  scale_color_manual(values = color_cosine_sim) +
  annotate(
    geom = "text", x = 2.25, y = 3, label = expression("w"["univar."]), color = color_cosine_sim["proj_uv_scaled"]
  ) +
  annotate(
    geom = "text", x = 16, y = 1, label = expression('w'['multivar.']), color = color_cosine_sim["proj_rda_scaled"]
  ) +
  theme(
    axis.line.y.right = element_line(color = color_cosine_sim["axis"]),
    axis.title.y.right = element_text(color = color_cosine_sim["axis"]),
    axis.text.y.right = element_text(color = color_cosine_sim["axis"]),
    axis.ticks.y.right = element_line(color = color_cosine_sim["axis"]),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Trial-level noise\ncomponent (log scale)",
  title = label_regions(region_noiseproj)
  )

## stats on all parcels*subjects
p_totalnoise <-
  noiseprojs[, c("region", "subj", "wave", "dimension", "proj_uv_absvar", "proj_mv_absvar")] %>%
  ## get total SD (sum SDs over dimensions):
  group_by(subj, region, wave) %>%
  summarize(
    proj_mv_absvar = mean(sqrt(proj_mv_absvar)),
    proj_uv_absvar = mean(sqrt(proj_uv_absvar))
  ) %>%
  ## convert to units SD then get log ratio: univariate - multivariate total var:
  mutate(delta_absvar = log(proj_uv_absvar / proj_mv_absvar)) %>%
  ## average over waves:
  group_by(subj, region) %>%
  summarize(delta_absvar = mean(delta_absvar)) %>%
  ggplot(aes(delta_absvar)) +
  geom_histogram(position = "identity", alpha = 0.5, color = "black", bins = 20) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  labs(
    x = "\U0394Total aligned trial-level variab.:\nunivar. \U2212 multivar. (log)",
    y = "Number of\nparcels*subjects",
    title = "All parcels*subjects"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_alignment <- p_alignment_eg / p_totalnoise

## association between SD(TRR) and reduction in trial-level noise ----

p_ratio <-
  d_ratio_vs_sdtrr %>%
  arrange(tplus) %>%
  ggplot(aes(delta_ratio, delta_sd)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")  +
  geom_abline(linetype = "dotted", color = "grey50")  +
  stat_ellipse(type = "t", linetype = "dotted", color = "grey50") +
  geom_point(aes(color = tplus), shape = 16, size = 2, stroke = 0, alpha = 1) +
  labs(
    x = "\u0394Variab. Ratio (log)",
    y = "\U0394SD(TRR)  (log)",
    color = bquote("t"^"+")
  ) +
  theme(legend.position = "none") +
  geom_text(
    data = . %>%
      summarize(r = cor(delta_ratio, delta_sd), rho = cor(delta_ratio, delta_sd, method = "spearman")) %>%
      mutate(label = paste0("r = ", round(r, 2), "\n", "\u03C1 = ", round(rho, 2))),
    aes(x = 1.2, y = -0.3, label = label),
    size = 3
  ) +
  scale_y_continuous(limits = c(-1.4, 1.4)) +
  scale_x_continuous(limits = c(-1.3, 1.3)) +
  annotate(
    geom = "text", x = 0, y = -1.2, size = 3,
    color = "grey50", fontface = "bold",
    label = expression(atop(bold(frac(trial,subject)), bold("variability")))
  ) +
  annotate(geom = "text", x = 1, y = -1.2, size = 3, label = "higher in\nunivariate", color = "grey50") +
  annotate(geom = "text", x = -1, y = -1.2, size = 3, label = "higher in\nmultivariate", color = "grey50") +
  annotate(
    geom = "text", x = -1, y = 0, size = 3, label = "TRR\nvariability", angle = 90, color = "grey50", fontface = "bold"
  ) +
  annotate(geom = "text", x = -1, y = 1.2, size = 3, label = "higher in\nunivariate", color = "grey50") +
  annotate(
    geom = "segment", x = -0.6, xend = 0.6, y = -1.2, yend = -1.2,
    arrow = arrow(ends = "both", type = "closed", length = unit(1, "mm")), color = "grey50"
  ) +
  annotate(
    geom = "segment", x = -1, xend = -1, y = 0.75, yend = -0.75,
    arrow = arrow(ends = "both", type = "closed", length = unit(1, "mm")), color = "grey50"
  ) +
  scale_color_continuous_diverging("Blue-Red 3", breaks = c(-6, 0, 6)) +
  guides(
    shape = "none", alpha = "none", color = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  theme(
    legend.position = c(0.375, 0.9),
    legend.direction = "horizontal",
    legend.key.height = unit(1 / 8, "cm"),
    legend.key.width = unit(1 / 4, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.margin = unit(c(6, 1, 6, 1), "mm")
  )

ggsave(
  file.path(path_figs_results, "variability_ratio_vs_precision.pdf"),
  p_ratio,
  dev = cairo_pdf, width = one_column_width, height = one_column_width*0.95,
  units = "mm"
)

## horizontal layout:

ggsave(
  file.path(path_figs_results, "noise_fig.pdf"),
  wrap_elements(p_alignment) + wrap_elements(p_ratio) + plot_layout(widths = c(1, 1.2)),
  dev = cairo_pdf, width = two_column_width, height = one_column_width*1.25,
  units = "mm"
)



## tables ----

## all parcels + stats table

table_all_stats <-
  posterior_summaries$trr[sum_fun %in% c("pointest", "map", "mean", "q05", "sd")] %>%
  bind_rows(summarystats$trr %>% rename(value = r) %>% mutate(sum_fun = "pointest", statistic = "trr")) %>%
  pivot_wider(
    id_cols = "region", names_from = c(response, sum_fun), values_from = "value", names_prefix = "trr_"
  ) %>%
  full_join(popef %>% select(popef_q05 = q05, popef_tplus = tplus, is_roi, region, network), by = "region") %>%
  mutate(is_dmcc35 = region %in% rois[dmcc35])
fwrite(table_all_stats, file.path(path_figs_results, "all_stats.csv"))

table_all_stats %>% filter(trr_rda_q05 > 0) %>% nrow
table_all_stats %>% filter(trr_uv_q05 > 0) %>% nrow
table_all_stats %>% filter(trr_rda_q05 > 0 & is_roi) %>% nrow


top_trr_rda <- table_all_stats %>% top_n(40, trr_rda_q05) %>% pull(region)
top_trr_uv <- table_all_stats %>% top_n(40, trr_uv_q05) %>% pull(region)

table_trr <-
  table_all_stats %>%
  mutate(
    region_lab = label_regions(region),
    region_lab = paste0(region_lab, case_when(is_dmcc35 ~ "*", TRUE ~ ""))
  ) %>%
  filter(region %in% top_trr_rda) %>%
  arrange(-trr_rda_q05) %>%
  select(
    region_lab, popef_tplus,
    trr_uv_map, trr_uv_q05, trr_uv_pointest,
    trr_rda_map, trr_rda_q05, trr_rda_pointest
  ) %>%
  mutate(across(where(is.numeric), function(x) round(x, 2)))

names(table_trr) <- c(
  "Parcel (Schaefer 400-17)",
  "$t^+$",
  "MAP",
  "5\\%ile",
  "ICC",
  "MAP",
  "5\\%ile",
  "ICC"
)
header <- c(" " = 2, "Univariate TRR" = 3, "Multivariate TRR" = 3)
latex_table <-
  table_trr %>%
  kable("latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(header) %>%
  ## remove first two and last lines, which define table environment
  ## do this in latex, so caption can be written there as well.
  strsplit("\n") %>%
  unlist %>%
  .[-c(1, 2, length(.))] %>%
  paste0(collapse = "\n")
writeLines(latex_table, file.path(path_figs_results, "table_trr.tex"))
