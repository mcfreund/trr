library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(data.table)
library(brms)
library(posterior)
library(ggplot2)
library(colorspace)
library(mfutils)
library(ggsegSchaefer)
library(patchwork)

source(here("code", "inferential", "_parameters_viz.R"))
source(here("code", "inferential", "_parameters_reliability.R"))
source(here("code", "inferential", "_utils_viz.R"))
source(here("code", "inferential", "_utils_reliability.R"))
source(here("code", "_atlas.R"))
source(here("code", "_constants.R"))
source(here("code", "_paths.R"))

## other constants, paths, etc...

session <- "baseline"
n_core <- 4
plan(multicore, workers = n_core)

dir.create(file.path(path_figs, "model_comparison", "supp"), showWarnings = FALSE, recursive = FALSE)

## read and subset data

posterior_samples <- load_summaries(
  data_type = "posterior_samples",
  prefixes = c("popef", "ratio", "trr"),
  session = session,
  regions = core32_nms,
  sum_fun = subset_and_order
)

posterior_summaries <- load_summaries(
  data_type = "posterior_summaries",
  prefixes = c("popef", "ratio", "trr"),
  sess = session,
  regions = core32_nms,
  sum_fun = subset_and_order
)
stopifnot(any(stat_names %in% unique(posterior_summaries$trr$sum_fun)))  ## check

summarystats <- load_summaries(
  data_type = "summarystat",
  prefixes = c("popef", "trr"),
  session = session,
  regions = core32_nms,
  sum_fun = subset_and_order,
  models_order = NULL
)
criteria <- subset_and_order(
  fread(file.path(path_reliab, "fit_criteria.csv")),
  regions = core32_nms,
  session = "baseline"
) ## model comparison info


# ELPD-LOO model comparison ----

# For rda and uv respectively
p_elpd <-
  criteria %>%
  filter(Term == "elpd_loo") %>%
  ggplot(aes(model_nm, Estimate)) +
  geom_line(aes(group = region)) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.1) +
  scale_y_continuous(
    trans = arsinh,
    breaks = c(0, -3000, -5000, -10000, -15000, -20000),
    labels = scales::scientific) +
  scale_x_discrete(labels = models_hbm_vals2nms) +
  theme(legend.position = "none") +
  labs(x = "HBM Model", y = "ELPD LOO") +
  facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate"))) +
  theme(axis.text.x = element_text(size = 7))

# Differences between reduced models and full model, excluding fixed_sigma
elpd_diff <- criteria %>%
  filter(Term == "elpd_loo") %>%
  pivot_wider(names_from = "model_nm", values_from = c("Estimate", "SE")) %>%
  mutate(
    `Homog. - Full`  = Estimate_fixed_sigma - Estimate_full,
    `ILS - Full`     = Estimate_no_lscov - Estimate_full,
    `ILS Sym - Full` = Estimate_no_lscov_symm - Estimate_full
  ) %>%
  select(!matches("^Estimate|^SE"))

# Scatterplot of differences
p_elpd_diff <- elpd_diff %>%
  ggplot(aes(x = `ILS Sym - Full`, y = `ILS - Full`)) +
  geom_point(shape = 21, color = "white", fill = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "none") +
  facet_wrap(vars(response), labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")))

ggsave(file.path(path_figs, "model_comparison", "elpd.pdf"), p_elpd + p_elpd_diff, width = 8, height = 2.5)



# Population-level Stroop effects ----

## format for plotting (to wide):
popef_w <- dcast(posterior_summaries$popef, ... ~ sum_fun, value.var = "value")

## Mean of population-level Stroop effect

# Bar plot
p_pop_bars_mu <- enlist(responses)
for (res in responses) {
  res_lab <- switch(res, uv = "Univariate Response", rda = "Multivariate Response")
  p_pop_bars_mu[[res]] <-
    popef_w[response == res] %>%
    ggplot(aes(model_nm, map, fill = model_nm)) +
    geom_hline(yintercept = 0) +
    geom_col(width = 0.5, color = "black") +
    geom_errorbar(aes(ymin = hdi95_lower, ymax = hdi95_upper), width = 0.1, color = "black") +
    scale_fill_manual(
      values = colors_models,
      guide = guide_legend(title = "models", title.position = "top", direction = "horizontal"),
      labels = names(models)
    ) +
    theme(legend.position = c(0.7, 0.07), legend.direction = "horizontal", axis.text.x = element_blank()) +
    facet_wrap(vars(region), labeller = as_labeller(label_regions)) +
    labs(title = res_lab, x = element_blank(), y = "Mean(Pop.-Level Stroop Effect) + 95% CI", fill = "models") +
    theme(strip.text = element_text(size = 6), legend.text = element_text(size = 6))
}

# Line plot
p_pop_box_mu <-
  popef_w %>%
  ggplot(aes(model_nm, map)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
  geom_line(aes(group = region), alpha = 0.5, color = "dodgerblue") +
  scale_color_viridis_d() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = models_hbm_vals2nms) +
  labs(x = "Model", y = "MAP(Pop.-Level Stroop Effect)") +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  facet_wrap(
    vars(response),
    labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")),
    scales = "free_y"
  )


## Standard error of population-level Stroop effect

# Bar plot
p_pop_bars_se <- enlist(responses)
for (res in responses) {
  res_lab <- switch(res, uv = "Univariate Response", rda = "Multivariate Response")
  p_pop_bars_se[[res]] <-
    popef_w[response == res] %>%
    ggplot(aes(model_nm, sd, fill = model_nm)) +
    geom_hline(yintercept = 0) +
    geom_col(width = 0.5, color = "black") +
    scale_fill_manual(
      values = colors_models,
      guide = guide_legend(title = "models", title.position = "top", direction = "horizontal"),
      labels = names(models)
    ) +
    theme(legend.position = c(0.7, 0.07), legend.direction = "horizontal", axis.text.x = element_blank()) +
    facet_wrap(vars(region), labeller = as_labeller(label_regions)) +
    labs(title = res_lab, x = element_blank(), y = "SE(Pop.-Level Stroop Effect)", fill = "models") +
    theme(strip.text = element_text(size = 6), legend.text = element_text(size = 7))
}

# Line plot
p_pop_box_se <-
  popef_w %>%
  ggplot(aes(model_nm, sd)) +
  geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
  geom_line(aes(group = region), alpha = 0.5, color = "dodgerblue") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = models_hbm_vals2nms) +
  labs(x = "Model", y = "SE(Pop.-Level Stroop Effect)") +
  theme(legend.position = "none", axis.text.x = element_text(size = 7)) +
  facet_wrap(
    vars(response),
    labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")),
    scales = "free_y"
  )

## Arrange and save
p_pop_bars <-
  (p_pop_bars_mu$uv + theme(legend.position = "none") + p_pop_bars_mu$rda + theme(legend.position = "none")) /
  (p_pop_bars_se$uv + theme(legend.position = "none") + p_pop_bars_se$rda)
p_pop_box <- p_pop_box_mu / p_pop_box_se
ggsave(file.path(path_figs, "model_comparison", "stroop_pop_bars.pdf"), p_pop_bars, width = 12, height = 12)
ggsave(file.path(path_figs, "model_comparison", "stroop_pop_box.pdf"), p_pop_box, width = 7, height = 6)



# Test-retest reliability ----

## Per-parcel posterior distributions

p_trr_dens <-
  posterior_samples$trr %>%
  ggplot(aes(value, color = model_nm, linetype = response)) +
  geom_density(linewidth = 0.5) +
  geom_vline(data = summarystats$trr, aes(xintercept = r, linetype = response), linewidth = 0.3) +
  facet_wrap(vars(region), labeller = as_labeller(label_regions)) +
  scale_color_manual(values = colors_models[models_hbm], labels = names(models), name = NULL) +
  scale_linetype_manual(values = linetype_response, labels = response_labels, name = NULL) +
  theme(
    legend.position = c(0.7, 0.05),
    legend.direction = "horizontal", legend.box = "vertical",
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  labs(
    title = "Posterior densities of TRR(Stroop Effect)",
    x = NULL,
    caption = "vertical lines indicate summary-stat TRR"
  )

ggsave(file.path(path_figs, "model_comparison", "stroop_trr_dens.pdf"), p_trr_dens, width = 7, height = 7)

## Statistics of posterior distributions

p_trr_stats <- enlist(unique(posterior_summaries$trr$sum_fun))
for (i in seq_along(p_trr_stats)) {
  stat_name <- stat_names[i]
  stat_label <- paste0(names(stat_name), "(Variab. Ratio)")
  stat_label <- paste0(names(stat_names)[i], "(TRR)")
  f1 <-
    posterior_summaries$trr[sum_fun == stat_name] %>%
    ggplot(aes(model_nm, value, color = region)) +
    geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
    geom_line(aes(group = region), alpha = 0.5, color = "dodgerblue") +
    scale_x_discrete(labels = models_hbm_vals2nms) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 7)) +
    labs(y = stat_label) +
    facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")), scales = "free_y")
  f2 <-
    posterior_summaries$trr[sum_fun == stat_name] %>%
    pivot_wider(names_from = response, values_from = value) %>%
    rename(model = model_nm, roi = region) %>%
    ggplot(aes(rda, uv)) +
    geom_point(shape = 21, color = "white", fill = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_viridis_d(guide = "none") +
    labs(x = paste0("Multivariate ", stat_label), y = paste0("Univariate\n", stat_label)) +
    facet_grid(cols = vars(model), labeller = labeller(model = setNames(names(models), models)))
  p_trr_stats[[i]] <- list(f1, f2)
}

pwalk(
  list(p_trr_stats, stat_names),
  \(l, nm) {
    ggsave(
      file.path(path_figs, "model_comparison", "supp", paste0("stroop_trr_", nm, ".pdf")),
      l[[1]] / l[[2]],
      width = 6, height = 5
    )
  }
)


# Variability ratio log(sd(trial)/sd(subj)) ----

posterior_samples$ratio[, model_nm := factor(model_nm, levels = models_hbm, ordered = TRUE)]

## Per-parcel posterior distributions

p_ratio_dens <-
  posterior_samples$ratio[statistic == "ratio"] %>%
  ggplot(aes(value, color = model_nm, linetype = response)) +
  geom_density(linewidth = 0.5) +
  facet_wrap(vars(region), labeller = as_labeller(label_regions)) +
  coord_cartesian(xlim = c(1, 6)) +
  scale_color_manual(values = colors_models[models_hbm], labels = names(models), name = NULL) +
  scale_linetype_manual(values = linetype_response, labels = response_labels, name = NULL) +
  theme(
    legend.position = c(0.7, 0.05),
    legend.direction = "horizontal", legend.box = "vertical",
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6)
  ) +
  labs(
    title = "Posterior densities of SD(trial)/SD(subj)",
    x = NULL,
    caption = "log(ratio)"
  )

ggsave(file.path(path_figs, "model_comparison", "stroop_ratio_dens.pdf"), p_ratio_dens, width = 7, height = 7)


## Statistics of posterior distributions

p_ratio_stats <- enlist(stat_names)
for (i in seq_along(stat_names)) {
  stat_name <- stat_names[i]
  stat_label <- paste0(names(stat_name), "(Variab. Ratio)")
  f1 <-
    posterior_summaries$ratio[sum_fun == stat_name & statistic == "ratio"] %>%
    ggplot(aes(model_nm, value, color = region)) +
    geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
    geom_line(aes(group = region), alpha = 0.5, color = "dodgerblue") +
    scale_x_discrete(labels = models_hbm_vals2nms) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 7)) +
    labs(y = stat_label) +
    facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")), scales = "free_y")
  f2 <-
    posterior_summaries$ratio[sum_fun == stat_name & statistic == "ratio"] %>%
    pivot_wider(names_from = response, values_from = value) %>%
    rename(model = model_nm, roi = region) %>%
    ggplot(aes(rda, uv)) +
    geom_point(shape = 21, color = "white", fill = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_viridis_d(guide = "none") +
    labs(x = paste0("Multivariate ", stat_label), y = paste0("Univariate\n", stat_label)) +
    facet_grid(cols = vars(model), labeller = labeller(model = setNames(names(models), models)))
  p_ratio_stats[[i]] <- list(f1, f2)
}

pwalk(
  list(p_ratio_stats, stat_names),
  \(l, nm) {
    ggsave(
      file.path(path_figs, "model_comparison", "supp", paste0("stroop_ratio_", nm, ".pdf")),
      l[[1]] / l[[2]],
      width = 6, height = 5
    )
  }
)

## arrange ----
lims_trr <- list(
  x = scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)),
  y = scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1))
)
lims_ratio <- list(
  x = scale_x_continuous(limits = c(1, 3.5), breaks = c(1, 2, 3, 4)),
  y = scale_y_continuous(limits = c(1, 3.5), breaks = c(1, 2, 3, 4))
)

p_stats <-
  (p_trr_stats$map[[1]] + lims_trr$y              + p_trr_stats$map[[2]] + lims_trr) /
  (p_trr_stats$q05[[1]] + lims_trr$y      + p_trr_stats$q05[[2]] + lims_trr) /
  (p_ratio_stats$map[[1]] + lims_ratio$y          + p_ratio_stats$map[[2]] + lims_ratio) /
  (p_ratio_stats$q05[[1]] + lims_ratio$y  + p_ratio_stats$q05[[2]] + lims_ratio)
ggsave(file.path(path_figs, "model_comparison", paste0("stroop_trratio_stats.pdf")), p_stats, width = 8, height = 7)
