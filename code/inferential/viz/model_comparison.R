library(here)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
library(brms)
library(posterior)
library(colorspace)
library(knitr)
library(mfutils)
library(ggsegSchaefer)
library(purrr)
library(furrr)
library(patchwork)

source(here("code", "_funs.R"))
source(here("code", "_constants.R"))
source(here("code", "inferential", "_plotting.R"))

## plot settings

theme_set(theme_minimal(base_size = 10))
theme_update(
  strip.background = element_rect(fill = "transparent", color = "transparent"),
  axis.line.y.left = element_line(),
  axis.line.x.bottom = element_line(),
  axis.ticks = element_line(),
  panel.grid = element_blank()
)

## other constants, paths, etc...

session <- "baseline"
n_core <- 4
plan(multicore, workers = n_core)

path_figs <- here("figs", "model_comparison")

## for reading HBM model res:
model_info <- expand.grid(model_nm = models_hbm, response = responses, roi_nm = core32_nms)

## read data for summarystat model, subset relevant rows/cols, get subject-level stats:

d_summarystat <-
  fread(file.path(path_out, "spatial", "projections__stroop__rda__n_resamples100__demean_run__cv_allsess.csv"))
d_summarystat <- d_summarystat %>% rename(ridge = value.ridge, rda = value.rda)
cols_keep <- c("roi", responses, "variable", "trial", "subj", "wave")
d_summarystat <- d_summarystat[test_session == session & roi %in% core32_nms, ..cols_keep] %>% na.omit()
d_summarystat <-
  melt(d_summarystat, id.vars = c("roi", "variable", "trial", "subj", "wave"), variable.name = "response")
s_summarystat_subj <- d_summarystat[, .(value = mean(value)), by = c("variable", "subj", "wave", "roi", "response")]

## read already-summarized data, inc. model comparison stats and diagnostics

misc <-
  readRDS(file.path(path_out, "inferential", atlas_nm, "core32_stats.rds")) %>%
  bind_rows(.id = "model__response__session") %>%
  separate(model__response__session, c("model", "response", "session"), sep = "__") %>%
  filter(session %in% .env$session) %>%
  mutate(model = factor(model, .env$models, ordered = TRUE)) %>%
  mutate(response = factor(response, .env$responses, ordered = TRUE)) %>%
  arrange(model, region, response) %>%
  rename(roi_nm = region, model_nm = model)



# ELPD-LOO model comparison ----

# For rda and uv respectively

p_elpd <- misc %>%
  filter(Term == "elpd_loo") %>%
  ggplot(aes(model_nm, Estimate)) +
  geom_line(aes(group = roi_nm)) +
  geom_errorbar(aes(ymin = Estimate - `Est.Error`, ymax = Estimate + `Est.Error`), width = 0.1) +
  scale_y_continuous(
    trans = arsinh, breaks = c(0, -3000, -5000, -10000, -15000, -20000),
    labels = scales::scientific) +
  scale_x_discrete(labels = names(models_hbm)) +
  theme(legend.position = "none") +
  labs(x = "HBM Model", y = "ELPD LOO") +
  facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate"))) +
  theme(axis.text.x = element_text(size = 7))

# Differences between reduced models and full model, excluding fixed_sigma
diff_dat <- get_diff_dat(misc %>% filter(Term == "elpd_loo"), name_term = "model_nm", id_term = c("roi_nm", "response"))

# Scatterplot of differences
p_elpd_diff <- diff_dat %>%
  filter(!model_nm %in% "fixed_sigma - full") %>%
  pivot_wider(names_from = model_nm, values_from = Estimate) %>%
  ggplot(aes(x = `no_lscov_symm - full`, y = `no_lscov - full`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "ILS Sym - Full", y = "ILS  - Full") +
  theme(legend.position = "none") +
  facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")))

ggsave(file.path(path_figs, "elpd.pdf"), p_elpd + p_elpd_diff, width = 8, height = 2.5)


# Population-level Stroop effects ----

## get summarystat model ests:
s_summarystat_pop <- s_summarystat_subj %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(hilo = hi - lo) %>%
  group_by(subj, roi, response) %>%  ## average over wave
  summarize(hilo = mean(hilo)) %>%
  group_by(roi, response) %>%  ## average over subj
  summarize(
    Estimate = mean(hilo),
    Est.Error = sd(hilo) / sqrt(n()),
    CI.Lower = Estimate - Est.Error * 1.96,
    CI.Upper = Estimate + Est.Error * 1.96
  ) %>%
  rename(roi_nm = roi) %>%
  mutate(model_nm = "summarystat")

## get HBM ests:
model_info_fixef <- model_info
model_info_fixef$hyp[model_info_fixef$model_nm %in% c("no_lscov_symm", "fixed_sigma")] <-
  "(hilo_wave1 + hilo_wave1)/2 = 0"
model_info_fixef$hyp[model_info_fixef$model_nm %in% c("no_lscov", "full")] <-
  "(hi_wave1 + hi_wave2)/2 - (lo_wave1 + lo_wave2)/2 = 0"
s_pop_hbm <- future_pmap_bind(
  model_info_fixef,
  read_summarize_hbm,
  session = session,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(x, hyp) {
    hypothesis(x, hyp)$hypothesis[, c("Estimate", "Est.Error", "CI.Lower", "CI.Upper")]
  }
)

## bind:
s_pop <- full_join(s_pop_hbm, s_summarystat_pop)
s_pop <- s_pop %>% mutate(model_nm = factor(model_nm, levels = models))  ## reorder factor for plotting
s_pop <- as.data.table(s_pop)

## Mean of population-level Stroop effect

# Bar plot
p_pop_bars_mu <- enlist(responses)
for (res in responses) {
  res_lab <- switch(res, uv = "Univariate Response", rda = "Multivariate Response")
  p_pop_bars_mu[[res]] <- s_pop %>%
    filter(response == res) %>%  ## for rda and uv
    ggplot(aes(model_nm, Estimate, fill = model_nm)) +
    geom_hline(yintercept = 0) +
    geom_col(width = 0.5, color = "black") +
    geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), width = 0.1, color = "black") +
    scale_fill_manual(
        values = colors_models,
        guide = guide_legend(title = "models", title.position = "top", direction = "horizontal"),
        labels = names(models)) +
    theme(legend.position = c(0.7, 0.07), legend.direction = "horizontal", axis.text.x = element_blank()) +
    facet_wrap(vars(roi_nm), labeller = labeller(roi_nm = function(x) gsub("17Networks_", "", x))) +
    labs(title = res_lab, x = element_blank(), y = "Mean(Population-Level Stroop Effect) + 95% CI", fill = "models") +
    theme(strip.text = element_text(size = 6), legend.text = element_text(size = 6))
}

# Line plot
p_pop_box_mu <- s_pop %>%
  ggplot(aes(model_nm, Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
  geom_line(aes(group = roi_nm), alpha = 0.5, color = "dodgerblue") +
  scale_color_viridis_d() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = names(models)) +
  labs(x = "Model", y = "Mean(Population-Level Stroop Effect) + 95% CI") +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")), scales = "free_y")

## Standard error of population-level Stroop effect

# Bar plot
p_pop_bars_se <- enlist(responses)
for (res in responses) {
  res_lab <- switch(res, uv = "Univariate Response", rda = "Multivariate Response")
  p_pop_bars_se[[res]] <- s_pop %>%
    filter(response == res) %>%  ## for rda and uv
    ggplot(aes(model_nm, Est.Error, fill = model_nm)) +
    geom_hline(yintercept = 0) +
    geom_col(width = 0.5, color = "black") +
    scale_fill_manual(
        values = colors_models,
        guide = guide_legend(title = "models", title.position = "top", direction = "horizontal"),
        labels = names(models)) +
    theme(legend.position = c(0.7, 0.07), legend.direction = "horizontal", axis.text.x = element_blank()) +
    # theme(legend.position = "none", axis.text.x = element_text(angle = axis_text_x_angle)) +
    facet_wrap(vars(roi_nm), labeller = labeller(roi_nm = function(x) gsub("17Networks_", "", x))) +
    labs(title = res_lab, x = element_blank(), y = "SE(Population-Level Stroop Effect) + 95% CI", fill = "models") +
    theme(strip.text = element_text(size = 6), legend.text = element_text(size = 6))
}

# Line plot
p_pop_box_se <- s_pop %>%
  ggplot(aes(model_nm, Est.Error)) +
  geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
  geom_line(aes(group = roi_nm), alpha = 0.5, color = "dodgerblue") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = names(models)) +
  labs(x = "Model", y = "SE(Population-Level Stroop Effect)") +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")), scales = "free_y")

p_pop_bars <-
    (p_pop_bars_mu$uv + theme(legend.position = "none") + p_pop_bars_mu$rda + theme(legend.position = "none")) /
    (p_pop_bars_se$uv + theme(legend.position = "none") + p_pop_bars_se$rda)
p_pop_box <- p_pop_box_mu / p_pop_box_se
ggsave(file.path(path_figs, "stroop_pop_bars.pdf"), p_pop_bars, width = 12, height = 12)
ggsave(file.path(path_figs, "stroop_pop_box.pdf"), p_pop_box, width = 7, height = 6)



# Test-retest reliability ----

## get summarystat model ests:

s_summarystat_subj_trr <- s_summarystat_subj %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(hilo = hi - lo) %>%
  pivot_wider(id_cols = c("subj", "roi", "response"), names_from = wave, values_from = hilo) %>%
  rename(roi_nm = roi) %>%
  mutate(model_nm = "summarystat")

## get summarystat model ests:

r_summarystat <- s_summarystat_subj_trr %>%
  group_by(roi_nm, model_nm, response) %>%
  summarize(r = cor(wave1, wave2))

## get HBM posteriors:

s_hbm_posterior <- future_pmap_bind(
  model_info,
  read_summarize_hbm,
  base_path = file.path(path_out, "inferential"),
  session = session,
  atlas_nm = atlas_nm,
  sum_fun = function(mdl) {
    mu <- ranef(mdl, summary = FALSE)$subj  ## for conditional TRR
    S <- VarCorr(mdl, summary = FALSE)$subj  ## for marginal TRR
    if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
      ## no ls cov symm and fixed sigma
      mu_stroop1 <- mu[, , "hilo_wave1"]
      mu_stroop2 <- mu[, , "hilo_wave2"]
      trr <- S$cor[, "hilo_wave1", "hilo_wave2"]
    } else {
      ## no ls cov and full
      conditions <- dimnames(S$cov)[[2]]
      n_resamples <- dim(S$cov)[1]
      W <- rbind(
        ("hi_wave1" == conditions) - ("lo_wave1" == conditions),
        ("hi_wave2" == conditions) - ("lo_wave2" == conditions)
      )
      trr <- rep(NA_real_, n_resamples)
      for (ii in seq_len(n_resamples)) {
        Sigma <- tcrossprod(tcrossprod(W, S$cov[ii, , ]), W)
        trr[ii] <- cov2cor(Sigma)[1, 2]
      }
    }
    data.table(trr = trr, resample = seq_along(trr))
  }
)

s_hbm_posterior <- s_hbm_posterior %>%
  mutate(model_nm = factor(model_nm, levels = models_hbm, ordered = TRUE))  ## reorder factor for plotting

## Per-parcel posterior distributions

p_trr_dens <- s_hbm_posterior %>%
  ggplot(aes(trr, color = model_nm, linetype = response)) +
  # geom_rect(xmin = 0.7, xmax = 1, ymin = 0, ymax = Inf, fill = "grey50", color = "grey50", alpha = 0.5) +
  geom_density(linewidth = 0.5) +
  geom_vline(data = r_summarystat, aes(xintercept = r, linetype = response), linewidth = 0.3) +
  facet_wrap(vars(roi_nm), labeller = labeller(roi_nm = function(x) gsub("17Networks_", "", x))) +
  scale_color_manual(values = colors_models[models_hbm], labels = names(models), name = NULL) +
  scale_linetype_discrete(labels = c("Univariate", "Multivariate"), name = NULL) +
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

ggsave(file.path(path_figs, "stroop_trr_dens.pdf"), p_trr_dens, width = 7, height = 5)

## Statistics of posterior distributions

stat_names <- c(map = "MAP", mean = "mean", sd = "sd", fifth = "5%-ile", nintyfifth = "95%-ile")
stat_funs <- c(max_aposteriori, mean, sd, function(x) quantile(x, 0.05), function(x) quantile(x, 0.95))
stat_ranges <- list(c(-1, 1), c(-1, 1), c(0, 0.7), c(-1, 1), c(-1, 1))
p_trr_stats <- enlist(stat_names)
add_sumstat <- FALSE
for (i in seq_along(stat_funs)) {
  tmp_dat <- s_hbm_posterior[, .(trr = stat_funs[[i]](trr)), by = c("model_nm", "roi_nm", "response")]
  if (stat_names[i] %in% c("MAP", "mean") & add_sumstat) {
    tmp_dat <- tmp_dat %>%
      full_join(r_summarystat %>% rename(trr = r)) %>%
      mutate(model_nm = factor(model_nm, levels = models, ordered = TRUE))
  } else {
    tmp_dat <- tmp_dat %>%
      mutate(model_nm = factor(model_nm, levels = models_hbm, ordered = TRUE))
  }
  f1 <- tmp_dat %>%
    ggplot(aes(model_nm, trr, color = roi_nm)) +
    geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
    geom_line(aes(group = roi_nm), alpha = 0.5, color = "dodgerblue") +
    scale_y_continuous(limits = stat_ranges[[i]]) +
    scale_x_discrete(labels = names(models)) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 7)) +
    labs(y = paste0(stat_names[i], "(TRR)")) +
    facet_wrap(~response, labeller = labeller(response = c(uv = "Univariate", rda = "Multivariate")), scales = "free_y")
  f2 <- tmp_dat %>%
    pivot_wider(names_from = response, values_from = trr) %>%
    rename(model = model_nm, roi = roi_nm) %>%
    ggplot(aes(rda, uv)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_viridis_d(guide = "none") +
    scale_x_continuous(limits = stat_ranges[[i]], breaks = c(-1, 0, 1)) +
    scale_y_continuous(limits = stat_ranges[[i]], breaks = c(-1, 0, 1)) +
    labs(x = "Multivariate", y = "Univariate") +
    facet_grid(cols = vars(model), labeller = labeller(model = setNames(names(models), models)))
  p_trr_stats[[i]] <- list(f1, f2)
}

# pwalk(
#   list(p_trr_stats, names(stat_names)),
#   \(plt, nm) ggsave(file.path(path_figs, paste0("stroop_trr_", nm, ".pdf")), plt, width = 8, height = 4)
# )



# Variability ratio log(sd(trial)/sd(subj)) ----

## get HBM posteriors:

s_hbm_subj_ratio <- future_pmap_bind(
  model_info,
  read_summarize_hbm,
  session = session,
  base_path = file.path(path_out, "inferential"),
  atlas_nm = atlas_nm,
  sum_fun = function(x) {
    ranefs <- ranef(x, summary = FALSE)$subj
    if ("sd_subj__hilo_wave1" %in% variables(x)) {
      stroop1 <- ranefs[, , "hilo_wave1"]
      stroop2 <- ranefs[, , "hilo_wave2"]
      if ("b_sigma_Intercept" %in% variables(x))  {
        hyp <- "exp(sigma_Intercept) = 0"
      } else {
        hyp <- "(exp(sigma_hilo_wave1) + exp(sigma_hilo_wave2) + exp(sigma_mean_wave1) + exp(sigma_mean_wave2)) / 4 = 0"
      }
    } else {
      stroop1 <- ranefs[, , "hi_wave1"] - ranefs[, , "lo_wave1"]
      stroop2 <- ranefs[, , "hi_wave2"] - ranefs[, , "lo_wave2"]
      hyp <- "(exp(sigma_hi_wave1) + exp(sigma_hi_wave2) + exp(sigma_lo_wave1) + exp(sigma_lo_wave2)) / 4 = 0"
    }
    sigma_stroop1 <- sqrt(Var(stroop1, 1))
    sigma_stroop2 <- sqrt(Var(stroop2, 1))
    sigma_stroop <- (sigma_stroop1 + sigma_stroop2)/2
    sigma_resid <- hypothesis(x, hyp)$samples$H1
    ratio <- sigma_resid / sigma_stroop

    data.table(
      sigma_resid = sigma_resid, 
      sigma_stroop = sigma_stroop, 
      ratio = sigma_resid / sigma_stroop,
      resample = seq_along(ratio)
    )

  }
)

s_hbm_subj_ratio <- s_hbm_subj_ratio %>%
  mutate(model_nm = factor(model_nm, levels = models_hbm, ordered = TRUE))

## Per-parcel posterior distributions

p_ratio_dens <- s_hbm_subj_ratio %>%
  ggplot(aes(log(ratio), color = model_nm, linetype = response)) +
  geom_density(linewidth = 0.5) +
  facet_wrap(vars(roi_nm), labeller = labeller(roi_nm = function(x) gsub("17Networks_", "", x))) +
  scale_color_manual(values = colors_models[models_hbm], labels = names(models), name = NULL) +
  scale_linetype_discrete(labels = c("Univariate", "Multivariate"), name = NULL) +
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

ggsave(file.path(path_figs, "stroop_ratio_dens.pdf"), p_ratio_dens, width = 7, height = 5)

## Statistics of posterior distributions

stat_names <- c(map = "MAP", mean = "mean", sd = "sd", fifth = "5%-ile", nintyfifth = "95%-ile")
stat_funs <- c(max_aposteriori, mean, sd, function(x) quantile(x, 0.05), function(x) quantile(x, 0.95))
stat_ranges <- list(c(0, 4.5), c(0, 4.5), c(0, 1), c(0, 4.5), c(0, 5))
p_ratio_stats <- enlist(stat_names)
for (i in seq_along(stat_funs)) {
  tmp_dat <- s_hbm_subj_ratio[, .(ratio = stat_funs[[i]](log(ratio))), by = c("model_nm", "roi_nm", "response")]
  tmp_dat <- tmp_dat %>%
    mutate(model_nm = factor(model_nm, levels = models_hbm, ordered = TRUE))
  f1 <- tmp_dat %>%
    ggplot(aes(model_nm, ratio, color = roi_nm)) +
    geom_boxplot(fill = "grey60", color = "black", width = 1/3) +
    geom_line(aes(group = roi_nm), alpha = 0.5, color = "dodgerblue") +
    scale_y_continuous(limits = stat_ranges[[i]]) +
    scale_x_discrete(labels = names(models)) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 7)) +
    labs(y = paste0(stat_names[i], "(Variab. Ratio)")) +
    facet_wrap(vars(response),
      labeller = labeller(
        response = c(
          uv = paste0("Univariate\n", stat_names[i], "(Variab. Ratio)"),
          rda = paste0("Multivariate ", stat_names[i], "(Variab. Ratio)")
        )
      ),
      scales = "free_y")
  f2 <- tmp_dat %>%
    pivot_wider(names_from = response, values_from = ratio) %>%
    rename(model = model_nm, roi = roi_nm) %>%
    ggplot(aes(rda, uv)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_viridis_d(guide = "none") +
    scale_x_continuous(limits = stat_ranges[[i]]) +
    scale_y_continuous(limits = stat_ranges[[i]]) +
    labs(x = "Multivariate", y = "Univariate") +
    facet_grid(cols = vars(model), labeller = labeller(model = setNames(names(models_hbm), models_hbm)))
  p_ratio_stats[[i]] <- list(f1, f2)
}

# pwalk(
#   list(p_ratio_stats, names(stat_names)),
#   \(plt, nm) ggsave(file.path(path_figs, paste0("stroop_ratio_", nm, ".pdf")), plt, width = 8, height = 4)
# )

## arrange ----

p_stats <-
  (p_trr_stats$MAP[[1]] + p_trr_stats$MAP[[2]]) /
  (p_trr_stats$`5%-ile`[[1]] + p_trr_stats$`5%-ile`[[2]]) /
  (p_ratio_stats$MAP[[1]] + p_ratio_stats$MAP[[2]]) /
  (p_ratio_stats$`5%-ile`[[1]] + p_ratio_stats$`5%-ile`[[2]])
ggsave(file.path(path_figs, paste0("stroop_trratio_stats.pdf")), p_stats, width = 8, height = 8)
