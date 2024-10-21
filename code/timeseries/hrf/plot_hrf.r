library(here)
library(rhdf5)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(colorspace)
library(patchwork)

source(here("code", "_constants.R"))
source(here("code", "_subjects.R"))
source(here("code", "_atlas.R"))
source(here("code", "_paths.R"))
source(here("code", "timeseries", "_utils_fmri.R"))
source(here("code", "inferential", "_utils_viz.R"))
source(here("code", "inferential", "_parameters_viz.R"))

theme_set(theme_default())


## functions ----

prep_xmat_for_plotting <- function(x, trs = target_trs$Stroop, n_trials = 3, pinv = FALSE) {
    if (pinv) {
        nms <- names(x)
        x <- as.data.table(cbind(x$tr, t(MASS::ginv(as.matrix(x[, -"tr"])))))
        names(x) <- nms
    }
    x <- melt(x, id.var = c("tr", "Run#1Pol#0"))
    x$is_target_tr <- FALSE
    for (trial in 1:n_trials) {
        x[
            variable == paste0("trial", trial, "#0") & tr %in% (trs + n_trials * (trial - 1)),
            is_target_tr := TRUE
            ]
    }
    x
}


## read ----

lss <- fread(here("out", "timeseries", "glm_comparison_lssep_1rpm-vs-selav.csv")) %>% rename(lss = glmall, r_spatial_lss = r_spatial, r_temporal_lss = r_temporal)
lsa <- fread(here("out", "timeseries", "glm_comparison_lsall_1rpm-vs-selav.csv")) %>% rename(lsa = glmall, r_spatial_lsa = r_spatial, r_temporal_lsa = r_temporal)
glmcomp <- lss[lsa, on = c("parcel", "network", "idx", "hemi", "subject", "wave", "session", "selav")]
glmcomp <- glmcomp[complete.cases(glmcomp)]
glmcomp <- glmcomp %>% mutate(repetition = factor(ifelse(wave == "wave1", "test", "retest"), levels = c("test", "retest")))

fir <- h5read(here::here("out", "timeseries", "deconvolved.h5"), "data")
fir <- as.data.table(lapply(fir, c))

x_trial <- fread(here("code", "timeseries", "hrf", "x_trial.1D.csv")) %>% 
    rename(tr = V1) %>% 
    mutate(tr = tr - 1) %>% 
    filter(tr < 13)
x_trial3 <- fread(here("code", "timeseries", "hrf", "x_3trials.1D.csv")) %>%
    rename(tr = V1) %>%
    mutate(tr = tr - 1) %>%
    filter(tr < 21)




## plot ----


## plot FIR model estimates

hilo <- fir[roi %in% core32, .(b = mean(b)), by = c("tr", "congruency", "session")]
hilo %>%
  ggplot(aes(as.factor(tr), b, linetype = congruency, color = session)) +
  geom_line(aes(group = interaction(congruency, session)), linewidth = 1) +
  scale_color_brewer(type = "qual") +
  labs(title = "deconvolved responses (FIR model)", x = "TR rel. stimulus onset", y = "mean(beta)")
ggsave(here("figs", "hrf", "deconvolved_responses.pdf"), width = 5, height = 3)

## plot FIR and fixed-shape regressors

hilo_w <- dcast(fir, roi + pc + tr + session + subj ~ congruency, value.var = "b")
hilo_w[, stroop_eff := incon - congr]
hilo_agg <- hilo_w[roi %in% core32, .(stroop_eff = mean(stroop_eff)), by = c("tr", "subj")]
scale_fir <- max(hilo_agg[, .(stroop_eff = mean(stroop_eff)), by = c("tr")]$stroop_eff)
hilo_agg[, stroop_eff := stroop_eff / scale_fir, by = "subj"]
hilo_agg <- rbind(hilo_agg, hilo_agg[tr == 1] %>% mutate(tr = 0, stroop_eff = 0))
hilo_agg <- hilo_agg[x_trial[, c("tr", "block#0")], on = "tr"]

p_decon <- hilo_agg %>%
    ggplot(aes(tr, stroop_eff)) +
    geom_vline(xintercept = target_trs$Stroop, linetype = "dotted", color = "grey30") +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.25) +
    stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.4) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(data = . %>% filter(tr %in% target_trs$Stroop), fun = mean, geom = "point", size = 1) +
    geom_line(data = x_trial, aes(y = `block#0`), color = "firebrick", linewidth = 1) +
    geom_point(data = x_trial[tr %in% target_trs$Stroop], aes(y = `block#0`), color = "firebrick", size = 1) +
    scale_x_continuous(breaks = 0:14) +
    labs(x = "TR rel. stimulus onset", y = "BOLD Response\n(max normalized)") +
    annotate(geom = "text", x = 6, y = 1.2, label = "Assumed HRF", hjust = 0, vjust = 1, color = "firebrick", size = 2) +
    annotate(geom = "text", x = 6, y = 1.2, label = "\nEstimated Stroop\nContrast", hjust = 0, vjust = 1, size = 2)
ggsave(here("figs", "hrf", "deconvolved_responses_stroop.pdf"), width = 3, height = 2)


## plot predictions and pseudoinverse for sequential trials

x_trial3_l <- prep_xmat_for_plotting(x_trial3)
p_x <- x_trial3_l %>%
    ggplot(aes(tr, value, color = variable, alpha = variable)) +
    geom_line(aes(group = variable)) +
    geom_point(data = . %>% filter(value > 0.8)) +
    scale_color_manual(values = c("trial1#0" = "grey20", "trial2#0" = "firebrick", "trial3#0" = "grey20")) +
    scale_alpha_manual(values = c("trial1#0" = 0.25, "trial2#0" = 1, "trial3#0" = 0.25)) +
    theme(legend.position = "none") +
    labs(x = NULL, y = expression(bold(X)))

pinv_trial3_l <- prep_xmat_for_plotting(x_trial3, pinv = TRUE)
p_pinv <- pinv_trial3_l %>%
    ggplot(aes(tr, value, color = variable, alpha = variable)) +
    geom_line(aes(group = variable)) +
    geom_point(data = . %>% filter(value > 0.25)) +
    scale_color_manual(values = c("trial1#0" = "grey20", "trial2#0" = "firebrick", "trial3#0" = "grey20")) +
    scale_alpha_manual(values = c("trial1#0" = 0.25, "trial2#0" = 1, "trial3#0" = 0.25)) +
    theme(legend.position = "none") +
    labs(x = "TR", y = expression((bold(X)^T * bold(X))^-1 * bold(X)^T))

p_xmat <- p_x / p_pinv
ggsave(here("figs", "hrf", "xmat_vs_pinv.pdf"), p_xmat, width = 3, height = 2.5)
pinv_trial3_l[, .(r = cor(is_target_tr, value)), by = "variable"]

p_hrf <- p_decon + p_xmat
ggsave(here("figs", "hrf", "hrf.pdf"), p_hrf, width = 5, height = 2.25)

## compare estimates to selav

p_temporal_lss <- glmcomp %>%
    ggplot(aes(r_temporal_lss)) +
    geom_histogram(position = "identity", bins = 10, fill = "grey40", color = "black") +
    labs(y = "Count\n(subj*sess*rep*parc)", x = expression("cross-trial cor"(beta[LSS], beta[selav]))) +
    coord_cartesian(xlim = c(0, 1))
p_temporal_lsa <- glmcomp %>%
    ggplot(aes(r_temporal_lsa)) +
    geom_histogram(position = "identity", bins = 10, fill = "grey40", color = "black") +
    labs(y = NULL, x = expression("cross-trial cor"(beta[LSA], beta[selav]))) +
    coord_cartesian(xlim = c(0, 1))
p_temporal_lss + p_temporal_lsa


p_box_lss_temporal <- glmcomp %>%
    ggplot(aes(network, r_temporal_lss)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-trial cor"(beta[LSS], beta[selav])), x = "Network (Schaefer-400-17)", title = "LS-Separate") +
    coord_flip(ylim = c(0.24, 1))
p_box_lsa_temporal <- glmcomp %>%
    ggplot(aes(network, r_temporal_lsa)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-trial cor"(beta[LSA], beta[selav])), x = "Network (Schaefer-400-17)", title = "GLM-All with ARMA(1,1)") +
    coord_flip(ylim = c(0.24, 1))
p_temporal <- p_box_lss_temporal + p_box_lsa_temporal + plot_layout(axis_title = "collect")
ggsave(here("figs", "hrf", "glm_comparison_temporal.pdf"), p_temporal, width = 7, height = 3)


p_box_lss_spatial <- glmcomp %>%
    ggplot(aes(network, r_spatial_lss)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-vertex cor"(stroop[LSS], stroop[selav])), x = "Network (Schaefer-400-17)", title = "GLM-Separate") +
    coord_flip(ylim = c(-0.6, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5))
p_box_lsa_spatial <- glmcomp %>%
    ggplot(aes(network, r_spatial_lsa)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-vertex cor"(stroop[LSA], stroop[selav])), x = "Network (Schaefer-400-17)", title = "GLM-All with ARMA(1,1)") +
    coord_flip(ylim = c(-0.6, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5))
p_spatial <- p_box_lss_spatial + p_box_lsa_spatial + plot_layout(axis_title = "collect")
ggsave(here("figs", "hrf", "glm_comparison_spatial.pdf"), p_spatial, width = 7, height = 3)


p_subj_lss <- glmcomp[, .(r = cor(lss, selav)), by = c("parcel", "network", "repetition", "session")] %>%
    ggplot(aes(network, r)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-subject cor"(stroop[LSS], stroop[selav])), x = "Network (Schaefer-400-17)", title = "LS-Separate") +
    coord_flip(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5))
p_subj_lsa <- glmcomp[, .(r = cor(lsa, selav)), by = c("parcel", "network", "repetition", "session")] %>%
    ggplot(aes(network, r)) +
    geom_boxplot(fill = "grey50", width = 0.5, linewidth = 1/2, outlier.size = 0.5, outlier.shape = 21, outlier.fill = "white") +
    facet_grid(vars(repetition), vars(session)) +
    labs(y = expression("cross-subject cor"(stroop[LSA], stroop[selav])), x = "Network (Schaefer-400-17)", title = "GLM-All with ARMA(1,1)") +
    coord_flip(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5))
p_stroop_subj <- p_subj_lss + p_subj_lsa + plot_layout(axis_title = "collect")
ggsave(here("figs", "hrf", "glm_comparison_subj.pdf"), p_stroop_subj, width = 7, height = 3)
p_spatial / p_temporal / p_stroop_subj + plot_layout(guides = "collect")


p_lsa <-
    (p_box_lsa_temporal + labs(title = NULL, x = NULL)) /
    (p_box_lsa_spatial + labs(title = NULL, x = NULL)) /
    (p_subj_lsa + labs(title = NULL, x = NULL)) +
    plot_annotation(title = "GLM-All with ARMA(1,1)", theme = theme(plot.title = element_text(hjust = 0.5)))
p_lss <-
    (p_box_lss_temporal + labs(title = NULL, x = NULL)) /
    (p_box_lss_spatial + labs(title = NULL, x = NULL)) /
    (p_subj_lss + labs(title = NULL, x = NULL)) +
    plot_annotation(title = "Least Squares-Separate", theme = theme(plot.title = element_text(hjust = 0.5)))
p_comparison <- wrap_elements(p_lss) + wrap_elements(p_lsa)
ggsave(here("figs", "hrf", "glm_comparison.pdf"), p_comparison, width = 7, height = 7)




## compare xmats to lsa ---

xmat <- fread(here("code", "timeseries", "hrf", "130518_wave1_Stroop_baseline_xmat.csv"))
cols <- names(xmat)[grepl("alltrials", names(xmat))]
xmat_sub <- xmat[, ..cols]
avg <- readRDS(here("out", "timeseries", "130518/RESULTS/Stroop/baseline_null_2rpm_wave1", "errts_averaging_matrix.rds"))[1:540, 1:108]
dim(avg)
dim(xmat)

vifs <- car::vif(lm(rnorm(nrow(xmat)) ~ ., data = xmat))
range(vifs[grepl("alltrials", names(vifs))])
pinv <- t(MASS::ginv(as.matrix(xmat)))
colnames(pinv) <- names(xmat)
pinv_cols <- pinv[, cols]
res <- vector("list", length(cols))
for (i in seq_along(cols)) {
    #fit <- lm(unlist(xmat_sub[, ..i]) ~  avg[, i])
    fit <- lm(pinv_cols[, i] ~  avg[, i])
    res[[i]] <- data.table(b = coef(fit)[2], r2 = summary(fit)$r.squared)
}
res <- rbindlist(res)

