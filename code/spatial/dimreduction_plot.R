library(here)
library(dplyr)
library(data.table)
library(abind)
library(doParallel)
library(foreach)
library(mfutils)
library(ggplot2)
library(colorspace)
library(patchwork)
source(here("code", "_constants.R"))
source(here("code", "_funs.R"))

sessions <- "baseline"
variable <- "hilo_all"
classes <- c("lo", "hi")  ## -, +
tasks <- "Stroop"
atlas_nm <- "schaefer2018_17_400_fsaverage5"
roi_col <- "parcel"  ## "parcel" or "network"
glm_nm <- "null_2rpm"
resid_type <- "errts"
n_cores <- 18
demean_run <- TRUE
do_waves <- c(1, 2)
subjs <- switch(toString(do_waves),
  "1, 2" = subjs_wave12_complete, "1, 3" = subjs_wave13_all, "2, 3" = subjs_wave23_all
)

## plot settings

theme_set(theme_bw(base_size = 14))
theme_update(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    legend.position = "none"
)

## atlas info and other constants
atlas <- get(atlas_nm)
rois <- unique(atlas$key[[roi_col]])
waves <- waves[do_waves]
n_classes <- length(classes)

## result file name
file_name <- paste0(
  "signalnoise__stroop",
  switch(demean_run + 1, "", "__demean_run"),
  "__baseline__wave", do_waves[1], do_waves[2], ".csv"
)

d <- fread(here("out", "spatial", file_name))[session %in% sessions]
d[, dimension := as.factor(dimension)]
d[, eigvals_scaled := eigvals/sum(eigvals), by = c("roi", "subj", "wave", "session")]
d[, proj_signal_noise_scaled := proj_signal_noise^2/ssq_signal]
#d[roi %in% rois[127] & subj == "130518" & wave == "wave1", .(tot = sum(proj_signal_noise_scaled))]
dl <- melt(d, id.vars = c("roi", "dimension", "subj", "wave", "session"))
dl <- dl[, .(value = mean(value)), by = c("roi", "dimension", "variable", "subj", "session")]  ## average across waves

# ## weights
# weights <- fread("/data/nil-external/ccp/freund/trr/out/spatial/weights__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12.csv")

networks <- rois %>% gsub("17Networks_.H_*", "\\1", .) %>% gsub("([A-z]*)_.*", "\\1", .) %>% gsub("([A-z]*)_.*", "\\1", .)
#key_rois <- rois[core32][1:4]
key_rois <- rois[networks %in% c("ContA", "ContB", "DorsAttnA", "DorsAttnB")]

## plot ----

## spectrum

p_var <-
  dl[
    as.numeric(dimension) < 21 & variable %in% c("eigvals_scaled", "proj_signal_noise_scaled"),
    .(value = mean(value)),
    by = c("roi", "dimension", "variable")
    ] %>%
  ggplot(aes(as.numeric(dimension), value, color = variable)) +
  geom_line(aes(group = paste0(roi, variable)), alpha = 0.05) +
  annotate(geom = "text", x = 2, y = 0.7, label = "noise", color = "#ea00ff", size = 8) +
  annotate(geom = "text", x = 2, y = 0.7, label = "\n\nsignal", color = "#0073ff", size = 8) +
  stat_summary(fun = mean, geom = "line", aes(group = paste0(variable)), size = 1) +
  labs(y = "var / total var", x = "PC") +
  scale_color_manual(values = c(eigvals_scaled = "#ea00ff", proj_signal_noise_scaled = "#0073ff")) +
  scale_x_continuous(trans = "log", breaks = c(1:10, 20))

p_snr <-
  dl[as.numeric(dimension) < 21 & variable %in% c("eigvals", "proj_signal_noise")] %>%
  dcast(... ~ variable, value.var = "value") %>%
  mutate(snr = proj_signal_noise^2 / eigvals) %>%
  .[, .(snr = mean(snr)), by = c("roi", "dimension")] %>%
  ggplot(aes(as.numeric(dimension), snr)) +
  geom_line(aes(group = roi), alpha = 0.1) +
  stat_summary(fun = mean, geom = "line", size = 1, color = "#31a354") +
  labs(y = "SNR", x = "PC") +
  scale_color_manual(values = c(eigvals_scaled = "#ea00ff", proj_signal_noise_scaled = "#0073ff")) +
  scale_x_continuous(trans = "log", breaks = c(1:10, 20))


noise_d1 <- dl[
  dimension == 1 & roi %in% key_rois & variable == "cossim_noise_unif",
  .(value = mean(value)), by = "roi"]

p_noise_unif <-
  dl[as.numeric(dimension) < 21 & roi %in% key_rois & variable %in% "cossim_noise_unif"] %>%
  ggplot(aes(as.numeric(dimension), value)) +
  stat_summary(fun = mean, geom = "line", aes(group = roi), alpha = 0.5, color = "#ea00ff") +
  labs(y = "cos(angle w. unity)", x = "PC") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(trans = "log", breaks = c(1:10, 20))

rois_sorted <-
  dl[roi %in% key_rois & variable %in% c("cossim_signal_unif"), .(value = mean(abs(value))), by = c("roi")] %>%
  arrange(-abs(value)) %>%
  pull(roi)

# p_signal_unif <-
#   dl[
#     dimension == 1 & roi %in% key_rois & variable %in% c("cossim_signal_unif"),
#     .(value = mean(abs(value))),
#     by = c("roi")
#     ] %>%
#   ggplot(aes(value)) +
#   geom_histogram(bins = 15, fill = "#0073ff", color = "black") +
#   labs(y = "num. ROI") +
#   theme(
#     axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y.left = element_blank(),
#     axis.title.y = element_blank()
#   ) +
#   coord_flip(xlim = c(0, 1)) +
#   scale_y_reverse(breaks = c(0, 7, 15))
  
p_hist <-  dl[
    dimension == 1 & roi %in% key_rois & variable %in% c("cossim_signal_unif", "cossim_weights_unif"),
    .(value = mean(abs(value))),
    by = c("roi", "variable")
    ] %>%
  ggplot(aes(value, fill = variable)) +
  geom_histogram(bins = 30, color = "black") +
  labs(y = "num. ROI") +
  scale_fill_manual(values = c(cossim_signal_unif = "#0073ff", cossim_weights_unif = "#31a354")) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y.left = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_flip(xlim = c(0, 1)) +
  scale_y_reverse(breaks = c(0, 20, 60))

p_spectrum <- p_var / p_snr
p_unif <- p_hist + p_noise_unif + plot_layout(ncol = 2, widths = c(1, 5))

ggsave(here("figs", "dimreduction", "spectrum.pdf"), p_spectrum, device = "pdf", width = 5, height = 6)
ggsave(here("figs", "dimreduction", "uniformity.pdf"), p_unif, device = "pdf", width = 5, height = 3)



# table(dl[roi %in% key_rois]$roi, dl[roi %in% key_rois]$dimension)

# dl[
#   roi %in% "17Networks_RH_ContA_PFCl_2" &
#   variable %in% c("eigvals_scaled", "proj_signal_noise_scaled", "proj_weights_noise_scaled")
#   #.(value = mean(value)),
#   #by = c("roi", "dimension", "variable")
#   ] %>%
#   ggplot(aes(as.numeric(dimension), value, color = variable, fill = variable)) +
#   #geom_line(aes(group = paste0(roi, variable)), alpha = 0.05) +
#   #annotate(geom = "text", x = 2, y = 0.7, label = "noise", color = "#ea00ff", size = 8) +
#   #annotate(geom = "text", x = 2, y = 0.7, label = "\n\nsignal", color = "#0073ff", size = 8) +
#   stat_summary(fun.data = mean_cl_boot, geom = "ribbon", aes(group = paste0(variable)), alpha = 0.1, color = "transparent") +
#   stat_summary(fun = mean, geom = "line", aes(group = paste0(variable)), size = 1) +
#   labs(y = "var / total var", x = "PC") +
#   scale_color_manual(values = c(eigvals_scaled = "#ea00ff", proj_signal_noise_scaled = "#0073ff")) +
#   scale_fill_manual(values = c(eigvals_scaled = "#ea00ff", proj_signal_noise_scaled = "#0073ff")) +
#   scale_x_continuous(trans = "log", breaks = c(1:10, 20))



# dl[
#   roi == "17Networks_LH_ContA_PFCl_1" & variable %in% c("eigvals_scaled", "proj_signal_noise_scaled", "proj_weights_noise_scaled"),
#   ] %>%
#   ggplot(aes(as.numeric(dimension), value, color = variable)) +
#   geom_line(aes(group = paste0(subj, variable)), alpha = 0.1) +
#   #annotate(geom = "text", x = 2, y = 0.7, label = "noise", color = "#ea00ff", size = 8) +
#   #annotate(geom = "text", x = 2, y = 0.7, label = "\n\nsignal", color = "#0073ff", size = 8) +
#   stat_summary(fun = mean, geom = "line", aes(group = paste0(variable)), size = 1) +
#   labs(y = "var / total var", x = "PC") +
#   #scale_color_manual(values = c(eigvals_scaled = "#ea00ff", proj_signal_noise_scaled = "#0073ff")) +
#   scale_x_continuous(trans = "log", breaks = c(1:10, 20))
