library(data.table)
library(dplyr)
library(ggplot2)
library(httpgd)
theme_set(theme_bw(base_size = 14))
hgd()

combo_paste <- function(a, b, sep = "_", ...) apply(expand.grid(a, b, ...), 1, paste0, collapse = sep)

## constants

dir_xmats <- "code/timeseries/hrf"
core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 
  340, 345, 349, 350, 351, 352, 354, 361, 365, 387
)
sessions <- c("baseline", "proactive", "reactive")
subjs <- c(
  "107321", "115825", "123117", "130114", "130518", "132017", "135730", "138837", "141422", "150423", "155938", "158136", "160830", "161832",
  "165032", "173738", "178243", "178647", "178950", "182436", "197449", "203418", "204319", "300618", "317332", "346945", "393550", "448347", 
  "580650", "594156", "601127", "672756", "765864", "814649", "849971", "877168", "DMCC1328342", "DMCC1596165", "DMCC1971064", "DMCC2442951", 
  "DMCC2609759", "DMCC2803654", "DMCC2834766", "DMCC3062542", "DMCC3963378", "DMCC4191255", "562345", "DMCC5195268", "DMCC5775387", "DMCC6418065",
  "DMCC6484785", "DMCC6627478", "DMCC6671683", "DMCC6705371", "DMCC6721369", "DMCC7921988", "DMCC8050964", "DMCC8078683", "DMCC8214059", "DMCC8260571", "DMCC9441378", "DMCC9478705"
)
to_drop <- c("233326", "DMCC6371570", "873968", "DMCC1971064", "DMCC9478705")
subjs <- setdiff(subjs, to_drop)
dir_results <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/"
ttypes <- c("biasInCon", "biasCon", "PC50InCon", "PC50Con")
# https://github.com/ccplabwustl/R01/blob/master/Jo/datasetQC/blocksTasksToDrop.txt


## read xmats ----

x_trial <- fread(file.path(dir_xmats, "x_trial.1D.csv"))
image(as.matrix(x_trial[, -"V1"]))


## read data ----

l <- setNames(vector("list", length(combo_paste(sessions, subjs))), combo_paste(sessions, subjs))
for (subj_i in seq_along(subjs)) {
  for (session_i in seq_along(sessions)) {
    
    fname <-  
      paste0(
        dir_results, subjs[subj_i], "/1TRpK_SURFACE_RESULTS/Stroop/", 
        sessions[session_i], "_Congruency_EVENTS_censored/",
        subjs[subj_i], "_timecourses_", sessions[session_i], "_", ttypes,
        "_Coef_tents_Schaefer2018_400Parcels_7Networks_order_10K"
    )
    
    fname_lh <- paste0(fname, "_L.txt")
    fname_rh <- paste0(fname, "_R.txt")
    
    d <- bind_rows(
      melt(
        rbindlist(lapply(fname_lh, fread)), 
        id.vars = c("File", "Sub-brick"), variable.name = "roi", value.name = "b"
        ),
      melt(
        rbindlist(lapply(fname_rh, fread)), 
        id.vars = c("File", "Sub-brick"), variable.name = "roi", value.name = "b"
      )
    )
    
    nm <- paste0(sessions[session_i], "_", subjs[subj_i])
    l[[nm]] <- d
    
    
  }
}
d <- rbindlist(l, idcol = "session_subj")
d <- tidyr::separate(d, session_subj, c("session", "subj"))
d <- as.data.table(d)

## wrangle

d$pc <- ifelse(grepl("bias", d$`Sub-brick`), "bias", "unbias")
d$congruency <- ifelse(grepl("InCon", d$`Sub-brick`), "incon", "congr")
d$tr <- as.numeric(gsub("([0-9].*)(\\[.*\\])", "\\1", d$`Sub-brick`)) + 1
core32_names <- as.character(unique(d$roi)[core32])


## plot ----


hilo_core32 <-
  d[roi %in% core32_names & congruency == "incon", .(b = mean(b)), by = "tr"]$b - 
  d[roi %in% core32_names & congruency == "congr", .(b = mean(b)), by = "tr"]$b
avg_core32 <-
  d[roi %in% core32_names, .(b = mean(b)), by = "tr"]$b + 
  d[roi %in% core32_names, .(b = mean(b)), by = "tr"]$b
hilo <- d[roi %in% core32_names, .(b = mean(b)), by = c("tr", "congruency", "session")]


hilo %>%
  ggplot(aes(as.factor(tr), b, linetype = congruency, color = session)) +
  geom_line(aes(group = interaction(congruency, session)), linewidth = 2) +
  scale_color_brewer(type = "qual")



mat <- cbind(
  hilo_core32,
  hilo_all,
  avg_core32, 
  avg_all, 
  X[, -"tr"]
  )
mat_s <- sweep(mat, 2, apply(mat, 2, max), "/")
mat_s <- rbind(0, mat_s)
l <- mat_s %>%
  as.data.table %>%
  cbind(tr = c(0, X$tr)) %>%
  melt(id.var = "tr")
l %>%
  ggplot(aes(tr, value)) +
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey50") +  
  geom_path(
    data = . %>% filter(variable %in% c("block_true_minus")) %>% mutate(tr = tr - 3.2/1.2),
    size = 1, 
    linetype = "dashed"
    ) +
  geom_path(
    data = . %>% filter(variable %in% c("block_true_minus")) %>% mutate(tr = tr + 3.2/1.2),
    size = 1, 
    linetype = "dashed"
    ) +
  geom_path(
    data = . %>% filter(variable %in% c("hilo_core32")), 
    size = 1, 
    color = "firebrick"
    ) +
  geom_point(
    data = . %>% filter(variable %in% c("hilo_core32"), tr < 6, tr > 2),
    size = 3,
    color = "white",
    fill = "firebrick",
    shape = 21
    ) +
  geom_point(
    data = . %>% filter(variable %in% c("block_true_minus"), tr < 6, tr > 2),
    size = 3,
    color = "white",
    fill = "black",
    shape = 21
    ) +
  geom_path(
    data = . %>% filter(variable %in% c("block_true_minus")), 
    size = 1,
    aes(linetype = variable)
    ) +
  scale_color_brewer(type = "qual", palette = 4) +
  scale_x_continuous(breaks = 0:14) +
  coord_cartesian(xlim = c(0, 14))
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "TR post event onset", y = "BOLD Response (max normalized)") +
  annotate(geom = "text", x = 8, y = 0.5, label = "Predicted HRF (fixed-shape)", hjust = 0, vjust = 1) +
  annotate(geom = "text", x = 8, y = 0.5, label = "\nEstimated Incongruent-Congruent\ncontrast (deconvolved)", color = "firebrick", hjust = 0, vjust = 1) +
  annotate(geom = "text", x = 4, y = 0.5, label = "Window of\nselective\naveraging", color = "grey30", vjust = 1)


## confirm pseudoinverse on adjacent trials weights the target TRs maximally
## - average n TR between trials?
## - 2, 3, 4 TRs



## old ----
dat <- cbind(
  hilo_core32,
  hilo_all,
  avg_core32, 
  avg_all, 
  X[, -"tr"]
  ) %>% 
  scale(center = FALSE) %>%
  as.data.table %>%
  cbind(tr = X$tr) %>%
  melt(id.var = "tr") %>%
  filter(!grepl("TENT|block2", variable))

dat %>%
  ggplot(aes(tr, value)) +
  geom_path(
    data = . %>% filter(variable %in% c("avg_core32", "avg_all", "hilo_core32")), 
    size = 2, 
    aes(color = variable)
    ) +
  geom_path(
    data = . %>% filter(variable %in% c("block_tr", "block_true", "block_true_minus")), 
    size = 2, 
    aes(linetype = variable)
    ) +
  scale_color_brewer(type = "qual", palette = 4) +
  scale_x_continuous(breaks = 1:14) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "tr post event onset", y = "value", color = "data", linetype = "model")


data <- scale(cbind(hilo_core32, avg_core32, avg_all), center = TRUE)
X_s <- as.matrix(X[, .SD, .SDcols = grep("block|spm", colnames(X))])

weights <- reshape2::melt(crossprod(data, X_s), varnames = c("data", "model"))

weights %>%
  ggplot(aes(model, value, color = data)) +
  geom_line(aes(group = data), size = 3) +
  geom_point(size = 5) +
  scale_color_brewer(type = "qual", palette = 4) +
  theme(legend.position = "top")
