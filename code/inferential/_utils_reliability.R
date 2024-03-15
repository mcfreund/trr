library(dplyr)
library(data.table)
library(furrr)
library(brms)

# max_aposteriori <- function(x) {
#   d <- density(x)
#   d$x[which.max(d$y)]
# }


add_names_and_bind <- function(res, .l) {
  for (i in seq_along(res)) {
    res[[i]] <- as.data.table(cbind(res[[i]], .l[i, ]))
  }
  rbindlist(res)
}


future_pmap_bind <- function(.x, .f, ...) add_names_and_bind(future_pmap(.x, .f, ...), .x)


get_network <- function(x, get_8 = FALSE) {
  res <- gsub("17Networks_.H_([[:alpha:]]*)_.*", "\\1", x)
  if(get_8) res <- sub(".$|Par|Cent|Peri", "", res)
  res
}


read_summarize_hbm <- function(model_nm, response, session, base_path, atlas_nm, roi_nm, sum_fun = identity, ...) {
  source(here::here("code", "inferential", "_parameters_reliability.R"))
  
  file_name <- paste0(get_model_info(model_nm, response, session)$model_prefix, ".rds")
  full_path <- file.path(base_path, atlas_nm, roi_nm, file_name)
  mod <- readRDS(full_path)
  sum_fun(mod, ...)

}


extract_popef <- function(mdl) {
  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
    ## no ls cov symm and fixed sigma
    hyps <- c(pop_stroop = "(hilo_wave1 + hilo_wave1) / 2 = 0")
  } else {
    ## no ls cov and full
    hyps <- c(pop_stroop = "(hi_wave1 + hi_wave2) / 2 - (lo_wave1 + lo_wave2) / 2 = 0")
  }

  d <- as.data.table(hypothesis(mdl, hyps)$samples)
  names(d) <- names(hyps)
  d[, resample_idx := seq_len(.N)]
  melt(d, id.vars = "resample_idx", value.name = "value", variable.name = "statistic")

}


extract_trr <- function(mdl) {
  #mu <- ranef(mdl, summary = FALSE)$subj  ## for conditional TRR
  s <- VarCorr(mdl, summary = FALSE)$subj  ## for marginal TRR

  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
    ## no ls cov symm and fixed sigma
    trr <- s$cor[, "hilo_wave1", "hilo_wave2"]
  } else {
    ## no ls cov and full
    conditions <- dimnames(s$cov)[[2]]
    n_resamples <- dim(s$cov)[1]
    w <- rbind(
      ("hi_wave1" == conditions) - ("lo_wave1" == conditions),
      ("hi_wave2" == conditions) - ("lo_wave2" == conditions)
    )
    trr <- rep(NA_real_, n_resamples)
    for (ii in seq_len(n_resamples)) {
      sigma <- tcrossprod(tcrossprod(w, s$cov[ii, , ]), w)
      trr[ii] <- cov2cor(sigma)[1, 2]
    }
  }

  d <- data.table(trr = trr, resample_idx = seq_along(trr))
  melt(d, id.vars = "resample_idx", value.name = "value", variable.name = "statistic")

}


extract_ratio <- function(mdl) {
  ranefs <- ranef(mdl, summary = FALSE)$subj

  if ("sd_subj__hilo_wave1" %in% variables(mdl)) {
    ## no ls cov symm and fixed sigma
    stroop1 <- ranefs[, , "hilo_wave1"]
    stroop2 <- ranefs[, , "hilo_wave2"]
    if ("b_sigma_Intercept" %in% variables(mdl))  {
      ## fixed sigma
      hyp <- "exp(sigma_Intercept) = 0"
    } else {
      ## no ls cov
      hyp <- "(exp(sigma_hilo_wave1) + exp(sigma_hilo_wave2) + exp(sigma_mean_wave1) + exp(sigma_mean_wave2)) / 4 = 0"
    }
  } else {
    ## no ls cov and full
    stroop1 <- ranefs[, , "hi_wave1"] - ranefs[, , "lo_wave1"]
    stroop2 <- ranefs[, , "hi_wave2"] - ranefs[, , "lo_wave2"]
    hyp <- "(exp(sigma_hi_wave1) + exp(sigma_hi_wave2) + exp(sigma_lo_wave1) + exp(sigma_lo_wave2)) / 4 = 0"
  }

  sigma_stroop1 <- sqrt(mfutils::Var(stroop1, 1))
  sigma_stroop2 <- sqrt(mfutils::Var(stroop2, 1))
  sigma_stroop <- (sigma_stroop1 + sigma_stroop2) / 2
  sigma_resid <- hypothesis(mdl, hyp)$samples$H1
  ratio <- log(sigma_resid / sigma_stroop)

  d <- data.table(
    sigma_resid = sigma_resid,
    sigma_stroop = sigma_stroop,
    ratio = ratio,
    resample_idx = seq_along(ratio)
  )

  melt(d, id.vars = "resample_idx", value.name = "value", variable.name = "statistic")

}


pull_criteria <- function(mdl, criteria_nms = c("loo", "waic")) {
  l <- lapply(criteria_nms, \(x) as.data.table(mdl$criteria[[x]]$estimates, keep.rownames = "Term"))
  rbindlist(l)
}


summarize_posterior_ <- function(x, sum_funs) {
  map_dbl(sum_funs, \(f, samples = x) f(samples))
}


vec2draws <- function(x, n_chains = 4) {

  # Transform vector/matrix/array to draws data frame

  if (length(dim(x)) < 3) {
    x <- as.data.frame(x)
    n_iter <- dim(x)[1] / n_chains
    if (dim(x)[1] %% n_chains) stop("Number of samples is not divisible by number of chains!")
    x[[".chain"]] <- rep(1:n_chains, each = n_iter)
  }
  x <- posterior::as_draws_rvars(x)

  x

}


reformat_for_downstream_analyses <- function(x) {
  x %>%
    rename(region = "roi_nm") %>% ## ggseg plotting functions expect this
    mutate(network = get_network(region))
}


summarize_posteriors <- function(
  data,
  sum_funs = list(
    mean = mean,
    median = median,
    map = \(x) as.numeric(bayestestR::map_estimate(x)),
    hdi95_lower = \(x) as.numeric(bayestestR::hdi(x, 0.95))[2],
    hdi95_upper = \(x) as.numeric(bayestestR::hdi(x, 0.95))[3],
    hdi89_lower = \(x) as.numeric(bayestestR::hdi(x, 0.89))[2],
    hdi89_upper = \(x) as.numeric(bayestestR::hdi(x, 0.89))[3],
    sd = sd,
    iqr = IQR,
    mad = mad,
    q95 = \(x) quantile(x, 0.95),
    q10 = \(x) quantile(x, 0.10),
    q05 = \(x) quantile(x, 0.05),
    rhat = \(x, n_chains_ = n_chains) posterior::rhat(vec2draws(x, n_chains_)[[1]]),
    ess_bulk = \(x, n_chains_ = n_chains) posterior::ess_bulk(vec2draws(x, n_chains_)[[1]]),
    ess_tail = \(x, n_chains_ = n_chains) posterior::ess_tail(vec2draws(x, n_chains_)[[1]])
  ),
  by = c("statistic", "model_nm", "response", "region", "session")
) {
  source(here::here("code", "inferential", "_parameters_reliability.R"))
  data[,
    list(value = summarize_posterior_(value, sum_funs), sum_fun = names(sum_funs)),
    by = c("statistic", "model_nm", "response", "region", "session")]

}


load_summaries <- function(
  data_type, prefixes,
  path = path_reliab,
  read_fun = NULL,
  sum_fun = I,
  ...
) {
  source(here::here("code", "_paths.R"))
  
  # Determine the file extension based on data_type
  ext <- switch(
    data_type,
    "posterior_samples" = ".RDS",
    "posterior_summaries" = ".RDS",
    "summarystat" = ".csv",
    stop("Unsupported data type")
  )
  
  file_names <- paste0(data_type, "_", prefixes, ext)
  
  # Set default read_fun based on file extension
  if (is.null(read_fun)) {
    read_fun <- switch(
      ext,
      ".RDS" = readRDS,
      ".csv" = data.table::fread,
      stop("Unsupported file extension")
    )
  }
  
  map(
    setNames(file.path(path, file_names), prefixes),
    ...,
    \(fn, ...) {
      x <- read_fun(fn)
      sum_fun(x, ...)
    }
  )

}


subset_and_order <- function(x, regions, session, models_order = models_hbm) {
  source(here::here("code", "inferential", "_parameters_reliability.R"))
  # subset:
  x_sub <- x[region %in% regions & session == .S, env = list(.S = I(session))]
  ## order factor cols for plotting:
  if (!is.null(models_order)) {
    x_sub[, model_nm := factor(model_nm, levels = models_order)]
  }

  x_sub

}
