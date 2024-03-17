## reliability model parameters

n_chains <- 4

models_hbm <- c("full", "no_lscov", "no_lscov_symm", "fixed_sigma")
models_hbm <- setNames(models_hbm, c("Full", "ILS", "ILS Sym", "Homog"))
models_hbm_vals2nms <- setNames(names(models_hbm), models_hbm)
models <- c(models_hbm, "Sum Stat" = "summarystat")
models_comparison <- c(ICC = "summarystat", HBM = "no_lscov_symm")

responses <- c("Univariate" = "uv", "Multivariate" = "rda")


get_model_info <- function(model_name, response_name, session_name) {

  if (model_name == "full") {

    res <- list(
      model_prefix =
        paste0(
          "hbm__response_", response_name, "__family_student__sigma_full__ranef_full__lscov_full__session_",
          session_name
        ),
      formula_string =
        paste0(
          " ~ ",
          "0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 + ",
          "(0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 | a | subj)"
        ),
      formula_sigma =
        formula(
          sigma ~ 0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 +
            (0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 | a | subj)
        )
    )

  } else if (model_name == "no_lscov") {

    res <- list(
      model_prefix =
        paste0(
          "hbm__response_", response_name, "__family_student__sigma_full__ranef_full__lscov_none__session_",
          session_name
        ),
      formula_string =
        paste0(
          " ~ ",
          "0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 + ",
          "(0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 | subj)"
        ),
      formula_sigma =
        formula(
          sigma ~ 0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 +
            (0 + lo_wave1 + lo_wave2 + hi_wave1 + hi_wave2 | subj)
        )
    )

  } else if (model_name == "no_lscov_symm") {

    res <- list(
      model_prefix =
        paste0(
          "hbm__response_", response_name, "__family_student__sigma_symm__ranef_symm__lscov_none__session_",
          session_name
        ),
      formula_string =
        paste0(
          " ~ ",
          "0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 + ",
          "(0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)"
        ),
      formula_sigma =
        formula(
          sigma ~ 0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 +
            (0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)
        )
    )

  } else if (model_name == "fixed_sigma") {

    res <- list(
      model_prefix =
        paste0(
          "hbm__response_", response_name, "__family_student__sigma_fixed__ranef_symm__lscov_none__session_",
          session_name
        ),
      formula_string =
        paste0(
          " ~ ",
          "0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 + ",
          "(0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)"
        ),
      formula_sigma = formula(sigma ~ 1)
    )

  }

  res

}

