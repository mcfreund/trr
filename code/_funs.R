read_gifti2matrix <- function(name){
    d <- gifti::read_gifti(name)$data
    matrix(unlist(d, use.names = FALSE), nrow = length(d), byrow = TRUE)
}

read_results <- function(waves, tasks, sessions, subjs, glmname, filename_fun, read_fun, n_cores = 20) {
  

  cl <- parallel::makeCluster(n_cores, type = "FORK")
  doParallel::registerDoParallel(cl)
  res <- 
    foreach::foreach(task_i = seq_along(tasks), .inorder = FALSE, .combine = "c") %:%
    foreach::foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %:%
    foreach::foreach(session_i = seq_along(sessions), .inorder = FALSE, .combine = "c") %:%
    foreach::foreach(subj_i = seq_along(subjs), .inorder = FALSE, .combine = "c") %dopar% {
      
      name_wave_i <- waves[wave_i]
      name_task_i <- tasks[task_i]
      name_session_i <- sessions[session_i]
      name_subj_i <- subjs[subj_i]
      
      filename <- filename_fun(
        name_wave_i = name_wave_i, name_task_i = name_task_i, name_session_i = name_session_i, 
        name_subj_i = name_subj_i, glmname = glmname
      )
      
      filename_full <- 
        here::here(
          "out", "timeseries", name_subj_i, "RESULTS", name_task_i, 
          paste0(name_session_i, "_", glmname, "_", name_wave_i), 
          filename
          )
      
      nm <- paste0(name_wave_i, "_", name_task_i, "_", name_session_i, "_", name_subj_i)
      setNames(list(read_fun(filename_full)), nm)
            
  }
  parallel::stopCluster(cl)
  
  res
  
}


get_subbrick_labels <- function(fname) {
  
  labs <- afni("3dinfo", paste0("-label ", fname))
  labs <- unlist(strsplit(labs, "\\|"))
  
}




get_filenames_results <- function(waves, tasks, sessions, subjs, glmname, filename_fun) {
  
  nms <- mikeutils::combo_paste(waves, tasks, sessions, subjs)
  v <- vector("character", length(nms))
  names(v) <- nms

  for (wave_i in seq_along(waves)) {
    
    name_wave_i <- waves[wave_i]
    
    for (task_i in seq_along(tasks)) {
      
      name_task_i <- tasks[task_i]
      
      for (session_i in seq_along(sessions)) {
        
        name_session_i <- sessions[session_i]
        
        for (subj_i in seq_along(subjs)) {
          
          name_subj_i <- subjs[subj_i]
          
          filename <- filename_fun(
            name_wave_i = name_wave_i, name_task_i = name_task_i, name_session_i = name_session_i, 
            name_subj_i = name_subj_i, glmname = glmname
          )
          
          filename_full <- 
            here::here(
              "out", "glms", name_subj_i, "RESULTS", name_task_i, 
              paste0(name_session_i, "_", glmname, "_", name_wave_i), 
              filename
            )
          
          nm <- paste0(name_wave_i, "_", name_task_i, "_", name_session_i, "_", name_subj_i)
          v[[nm]] <- filename_full
          
        }
        
      }
      
    }
    
  }
  
  v
  
}





read_betas_dmcc <- function(
  .subjs,
  .task,
  .glm,
  .dir
) {
  # .subjs = "130518"
  # .task = "Axcpt"
  # .glm = "baseline_Cues_EVENTS_censored"
  # .dir = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS"
  
  ## initialize array
  
  pick.a.file <- 
    file.path(.dir, .subjs[1], "1TRpK_SURFACE_RESULTS",  .task, paste0(.glm), paste0("STATS_", subjs[1], "_REML_L.func.gii"))
  labs <- afni("3dinfo", paste0("-label ", pick.a.file))
  labs <- unlist(strsplit(labs, "\\|"))
  is.reg <- !grepl("Full|block|Tstat|Fstat", labs)
  tab <- do.call(rbind, strsplit(gsub("_Coef", "", labs[is.reg]), "#"))
  trs <- as.numeric(unique(tab[, 2])) + 1
  regs <- unique(tab[, 1])
  
  n.vertex <- 10242
  n.tr <- length(trs)
  n.reg <- length(regs)
  n.subj <- length(subjs)
  
  betas <- array(
    NA,
    dim = c(n.vertex*2, n.reg, n.tr, n.subj),
    dimnames = list(vertex = NULL, reg = regs, tr = NULL, subj = subjs)
  )
  
  vertex.inds <- cbind(L = 1:n.vertex, R = (n.vertex + 1):(n.vertex * 2))
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 41; hemi.i = "L"
    
    for (hemi.i in c("L", "R")) {
      # hemi.i = "R"
      
      inds <- vertex.inds[, hemi.i]
      
      fname <- file.path(
        .dir, .subjs[subj.i], "1TRpK_SURFACE_RESULTS",  .task, paste0(.glm),  
        paste0("STATS_", subjs[subj.i], "_REML_", hemi.i, ".func.gii")
      )
      
      if (!file.exists(fname)) next
      
      B <- mikeutils::read_gifti2matrix(fname)[is.reg, ]
      
      is.ok.i <- isTRUE(all.equal(dim(B), c(n.reg * n.tr, n.vertex)))
      if (!is.ok.i) stop("mismatched beta array")
      
      
      for (reg.i in seq_len(n.reg)) {
        # reg.i = 1
        
        is.reg.i <- grepl(paste0("^", regs[reg.i], "#"), labs[is.reg])
        B.reg.i <- t(B[is.reg.i, ])
        
        is.ok.ii <- isTRUE(all.equal(dim(betas[inds, reg.i, , subj.i]), dim(B.reg.i)))
        if (!is.ok.ii) stop("mismatched regressor array")
        
        betas[inds, reg.i, , subj.i] <- B.reg.i
        
      }
      
    }
    
  }
  
  betas
  
}

afni <- function (f, args, afni_path = "/usr/local/pkg/linux_openmp_64/", ...) {
    system2(command = paste0(afni_path, f), args = args, stdout = TRUE, ...)
}



## model_name
## object: filename_prefix, formula_response, formula_sigma
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