# Test-Retest Reliability

Test-Retest Reliability of fMRI measures in Dual Mechanisms of Cognitive Control dataset

## Repo Organization

- all analysis scripts in `code/`
- scripts in `code` write to `out/`
- `in/` contains behavioral data files, subject lists
- files in `in/` are not modified by any scripts

## Key Factors in DMCC Dataset

- trialtype: "high demand" (difficult trials), "low demand" (easy trials)
- run: run1, run2
- task: Axcpt, Cuedts, Stern, Stroop
- session: baseline, proactive, reactive
- wave: wave1, wave2 (wave1 = test, wave2 = retest)

## Relevant Links, Paths

- [test-retest subject list and QC grades](https://3.basecamp.com/3758557/buckets/3792852/messages/4106700214)
- [first-pass reliability analyses](https://3.basecamp.com/3758557/buckets/3792852/messages/3983554628)
- minimally preprocessed BOLD timeseries (TR by vertex)
  - wave 1: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS
  - wave 2: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase3/fMRIPrep_AFNI_ANALYSIS
- trial onset times
  - /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS

## General Pipeline

- Get the trial-level activation of each voxel:
  - Run `code/timeseries/null_2rpm/run_deconvolve.sh` to detrend the timeseries. The script will read the subject list in `in/` and save the detrended timeseries to `out/timeseries/`.
  - Run `code/timeseries/selav/selav_single_trials.R` to get the trial-level statistics (selective averaging). The script will read the `.gii` files in the GLM output folders in `out/timeseries/`, e.g., `out/timeseries/130518/RESULTS/Stroop/baseline_null_2rpm_wave1/errts_L.gii` (along with `errts_R.gii`), and save the results in `errts_trials_target_epoch.RDS` under the same folder.
- Calculate the univariate or multivariate trial-level statistics for each parcel:
  - ~~Run `code/spatial/univariate_means.R` to get the univariate statistics for each parcel (averaging across all voxels). The script will read the `errts_trials_target_epoch.RDS` files generated in the last step and save the results in `trial-means_schaefer2018_17_400_fsaverage5-parcel_resid-errts.csv` in the same folder.~~
  - Run `in/behav/write_behav_subjs.R` to generate trial-level behavioral data from `in/behav/orig/` to CSV files like `behavior-and-events_wave13_Stroop.csv` under `in/behav`.
  - Run `code/spatial/multivariate_projections_within-task_cv-sessions.R` to generate univariate and multivariate statistics for reliability analysis. The script reads the behavior CSVs and outputs to a csv like `out/spatial/projections__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12.csv`
  - Recode "wave2/3" and "wave1/3" to "wave1/2" and gather all three CSVs into `out/spatial/projections__stroop__rda__n_resamples100__divnorm_vertex__cv_allsess.csv` by `code/spatial/gather_csvs.R`
- Fit hierarchical models for the statistics:
  - Warning: fitting the bayesian model is VERY SLOW and usually some jobs fail at the first run!
  - Run `code/inferential/estimate_reliability.R` to fit the models. The fitted models will be saved to `out/inferential/schaefer2018_17_400_fsaverage5/NAMES_OF_ROI/`.
    - Each model is about 36M and requires about an hour to train (with 4 chains on 4 cores) without further parallelization.
    - With `mclapply()` over the ROIs we can fit 32 models in about 15 hours. 
    - Currently, we only focus on "baseline" session and the core32 ROIs.
- Generate summary statistics and plots:
  - Run `code/inferential/summarize_models.R` to extract the statistics of all models into `out/inferential/schaefer2018_17_400_fsaverage5/core32_stats.rds`, which is a list named by `(model_name, "__", response_name, "__", session)`. Each element in the list is a tibble, containing the statistics.
  - `knitr::knit()` this file `reports/model_comparison.rmd` to generate the report `reports/model_comparison.md` and relevant figures under `reports/figure`

## Method

### Subjects, tasks, data collection, parcellation and definition of ROIs, etc.

(similar to previous DMCC papers)

### fMRI preprocessing and single trial level activation

- Detrending
- Selective averaging
  - Rationale
  - Identifying the frames-of-interest

### Summary and multivariate statistics

- Removing the mean of class means per vertex from each run
- Exploring univariate pre-whitening (maybe we should move this part to supplementary?)
- Bad trial/vertex removal
- Stratefied sampling
- RDA:
  - Rationale
  - Implementation (package and hyperparameters)
- Ridge regression (move to supplementary?)

### Hierarchical Bayesian Modeling

- Models:
  - full
  - no_lscov
  - no_lscov_symm
  - fixed_sigma
- Implementation
- Model comparison and diagnostic statistics

### Statistics

- Population-level effects
- Variance ratio
- TRR:
  - Pearson correlation (ICC(3, 1))
  - Posterior summary statistics: mean / median / mode