# Test-Retest Reliability

[![DOI](https://zenodo.org/badge/460533866.svg)](https://zenodo.org/badge/latestdoi/460533866)

Test-Retest Reliability of fMRI measures in Dual Mechanisms of Cognitive Control dataset.
This code is capable of reproducing all results and figures within [Freund, Chen, Chen, and Braver (2024; DOI: 10.1101/2024.04.24.591032)](https://doi.org/10.1101/2024.04.24.591032).
For questions or concerns please raise an issue on GitHub.


## Setup and downloading required dataset

To reproduce the results and figures in the paper, you will need to:
1. Clone this repository.
2. Download the data derivatives from [https://doi.org/10.5281/zenodo.14043318](https://doi.org/10.5281/zenodo.14043318).
3. Unpack the zipped file of data derivatives into the `out` directory of this repository.
4. Skip down to the `Running Analyses` section of this README -- in particular, work through the `General Pipeline` subsection.

NB: **fMRI data derivatives** includes the output of our:
- timeseries models (vertex by trial matrices of activation estimates per subject and session)
- spatial models (a long-form table of trial-level univariate and multivariate contrasts for each region, session, and subject)
- reliability models (RDS files containing hierarchical model estimates)

The raw or minimally preprocessed fMRI timeseries data (i.e., inputs to timeseries models) will soon be made available on OpenNeuro, organized as a subcomponent of a larger release of the [Dual Mechanisms of Cognitive Control dataset](https://sites.wustl.edu/dualmechanisms/). **If you would like these data before the OpenNeuro release, please contact us and we will fulfill your request directly.**


## Running Analyses

### Repo Organization

- all analysis scripts in `code/`
- scripts in `code` read and write to `out/`
- `in/` contains behavioral data files, subject lists
- files in `in/` are not modified by any scripts

### Key Factors in DMCC Dataset

- trialtype: "high demand" (difficult trials), "low demand" (easy trials)
- run: run1, run2
- task: Axcpt, Cuedts, Stern, Stroop
- session: baseline, proactive, reactive
- wave: wave1, wave2 (wave1 = test, wave2 = retest)

### General Pipeline

#### 1. Get the trial-level activation of each voxel

Run `code/timeseries/null_2rpm/run_deconvolve.sh` to detrend the timeseries.
  - Reads: the subject list in `in/`
  - Writes: the detrended timeseries to `out/timeseries/`

Run `code/timeseries/selav/selav_single_trials.R` to get the trial-level statistics (selective averaging).
  - reads: the `.gii` files in the GLM output folders in `out/timeseries/`
    - e.g., `out/timeseries/130518/RESULTS/Stroop/baseline_null_2rpm_wave1/errts_L.gii` (along with `errts_R.gii`)
  - writes: the results in `errts_trials_target_epoch.RDS` under the same folder as the input

#### 2. Calculate the univariate or multivariate trial-level statistics for each parcel:

Run `code/behav/write_behav_subjs.R` to generate trial-level behavioral data
  - reads: `in/behav/orig/` to CSV files
  - writes: behavioral data csvs to `in/behav`, e.g., `behavior-and-events_wave13_Stroop.csv` .

Run `code/spatial/multivariate_projections_within-task_cv-sessions.R` to generate univariate and multivariate statistics for reliability analysis
  - reads: the behavior CSVs
  - writes: to a csv like `out/spatial/projections__stroop__rda__n_resamples100__demean_run__cv_allsess_wave12.csv`

Run `code/spatial/gather_csvs.R` to recode "wave2/3" and "wave1/3" to "wave1/2" and gather all three CSVs into one
- reads: the behavior CSVs
- writes: into `out/spatial/`
  1. `projections_*` contains trial-level scalar contrast values per univariate and multivariate model
  2. `weights_*` contains vertex-level weights for multivariate decoder
  3. `noise_projs_*` contains noise PCs (vertex-level weights that account for trial-level variability)

#### 3. Fit hierarchical models for the statistics

Run `code/inferential/estimate_reliability.R` to fit the reliability models.
  - reads: `projections_*`
  - writes: the fitted models to `out/inferential/schaefer2018_17_400_fsaverage5/NAMES_OF_ROI/`
  - Warning: fitting the bayesian model is VERY SLOW and usually some jobs fail at the first run!
    - Each model is about 36M and requires about an hour to train (with 4 chains on 4 cores) without further parallelization.
    - With `mclapply()` over the ROIs we can fit 32 models in about 15 hours.
    - Currently, we only focus on "baseline" session.

#### 4. Generate summary statistics and plots

Run `code/inferential/summarize_model_output.R` to extract the fitted parameters of all models
  - reads: the fitted models from step 3.
  - writes: to `out/inferential`, the RDS files `posterior_*` and CSV files `summarystat_*`

Run `code/inferential/plot_manuscript_figures.R` to plot all figures
  - reads: from `out/inferential`, the RDS files `posterior_*` and CSV files `summarystat_*`
  - writes: to `figs`



---

## Relevant Links, Paths (for internal use)

- [test-retest subject list and QC grades](https://3.basecamp.com/3758557/buckets/3792852/messages/4106700214)
- [first-pass reliability analyses](https://3.basecamp.com/3758557/buckets/3792852/messages/3983554628)
- minimally preprocessed BOLD timeseries (TR by vertex)
  - wave 1: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS
  - wave 2: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase3/fMRIPrep_AFNI_ANALYSIS
- trial onset times
  - /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS
