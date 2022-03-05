# trr

Test--Retest Reliability of fMRI measures in Dual Mechanisms of Cognitive Control dataset

## repo organization

- all analysis scripts in `code/`
- scripts in `code` write to `out/`
- `in/` contains behavioral data files, subject lists
- files in `in/` are not modified by any scripts

## key factors in DMCC dataset

- trialtype: "high demand" (difficult trials), "low demand" (easy trials)
- run: run1, run2
- task: Axcpt, Cuedts, Stern, Stroop
- session: baseline, proactive, reactive
- wave: wave1, wave2 (wave1 = test, wave2 = retest)

## relevant links, paths

- [test-retest subject list and QC grades](https://3.basecamp.com/3758557/buckets/3792852/messages/4106700214)
- [first-pass reliability analyses](https://3.basecamp.com/3758557/buckets/3792852/messages/3983554628#__recording_4046665465)
- minimally preprocessed BOLD timeseries (TR by vertex)
  - wave 1: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS
  - wave 2: /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase3/fMRIPrep_AFNI_ANALYSIS
- trial onset times
  - /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/EVTS