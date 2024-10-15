To deconvolve the HRF associated with Stroop trial onset, we build an FIR design in which each TR post-stimulus onset is modeled as a stick function.
We do this in afni with `'TENTzero(0,16.8,15)'`.

To compare how these stick functions align with the expected (convolved) HRF, we assume a fixed-shape HRF, with `'BLOCK(1,1)'`, as used in Freund Bugg & Braver 2020 JNeurosci.
We assume that this HRF begins at the "effective" onset time of the stimulus.
The effective onset time of the stimulus is the true onset time adjusted to account for slice-timing correction.

The stimulus onset and TR onset were not completely in phase with each other, however their phase relationship was highly consistent across trials. 
Specifically, the true onset time was always `+0.71 s` after the onset of the TR.

Because of slice-timing correction, all slices were aligned to the middle of each TR.
To get the effective stimulus onset time, this adjustment is required to be subtracted: `0.59 s = ((n - 1) * tr) / (2 * n)`, where n = 60 is the number of slices per TR, and tr = 1.2 s.

Therefore, the effective stimulus onset time is `0.12 = 0.71 - 0.59` after the onset of the nearest TR.


```bash
bash
conda deactivate
path="code/timeseries/hrf"

## build the FIR model and the fixed-shape model for comparison:
3dDeconvolve -nodata 26 1.2 -polort -1 -num_stimts 3   \
	-stim_times 1 '1D: 0' 'TENTzero(0,16.8,15)' -stim_label 1 TENTzero \
    -stim_times 2 '1D: 0.12' 'BLOCK(1,1)' -stim_label 2 block \
	-stim_times 3 '1D: 0.12' 'SPMG1(1)' -stim_label 3 spm \
	-x1D "${path}/x_trial.1D"
```

To build expected HRFs for consecutive trials, we use the `BLOCK` again.

Following each trial, was an inter-trial interval of either 1, 2, of 3 TRs (1.2, 2.4, or 3.6 s).
After the ITI was a 0.4 second warning sign (flicker of fixation).
Following the warning sign, the stimulus itself was presented for 2 s.

So, if a stimulus onset occurred at time = 0, the soonest a stimulus could be presented would be t = 2 + 0.4 + 1.2 s, or 3 TRs.

Here we build 3 trials that are as close in time as possible (i.e., only 3 TRs interveining).

```bash
## build three sequential trials:
3dDeconvolve -nodata 26 1.2 -polort 0 -num_stimts 3   \
	-stim_times 1 '1D: 0' 'BLOCK(1,1)' -stim_label 1 trial1 \
	-stim_times 2 '1D: 3.6' 'BLOCK(1,1)' -stim_label 2 trial2 \
	-stim_times 3 '1D: 7.2' 'BLOCK(1,1)' -stim_label 3 trial3 \
	-x1D "${path}/x_3trials.1D"
```


Now I just read and write these xmats into R, to save them in an easily-readable csv format.

```R
path <- "code/timeseries/hrf"
names <- c("x_trial.1D", "x_3trials.1D")
file_names <- file.path(path, names)
lapply(
	file_names,
	function(nm) {
		x <- mikeutils::read_xmat(nm)
		write.csv(x, paste0(nm, ".csv"))
	}
)
```