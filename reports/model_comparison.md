---
title: "Model Comparison"
output: 
  html_document:
    toc: false
    toc_depth: 2
    toc_float: false
    number_sections: false
    df_print: paged
---



# Model Comparison

## Fitting Statistics

First we compare the Bayesian R-squared of the models (see [here](https://paul-buerkner.github.io/brms/reference/bayes_R2.brmsfit.html)). Each line is for one ROI. The higher the value, the better the fit.

![plot of chunk Bayes-R2](figure/Bayes-R2-1.png)

```
## # A tibble: 4 x 2
##   model         mean_R2
##   <ord>           <dbl>
## 1 fixed_sigma   0.0108 
## 2 no_lscov_symm 0.0103 
## 3 no_lscov      0.00982
## 4 full          0.00951
```

We can see that generally, fixed_sigma > no_lscov_symm > no_lscov > full.

Next we show the difference in the expected log predictive density (elpd) estimated by leave-one-out (loo) cross-validation across models (see [here](https://mc-stan.org/loo/reference/loo.html)). The higher the value, the better the fit. Here "full" is used as the base level.

![plot of chunk diff-elpd-loo](figure/diff-elpd-loo-1.png)

```
## [1] "Range of the difference in elpd_loo between fixed_sigma and full models: (-1208.55845425573, -45.5684478145631)"
```

We can see that "no_lscov_symm" is the best and "fixed_sigma" is much worse than others. Results using [waic](https://mc-stan.org/loo/reference/waic.html) instead of `loo` are almost identical.



## Test-Retest Reliability

### Brain Plots

Here we compare the test-retest reliability (TRR) of high-low control demand contrast among response variables and models:

![plot of chunk TRR-brain](figure/TRR-brain-1.png)

Here we plot the difference in TRR between different types of response:

![plot of chunk TRR-res-diff-brain](figure/TRR-res-diff-brain-1.png)

and between different models:

![plot of chunk TRR-mdl-diff-brain](figure/TRR-mdl-diff-brain-1.png)

### Distribution Plots

Here we plot the distribution of TRR over 32 ROIs for each response type and model:

![plot of chunk TRR-hist](figure/TRR-hist-1.png)

Then we plot the distribution of the difference in TRR between different response types for each model.

![plot of chunk TRR-res-diff-hist](figure/TRR-res-diff-hist-1.png)

This is a similar plot but for the difference between models for each response type:

![plot of chunk TRR-mdl-diff-hist](figure/TRR-mdl-diff-hist-1.png)

## Comparison with Previous Results

Here we plot the TRR over core32 parcels and compare it with our previous results. The old method is different from the current one in the following ways:

- Subjects: only 18 instead of 27
- Preprocessing: no divisive normalization for each session
- Response variables:
  - "ridge" is the same as the current one, using `mda::fda()`("ridge")
  - "rda_full" and "rda_diag" use `sparsediscrim::lda_schafer()` instead of `klaR()` to generate the response variable, with either the full covariance matrix (including noise correlation) or only its diagonal (the signal).
  - "uv" is not centered nor scaled
- The TRR model is similar to "fix_sigma" model except that:
  - The model is Gaussian instead of t-distributed
  - The fixed effects are treatment-coded

![plot of chunk TRR-hist-core32](figure/TRR-hist-core32-1.png)

And the distribution of the difference between response types:

![plot of chunk TRR-diff-hist-core32](figure/TRR-diff-hist-core32-1.png)
