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

First we check the outliers in the data by the fitting statistics:

![plot of chunk outliers](figure/outliers-1.png)

There seems to be one outlier for rda, so we exclude it from the future analysis. We will fit another model for it again.


```r
outlier_dat <- dat %>% filter(Term == "elpd_loo", Estimate < -60000)
outlier_dat
```

```
## # A tibble: 1 x 17
##   model response session Term  Estimate Est.Error CI.Lower CI.Upper Q.Lower Q.Upper   MAP  rhat ess_bulk ess_tail Grouping region is_core32
##   <ord> <ord>    <chr>   <chr>    <dbl>     <dbl>    <dbl>    <dbl>   <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl> <chr>    <chr>  <lgl>    
## 1 no_l… rda      baseli… elpd…  -70322.      573.       NA       NA      NA      NA    NA    NA       NA       NA <NA>     17Net… FALSE
```

```r
outliers <- outlier_dat$region
dat <- dat %>% filter(!region %in% .env$outliers)
```

~~Then we compare the Bayesian R-squared of the models (see [here](https://paul-buerkner.github.io/brms/reference/bayes_R2.brmsfit.html)). Each line is for one ROI. The higher the value, the better the fit.~~



~~Next we show the difference in the expected log predictive density (elpd) estimated by leave-one-out (loo) cross-validation across models (see [here](https://mc-stan.org/loo/reference/loo.html)). The higher the value, the better the fit. Here "no_lscov_symm" is used as the base level.~~



~~We can see that "no_lscov_symm" is the best and "fixed_sigma" is much worse than others. Results using [waic](https://mc-stan.org/loo/reference/waic.html) instead of `loo` are almost identical.~~



## Test-Retest Reliability

### Distribution Plots

Here we plot the distribution of TRR over 32 ROIs for each response type and model:

![plot of chunk TRR-hist](figure/TRR-hist-1.png)

Then we plot the distribution of the difference in TRR between different response types for each model.

![plot of chunk TRR-res-diff-hist](figure/TRR-res-diff-hist-1.png)

~~This is a similar plot but for the difference between models for each response type:~~



### Brain Plots

Here we compare the test-retest reliability (TRR) of high-low control demand contrast among response variables and models:

![plot of chunk TRR-brain](figure/TRR-brain-1.png)

Here we plot the difference in TRR between different types of response:

![plot of chunk TRR-res-diff-brain](figure/TRR-res-diff-brain-1.png)

~~and between different models:~~



## Comparison with Previous Results (Skipped)

Here we plot the TRR over core32 parcels and compare it with our previous results. The old method is different from the current one in the following ways:

- Subjects: only 18 instead of 27
- Preprocessing: no divisive normalization for each session
- Response variables:
  - "ridge" uses `mda::fda()`, which is similar to `"rda"` but implemented using ridge regression
  - "rda_full" and "rda_diag" use `sparsediscrim::lda_schafer()` instead of `klaR` to generate the response variable, with either the full covariance matrix (including noise correlation) or only its diagonal (the signal).
  - "uv" is not demeaned for each run
- The TRR model is similar to "fix_sigma" model except that:
  - The model is Gaussian instead of t-distributed
  - The fixed effects are treatment-coded



And the distribution of the difference between response types:


