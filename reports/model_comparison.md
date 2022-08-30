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

Apparently there are two outliers with unexpectedly high Bayes_R2 but unexpectedly low elpd_loo (which is strange - need to check the data). For the sake of consistency, we decided to remove the two regions associated with the outliers from the data:


```r
outlier_dat <- dat %>% filter(Term == "elpd_loo", Estimate < -60000)
outlier_dat
```

```
## # A tibble: 2 x 16
##   model         response session  Term     Estimate Est.Error CI.Lower CI.Upper Q.Lower Q.Upper   MAP  rhat ess_bulk ess_tail Grouping region                     
##   <ord>         <ord>    <chr>    <chr>       <dbl>     <dbl>    <dbl>    <dbl>   <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl> <chr>    <chr>                      
## 1 no_lscov_symm rda      baseline elpd_loo  -65709.      702.       NA       NA      NA      NA    NA    NA       NA       NA <NA>     17Networks_LH_ContA_PFClv_2
## 2 no_lscov_symm rda      baseline elpd_loo  -77219.      556.       NA       NA      NA      NA    NA    NA       NA       NA <NA>     17Networks_RH_ContA_PFCl_3
```

```r
outliers <- outlier_dat$region
dat <- dat %>% filter(!region %in% .env$outliers)
```

Then we compare the Bayesian R-squared of the models (see [here](https://paul-buerkner.github.io/brms/reference/bayes_R2.brmsfit.html)). Each line is for one ROI. The higher the value, the better the fit.


```
## # A tibble: 4 x 2
##   model         mean_R2
##   <ord>           <dbl>
## 1 fixed_sigma   0.00722
## 2 no_lscov_symm 0.00636
## 3 no_lscov      0.00544
## 4 full          0.00524
```

![plot of chunk Bayes-R2](figure/Bayes-R2-1.png)

We can see that generally, fixed_sigma > no_lscov_symm > no_lscov > full.

Next we show the difference in the expected log predictive density (elpd) estimated by leave-one-out (loo) cross-validation across models (see [here](https://mc-stan.org/loo/reference/loo.html)). The higher the value, the better the fit. Here "full" is used as the base level.


```
## [1] "Range of the difference in EPLD_LOO between different models:"
```

```
##                                 Q0          Q5         Q25         Q50         Q75         Q95       Q100
## no_lscov - full         -0.2954626    1.930878    3.140695    4.526943    5.541599    8.118259    8.74665
## no_lscov_symm - full     3.1950237    6.328993    9.974277   13.406579   16.213133   18.726019   21.18368
## fixed_sigma - full   -1091.5676938 -841.875589 -537.705544 -460.460139 -360.815230 -243.393282 -211.90811
```

```
## Warning: Removed 30 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 30 rows containing missing values (geom_point).
```

```
## Warning: Removed 32 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 32 rows containing missing values (geom_point).
```

![plot of chunk diff-elpd-loo](figure/diff-elpd-loo-1.png)

We can see that "no_lscov_symm" is the best and "fixed_sigma" is much worse than others. Results using [waic](https://mc-stan.org/loo/reference/waic.html) instead of `loo` are almost identical.



## Test-Retest Reliability

### Distribution Plots

Here we plot the distribution of TRR over 32 ROIs for each response type and model:

![plot of chunk TRR-hist](figure/TRR-hist-1.png)

Then we plot the distribution of the difference in TRR between different response types for each model.

![plot of chunk TRR-res-diff-hist](figure/TRR-res-diff-hist-1.png)

This is a similar plot but for the difference between models for each response type:

![plot of chunk TRR-mdl-diff-hist](figure/TRR-mdl-diff-hist-1.png)

### Brain Plots

Here we compare the test-retest reliability (TRR) of high-low control demand contrast among response variables and models:

![plot of chunk TRR-brain](figure/TRR-brain-1.png)

Here we plot the difference in TRR between different types of response:

![plot of chunk TRR-res-diff-brain](figure/TRR-res-diff-brain-1.png)

and between different models:

![plot of chunk TRR-mdl-diff-brain](figure/TRR-mdl-diff-brain-1.png)

## Comparison with Previous Results

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

![plot of chunk TRR-hist-core32](figure/TRR-hist-core32-1.png)

And the distribution of the difference between response types:

![plot of chunk TRR-diff-hist-core32](figure/TRR-diff-hist-core32-1.png)
