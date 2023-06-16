---
title: "Test-retest reliability model results for all parcels"
output: 
    html_document:
        toc: false
        toc_depth: 2
        toc_float: false
        number_sections: false
        df_print: paged
---



(New) contingency table:


```
## # A tibble: 6 × 7
##   region                      ICC_rda ICC_uv HBM_MAP_rda HBM_MAP_uv uncertainty_rda uncertainty_uv
##   <chr>                         <dbl>  <dbl>       <dbl>      <dbl>           <dbl>          <dbl>
## 1 17Networks_LH_ContA_Cingm_1 -0.101   0.293     -0.0852      0.210           0.565          0.540
## 2 17Networks_LH_ContA_IPS_1    0.531   0.592      0.612       0.864           0.193          0.144
## 3 17Networks_LH_ContA_IPS_2   -0.0428  0.482      0.462       0.631           0.423          0.352
## 4 17Networks_LH_ContA_IPS_3    0.537   0.397      0.607       0.699           0.192          0.261
## 5 17Networks_LH_ContA_IPS_4    0.639   0.621      0.780       0.829           0.160          0.177
## 6 17Networks_LH_ContA_IPS_5    0.712   0.537      0.889       0.826           0.107          0.180
```

Highly reliable regions:

![plot of chunk high-trr](figure/high-trr-1.png)

We take a look at the fixed effect:

![plot of chunk scatters-t](figure/scatters-t-1.png)

Here we plot the TRRs estimated by Pearson correlation and hierarchical Bayesian model (no_lscov_symm).

![plot of chunk scatters-uvmv](figure/scatters-uvmv-1.png)

We fit two ellipses for the distribution over control and non-control regions respectively.

![plot of chunk fit-ellipses](figure/fit-ellipses-1.png)

Then we investigate the improvements by replacing Pearson correaltion with Bayesian model:

<img src="figure/scatter-bias-1.png" alt="plot of chunk scatter-bias" width="50%" style="display: block; margin: auto;" />

We seperate this plot for each network to investiagte which network benefits most from hierarchical Bayesian modeling:

![plot of chunk scatters-bias](figure/scatters-bias-1.png)

We then take a look at whether the improvement by HBM correlates with Pearson(TRR).

![plot of chunk scatter-trr-bias](figure/scatter-trr-bias-1.png)

Then we plot the difference between multivariate and univariate methods for each TRR model (Pearson correlation or HBM):

<img src="figure/scatter-uvmv-diff-1.png" alt="plot of chunk scatter-uvmv-diff" width="50%" style="display: block; margin: auto;" />

We seperate this plot for each network to investiagte which network benefits most from multivariate methods:

![plot of chunk scatters-uvmv-diff](figure/scatters-uvmv-diff-1.png)
