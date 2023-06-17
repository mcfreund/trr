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
## # A tibble: 8 × 4
## # Groups:   mv_HBM_vs_uv_HBM, uv_HBM_vs_uv_ICC [4]
##   mv_HBM_vs_uv_HBM uv_HBM_vs_uv_ICC uv_uncty_vs_mv_uncty     n
##              <dbl>            <dbl>                <dbl> <int>
## 1               -1               -1                   -1     6
## 2               -1               -1                    1    11
## 3               -1                1                   -1    83
## 4               -1                1                    1   118
## 5                1               -1                   -1    19
## 6                1               -1                    1    61
## 7                1                1                   -1     9
## 8                1                1                    1    93
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
