---
title: "Figures for draft of TRR manuscript modified by Ruiqi"
date: "2023-07-10"
output: 
    html_document:
        toc: true
        toc_depth: 3 
        toc_float: true
        number_sections: true
        df_print: paged
        code_folding: hide
---

# Modifed figures for TRR manuscript



## Load data


```
## # A tibble: 400 × 15
##    region uv_po…¹ is_roi  mv_ICC uv_ICC mv_tr…² uv_tr…³ mv_trr…⁴ uv_tr…⁵ mv_tr…⁶
##    <chr>    <dbl> <lgl>    <dbl>  <dbl>   <dbl>   <dbl>    <dbl>   <dbl>   <dbl>
##  1 17Net…    6.40 TRUE    0.438   0.583  0.624    0.935  5.34e-1   0.699  0.153 
##  2 17Net…    6.34 TRUE    0.295   0.534  0.498    0.874  3.99e-1   0.480 -0.0446
##  3 17Net…    5.73 TRUE    0.537   0.397  0.679    0.920  6.07e-1   0.699  0.260 
##  4 17Net…    5.46 TRUE    0.463   0.572  0.634    0.965  5.62e-1   0.838  0.183 
##  5 17Net…    5.45 TRUE    0.712   0.537  0.973    0.961  8.89e-1   0.826  0.675 
##  6 17Net…    5.37 TRUE   -0.0643  0.298  0.0467   0.878 -2.03e-4   0.643 -0.845 
##  7 17Net…    5.13 TRUE    0.184   0.113  0.476    0.423  4.13e-1   0.339 -0.101 
##  8 17Net…    4.98 TRUE    0.365   0.357  0.728    0.849  5.84e-1   0.412  0.0688
##  9 17Net…    4.95 TRUE    0.639   0.621  0.900    0.960  7.80e-1   0.829  0.465 
## 10 17Net…    4.88 TRUE    0.579   0.715  0.846    0.981  7.62e-1   0.919  0.490 
## # … with 390 more rows, 5 more variables: uv_trr_05 <dbl>, mv_trr_sd <dbl>,
## #   uv_trr_sd <dbl>, mv_ratio <dbl>, uv_ratio <dbl>, and abbreviated variable
## #   names ¹​uv_pop_tstat_div, ²​mv_trr_map, ³​uv_trr_map, ⁴​mv_trr_mean,
## #   ⁵​uv_trr_mean, ⁶​mv_trr_05
```

## Contingency table


Table: Contingency table for the three comparisons of interest.

|  MV vs UV TRR   | UV TRR < UV ICC, MV precision < UV precision | UV TRR < UV ICC, MV precision > UV precision | UV TRR > UV ICC, MV precision < UV precision | UV TRR > UV ICC, MV precision > UV precision |
|:---------------:|:--------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|
| MV TRR < UV TRR |                      4                       |                      11                      |                      82                      |                     116                      |
| MV TRR > UV TRR |                      21                      |                      64                      |                      10                      |                      92                      |

There are two results related to the comparison of interest (i.e., the comparison between the multivariate and univariate TRR estimates):

1. When MV TRR is greater than UV TRR, it will also be more precise than UV TRR; similarly if UV TRR is greater it will be more precise than MV too.

2. **In most regions where UV TRR is smaller than UV ICC, MV TRR is greater than UV TRR.** On the contrary, in most regions where UV TRR is greater than UV ICC, MV TRR is smaller than UV TRR.

### Relationship between MV TRR and UV TRR in parcels with positive or negative ICC

<img src="figure/mv-vs-uv-TRR-sgn-uv-ICC-1.png" alt="plot of chunk mv-vs-uv-TRR-sgn-uv-ICC" width="100%" />

When using MAP TRR, MV TRR tended to be higher than UV TRR in parcels with negative UV ICC, but lower than UV TRR in parcels with positive UV ICC. When using mean TRR, MV TRR were still higher than UV TRR in parcels with negative UV ICC, but there was no differences in parcels with positive UV ICC.

## Scatter plots

For uv_MAP_vs_ICC over mv_vs_uv_prcs, separated by the sign of mv_vs_uv_MAP, colored by ROI or non ROI:

<img src="figure/scatter0-1.png" alt="plot of chunk scatter0" width="100%" />

For uv_MAP_vs_ICC over mv_vs_uv_prcs, colored by mv_vs_uv_MAP:

<img src="figure/scatter1-1.png" alt="plot of chunk scatter1" width="100%" />

For uv_MAP_vs_ICC over mv_vs_uv_MAP:

<img src="figure/scatter2-1.png" alt="plot of chunk scatter2" width="100%" />

For mv_vs_uv_prcs over mv_vs_uv_MAP:

<img src="figure/scatter3-1.png" alt="plot of chunk scatter3" width="100%" />

## Relationship between uncertainty, noise ratio, and ICC/TRR

<img src="figure/scatter4-1.png" alt="plot of chunk scatter4" width="100%" />

Obviously, the uncertainty of TRR is positively correlated with the noise ratio (with a Spearman correlation above .9), and negatively correlated with ICC (Spearman correlation around -.6).

### Correlation between statistics

<img src="figure/correlation-1.png" alt="plot of chunk correlation" width="100%" />

### Relationship between mv_TRR_05 and uv_TRR_05

<img src="figure/mv-vs-uv-05-1.png" alt="plot of chunk mv-vs-uv-05" width="100%" />

Obviously, MVPA generates higher lower bound of TRR in most parcels, and more parcels have a lower bound above zero when using MVPA.

## Examine some ROIs


