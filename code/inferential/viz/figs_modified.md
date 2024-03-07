---
title: "Figures for draft of TRR manuscript modified by Ruiqi"
date: "2023-07-19"
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

1. In any cases, MV precision tends to be better than UC precision. However, in parcels where MV TRR < UV TRR, MV precision is more likely to be lower than UV precision.

2. **In most regions where UV TRR is smaller than UV ICC, MV TRR is greater than UV TRR.** On the contrary, in most regions where UV TRR is greater than UV ICC, MV TRR is smaller than UV TRR.

### "Weird" parcels with TRR < ICC

Raw TRR values from UV and MV were shown below, separated by the sign of uv_TRR - uv_ICC:

<div class="figure">
<img src="figure/mv-vs-uv-TRR-sgn-uv-ICC-1.png" alt="plot of chunk mv-vs-uv-TRR-sgn-uv-ICC" width="100%" />
<p class="caption">plot of chunk mv-vs-uv-TRR-sgn-uv-ICC</p>
</div>

Obviously, these "weird" parcels where TRR < ICC usually have very low reliability, sometimes even negative. As we will show below, they also have very high uncertainty:

<div class="figure">
<img src="figure/weird-parcels-1.png" alt="plot of chunk weird-parcels" width="100%" />
<p class="caption">plot of chunk weird-parcels</p>
</div>

We can see a negative correlation between TRR-ICC and TRR uncertainty, and a positive correlation between TRR-ICC and ICC itself. In other words, those "weird" parcels with TRR < ICC usually has low reliability with very high uncertainty.

Next, we examine whether MVPA helps reduce uncertainty and improve reliability in these parcels:

<div class="figure">
<img src="figure/weird-parcels-2-1.png" alt="plot of chunk weird-parcels-2" width="100%" />
<p class="caption">plot of chunk weird-parcels-2</p>
</div>

We can see that, in "weird" parcels, MVPA did help improve both the reliability and the precision of TRR estimates. In normal parcels, MVPA did not improve the reliability (in fact it was lower), but still improved the precision.

This plot also shows the other type of "weird" parcels, i.e., those with higher uncertainty (less precision) using MV. From this plot we can see that these parcels already had a very high and precise TRR estimate using UV (though some of this kind of parcels still get even more precise estimate with MVPA).

## Relationship between uncertainty, noise ratio, and ICC/TRR

<div class="figure">
<img src="figure/scatter4-1.png" alt="plot of chunk scatter4" width="100%" />
<p class="caption">plot of chunk scatter4</p>
</div>

Obviously, the uncertainty of TRR is positively correlated with the noise ratio (with a Spearman correlation above .9), and negatively correlated with ICC (Spearman correlation around -.6).

Therefore, it is interesting to examine the parcels where MVPA can and cannot reduce noise ratio:

<div class="figure">
<img src="figure/uv_vs_mv_ratio_brain-1.png" alt="plot of chunk uv_vs_mv_ratio_brain" width="100%" />
<p class="caption">plot of chunk uv_vs_mv_ratio_brain</p>
</div>

```
## [1] "Parcels where UV ratio is lowest compare to MV (i.e., weird ones):"
```

```
## # A tibble: 400 × 2
##    region                               uv_vs_mv_ratio
##    <chr>                                         <dbl>
##  1 17Networks_RH_LimbicA_TempPole_4             -0.945
##  2 17Networks_RH_SomMotA_5                      -0.770
##  3 17Networks_RH_VisCent_ExStr_2                -0.619
##  4 17Networks_RH_LimbicB_OFC_6                  -0.574
##  5 17Networks_LH_LimbicB_OFC_5                  -0.501
##  6 17Networks_RH_TempPar_5                      -0.479
##  7 17Networks_RH_SomMotB_Aud_1                  -0.442
##  8 17Networks_LH_ContC_pCun_1                   -0.437
##  9 17Networks_LH_SomMotA_5                      -0.430
## 10 17Networks_LH_VisPeri_ExStrSup_5             -0.420
## 11 17Networks_LH_LimbicA_TempPole_4             -0.362
## 12 17Networks_RH_SalVentAttnA_ParOper_3         -0.349
## 13 17Networks_RH_ContC_Cingp_2                  -0.331
## 14 17Networks_LH_LimbicA_TempPole_5             -0.313
## 15 17Networks_LH_VisPeri_ExStrInf_2             -0.305
## 16 17Networks_RH_SomMotB_S2_6                   -0.303
## 17 17Networks_RH_DorsAttnB_FEF_1                -0.301
## 18 17Networks_LH_DefaultB_Temp_1                -0.296
## 19 17Networks_LH_LimbicA_TempPole_1             -0.271
## 20 17Networks_RH_VisPeri_ExStrInf_3             -0.255
## 21 17Networks_RH_SomMotA_9                      -0.252
## 22 17Networks_LH_DefaultA_pCunPCC_4             -0.251
## 23 17Networks_RH_SomMotB_Ins_1                  -0.249
## 24 17Networks_LH_SomMotA_11                     -0.247
## 25 17Networks_LH_DefaultA_PFCm_1                -0.244
## 26 17Networks_LH_ContB_Temp_1                   -0.244
## 27 17Networks_RH_SalVentAttnA_FrOper_1          -0.243
## 28 17Networks_LH_SomMotA_19                     -0.239
## 29 17Networks_LH_LimbicA_TempPole_2             -0.231
## 30 17Networks_LH_VisPeri_ExStrSup_2             -0.230
## 31 17Networks_RH_SomMotA_10                     -0.229
## 32 17Networks_RH_SalVentAttnA_ParMed_2          -0.215
## 33 17Networks_RH_VisPeri_ExStrInf_5             -0.211
## 34 17Networks_LH_SomMotB_Cent_3                 -0.206
## 35 17Networks_LH_SalVentAttnB_Ins_2             -0.201
## 36 17Networks_RH_SomMotA_19                     -0.195
## 37 17Networks_RH_SomMotA_14                     -0.193
## 38 17Networks_RH_DorsAttnA_TempOcc_1            -0.191
## 39 17Networks_LH_VisCent_ExStr_1                -0.190
## 40 17Networks_RH_LimbicB_OFC_2                  -0.190
## # … with 360 more rows
```

### Correlation between statistics

<div class="figure">
<img src="figure/correlation-1.png" alt="plot of chunk correlation" width="100%" />
<p class="caption">plot of chunk correlation</p>
</div>

### Relationship between mv_TRR_05 and uv_TRR_05

<div class="figure">
<img src="figure/mv-vs-uv-05-1.png" alt="plot of chunk mv-vs-uv-05" width="100%" />
<p class="caption">plot of chunk mv-vs-uv-05</p>
</div>

Obviously, MVPA generates higher lower bound of TRR in most parcels, and more parcels have a lower bound above zero when using MVPA.

## Examine some ROIs



## More scatter plots

For uv_MAP_vs_ICC over mv_vs_uv_prcs, separated by the sign of mv_vs_uv_MAP, colored by ROI or non ROI:



For uv_MAP_vs_ICC over mv_vs_uv_prcs, colored by mv_vs_uv_MAP:



For uv_MAP_vs_ICC over mv_vs_uv_MAP:

<div class="figure">
<img src="figure/scatter2-1.png" alt="plot of chunk scatter2" width="100%" />
<p class="caption">plot of chunk scatter2</p>
</div>

We can see a strong negative correlation. In other words, the more gain from replacing ICC with HBM, the less gain from further replacing univariate with multivariate.

For mv_vs_uv_prcs over mv_vs_uv_MAP:

<div class="figure">
<img src="figure/scatter3-1.png" alt="plot of chunk scatter3" width="100%" />
<p class="caption">plot of chunk scatter3</p>
</div>

We can see a positive correlation between the gain of TRR from MV and the gain of precision from MV, along with the negative correlation mentioned above. It seems like that the parcels where MV precision is much worse than UV all have positive UV ICC with a further boost from UV HBM. In these parcels, MV is worse both in TRR and in precision.

## Calculate some numbers needed in the results section


```
## # A tibble: 2 × 3
##   is_roi mean_mv_ICC mean_uv_ICC
##   <lgl>        <dbl>       <dbl>
## 1 FALSE        0.186       0.146
## 2 TRUE         0.336       0.354
```

```
## # A tibble: 1 × 4
##   n_mv_ICC_gt_07 n_uv_ICC_gt_07 p_mv_ICC_gt_07 p_uv_ICC_gt_07
##            <int>          <int>          <dbl>          <dbl>
## 1              5              3           1.25           0.75
```

```
## # A tibble: 2 × 3
##   is_roi mean_mv_TRR mean_uv_TRR
##   <lgl>        <dbl>       <dbl>
## 1 FALSE        0.387       0.346
## 2 TRUE         0.563       0.762
```

```
## # A tibble: 2 × 3
##   is_roi mv_TRR_vs_ICC uv_TRR_vs_ICC
##   <lgl>          <dbl>         <dbl>
## 1 FALSE          0.201         0.201
## 2 TRUE           0.227         0.408
```

```
## # A tibble: 1 × 4
##   n_mv_trr_gt_07 n_uv_trr_gt_07 p_mv_trr_gt_07 p_uv_trr_gt_07
##            <int>          <int>          <dbl>          <dbl>
## 1            165            186           41.2           46.5
```

```
## # A tibble: 2 × 3
##   is_roi mean_mv_precision mean_uv_precision
##   <lgl>              <dbl>             <dbl>
## 1 FALSE               2.70              2.19
## 2 TRUE                4.09              3.21
```

```
## # A tibble: 2 × 5
##   is_roi n_mv_more_precise n_uv_more_precise p_mv_more_precise p_uv_more_precise
##   <lgl>              <int>             <int>             <dbl>             <dbl>
## 1 FALSE                256               104              71.1              28.9
## 2 TRUE                  27                13              67.5              32.5
```

```
## # A tibble: 1 × 2
##   n_mv_precise_and_less_noise p_mv_precise_and_less_noise
##                         <int>                       <dbl>
## 1                         259                        64.8
```
