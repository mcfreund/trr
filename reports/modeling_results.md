---
title: "Test-retest reliability model results"
output: 
    html_document:
        toc: false
        toc_depth: 2
        toc_float: false
        number_sections: false
        df_print: paged
---



# Relilability Modeling Results

This document outlines the test-retest reliability (TRR) models considered in this project, and relavant results. The TRRs are calculated in "baseline" session only.

## Design and Notation


- trial $t \in 1, ..., T$
- condition $c \in \{\text{incongruent}, \text{congruent}\}$
- task $k \in \{\text{axcpt}, \text{stroop}\}$
- session $s \in \{\text{baseline}, \text{proactive}, \text{reactive}\}$
- repetition $r \in \{\text{test}, \text{retest}\}$
- participant $p \in 1, ..., P$
- Response variable: $y$, considered here for a single ROI and for a single type of spatial model (i.e., univariate or multivariate model)

The following indicator variables will be used to denote a $\textbf{dummy or treatment-coding}$ scheme:

$$
\begin{equation*}
\begin{split}
	&\text{incon} = 
	\begin{cases}
		1 & \text{if condition = incongruent} \\
		0 & \text{if condition = congruent}
	\end{cases},\\
	&\text{retest} = 
	\begin{cases}
		1 & \text{if repetition = retest} \\
		0 & \text{if repetition = test}
	\end{cases},\\
	&\text{test} = 1 - \text{retest}
\end{split}
\end{equation*}
$$


The following expression will be used to denote a $\textbf{contrast-coding}$ scheme for the condition factor:

$$
  \text{stroop} = 
  \begin{cases}
    \frac{1}{2} & \text{if condition = incongruent} \\
    -\frac{1}{2} & \text{if condition = congruent}
  \end{cases}
$$

## Model 0: Pearson correlations

$$
z_{rp} = {\langle y_{tcrp} \rangle}_{t,\ c = \text{Incongruent}} - {\langle y_{tcrp} \rangle}_{t,\ c = \text{Congruent}} \\
\text{TRR} = \text{Cor}(z_{1p}, z_{2p})
$$

First we plot the TRR estimated from univariate or multivariate methods for each region:


```
## [1] "TRRs estimated from Pearson correlation:"
```

```
## # A tibble: 400 x 4
##    region                          rda     uv is_core32
##    <chr>                         <dbl>  <dbl> <lgl>    
##  1 17Networks_LH_ContA_Cingm_1 -0.101  0.293  FALSE    
##  2 17Networks_LH_ContA_IPS_1    0.531  0.592  FALSE    
##  3 17Networks_LH_ContA_IPS_2   -0.0428 0.482  FALSE    
##  4 17Networks_LH_ContA_IPS_3    0.537  0.397  FALSE    
##  5 17Networks_LH_ContA_IPS_4    0.639  0.621  FALSE    
##  6 17Networks_LH_ContA_IPS_5    0.712  0.537  FALSE    
##  7 17Networks_LH_ContA_PFCd_1   0.534  0.0817 TRUE     
##  8 17Networks_LH_ContA_PFCl_1   0.432  0.564  TRUE     
##  9 17Networks_LH_ContA_PFCl_2   0.438  0.583  TRUE     
## 10 17Networks_LH_ContA_PFCl_3   0.295  0.534  TRUE     
## # … with 390 more rows
```

![plot of chunk Pearson-scatter](figure/Pearson-scatter-1.png)

And the distribution of the difference between the result of "rda" and "uv":

![plot of chunk Pearson-diff](figure/Pearson-diff-1.png)

## Model 1: “fixed_sigma"

This model is t-distributed, where the sigma term is independent from the conditions:


```r
formula_string <-
  paste0(
    " ~ ",
    "0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 + ",
    "(0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)"
  )
formula_sigma <- formula(sigma ~ 1)
```

We made the same plots, but this time we only include core32 regions, and show the mean of the data by a red dot (same below).


```
## [1] "TRRs estimated from 'fixed_sigma' model:"
```

```
## # A tibble: 32 x 3
##    region                           uv    rda
##    <chr>                         <dbl>  <dbl>
##  1 17Networks_LH_ContA_PFCd_1   0.135  0.786 
##  2 17Networks_LH_ContA_PFCl_1   0.772  0.655 
##  3 17Networks_LH_ContA_PFCl_2   0.727  0.539 
##  4 17Networks_LH_ContA_PFCl_3   0.680  0.342 
##  5 17Networks_LH_ContA_PFClv_2  0.461  0.610 
##  6 17Networks_LH_ContB_IPL_2    0.339  0.646 
##  7 17Networks_LH_ContB_PFClv_1  0.719  0.604 
##  8 17Networks_LH_ContB_PFClv_2 -0.0692 0.205 
##  9 17Networks_LH_ContB_PFClv_3  0.0971 0.567 
## 10 17Networks_LH_ContC_Cingp_2  0.0554 0.0171
## # … with 22 more rows
```

![plot of chunk fixed-sigma-scatter](figure/fixed-sigma-scatter-1.png)

Here we compare the results if we normalize the data within each vertex before fitting the spatial model.

![plot of chunk w-wo-divnorm-fs](figure/w-wo-divnorm-fs-1.png)

We can see that divisive normalization seems to reduce TRR for "uv" but did not make big difference for "rda". Therefore, divisive normalization might be biased towards multivariate methods without improving the overall result, and we didn't use it in all other analysis in this report.


## Model 2: "no_lscov_symm"

This model is similar to model 1 but the sigma term is no longer independent from the conditions. Instead, it is also predicted with a same formula as that for the mean:


```r
formula_string <-
  paste0(
    " ~ ",
    "0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 + ",
    "(0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)"
  )
formula_sigma <-
  formula(
    sigma ~ 0 + mean_wave1 + mean_wave2 + hilo_wave1 + hilo_wave2 +
    (0 + mean_wave1 + mean_wave2 | subj) + (0 + hilo_wave1 + hilo_wave2 | subj)
  )
```

The same scatter plot:


```
## [1] "TRRs estimated from 'no_lscov_symm' model:"
```

```
## # A tibble: 32 x 3
##    region                          uv    rda
##    <chr>                        <dbl>  <dbl>
##  1 17Networks_LH_ContA_PFCd_1  0.312  0.766 
##  2 17Networks_LH_ContA_PFCl_1  0.804  0.708 
##  3 17Networks_LH_ContA_PFCl_2  0.721  0.565 
##  4 17Networks_LH_ContA_PFCl_3  0.500  0.419 
##  5 17Networks_LH_ContA_PFClv_2 0.353  0.597 
##  6 17Networks_LH_ContB_IPL_2   0.462  0.661 
##  7 17Networks_LH_ContB_PFClv_1 0.667  0.686 
##  8 17Networks_LH_ContB_PFClv_2 0.0592 0.0994
##  9 17Networks_LH_ContB_PFClv_3 0.0671 0.512 
## 10 17Networks_LH_ContC_Cingp_2 0.114  0.0840
## # … with 22 more rows
```

![plot of chunk no-lscov-symm-scatter](figure/no-lscov-symm-scatter-1.png)

Violin plot:

![plot of chunk no-ls-cov-symm-violin](figure/no-ls-cov-symm-violin-1.png)

## Comparison

![plot of chunk comparison](figure/comparison-1.png)
