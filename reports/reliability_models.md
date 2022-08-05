---
title: "Test-retest reliability models"
output: 
    html_document:
        toc: false
        toc_depth: 2
        toc_float: false
        number_sections: false
        df_print: paged
---

This document outlines the test-retest reliability models considered in this project.






# Design and Notation


- trial $t \in 1, ..., T$
- condition $c \in \{\text{incongruent}, \text{congruent}\}$
- task $k \in \{\text{axcpt}, \text{stroop}\}$
- session $s \in \{\text{baseline}, \text{proactive}, \text{reactive}\}$
- repetition $r \in \{\text{test}, \text{retest}\}$
- participant $p \in 1, ..., P$
- Response variable: $y$, considered here for a single ROI and for a single type of spatial model (i.e., univariate or multivariate model)

The following indicator variables will be used to denote a \textbf{dummy or treatment-coding} scheme:

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


The following expression will be used to denote a \textbf{contrast-coding} scheme for the condition factor:

$$
  \text{stroop} = 
  \begin{cases}
    \frac{1}{2} & \text{if condition = incongruent} \\
    -\frac{1}{2} & \text{if condition = congruent}
  \end{cases}
$$


# Model 1: single session, Chen (2021)-style random effects

## notation

Following Laird and Ware notation (see [here](https://stefvanbuuren.name/fimd/sec-threeformulations.html) for an example).

$$
\begin{equation}
	\begin{split}
	  & \mathbf{X}_p = (1, \mathrm{retest}, \mathrm{incon}, \mathrm{retest}\cdot\mathrm{incon})_{tcr} \\
	  & \mathbf{Z}_p = (\mathrm{test}, \mathrm{retest}, \mathrm{test} \cdot \mathrm{stroop}, \mathrm{retest} \cdot \mathrm{stroop})_{tcr} \\
		& \mathbf{y}_{p} = \mathbf{X}_p\mathbf{a} + \mathbf{Z}_p \boldsymbol{\alpha}_p + \boldsymbol{\epsilon}_{p} \\
		& \boldsymbol{\alpha}_p \sim \mathcal{N}(\textbf{0}, \textbf{S}) \\
		& \mathbf{S} = 
		  \begin{bmatrix}
        \tau_{11} & \tau_{12} &  &  \\
        \tau_{21} & \tau_{22} &  &  \\
         &  & \tau_{33} & \tau_{34} \\
         &  & \tau_{43} & \tau_{44}
      \end{bmatrix} \\
		& \epsilon_{tcrp} \sim \mathcal{N}(0, \sigma^2) \\
	\end{split}
\end{equation}
$$


Notes:

- covariances between repetition main effect and interaction terms in random effects were set to zero
- in the fixed effects, repetition and condition are treatment coded, with test and congruent as the reference levels
- in the random effects, the condition factor is contrast coded, while the repetition factor is treatment coded


## example code

```r
m1 <- brm(
  y ~ repetition * condition + 
    (0 + mean_test + mean_retest | subj) + 
    (0 + stroop_test + stroop_retest | subj),
    data = d
  )
```

Where `mean_test`, `mean_retest`, `stroop_test`, `stroop_retest` are numeric vectors capturing the mean across 
conditions (incongruent + congruent), and the stroop effect contrast (incongruent - congruent),  per repetition:



```r
cbind(trial = rep(1:10, 2), mean_test, mean_retest, stroop_test, stroop_retest)
```

```
##    trial mean_test mean_retest stroop_test stroop_retest
## 1      1         1           0         0.5           0.0
## 2      2         1           0        -0.5           0.0
## 3      3         1           0         0.5           0.0
## 4      4         1           0        -0.5           0.0
## 5      5         1           0         0.5           0.0
## 6      6         1           0        -0.5           0.0
## 7      7         1           0         0.5           0.0
## 8      8         1           0        -0.5           0.0
## 9      9         1           0         0.5           0.0
## 10    10         1           0        -0.5           0.0
## 11     1         0           1         0.0           0.5
## 12     2         0           1         0.0          -0.5
## 13     3         0           1         0.0           0.5
## 14     4         0           1         0.0          -0.5
## 15     5         0           1         0.0           0.5
## 16     6         0           1         0.0          -0.5
## 17     7         0           1         0.0           0.5
## 18     8         0           1         0.0          -0.5
## 19     9         0           1         0.0           0.5
## 20    10         0           1         0.0          -0.5
```


# Model 2: single session, full model


## notation

Following [Gang Chen's](https://github.com/afni-gangc/afni-gangc.github.io/blob/main/_posts/2022-05-20-blog-post-title-from-file-name.md)
notation.

$$
\begin{equation}
	\begin{split}
		& y_{tcrp} \sim \mathcal{T}(\alpha_{crp}, \nu_{crp}) \\
		& \alpha_{crp} = \mathbf{m} + \boldsymbol{\mu}_p \\
		& \text{log } \nu_{crp} = \boldsymbol{\gamma} + \boldsymbol{\tau}_{p} \\
		& (\boldsymbol{\mu}, \boldsymbol{\tau})_p^{\intercal} \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma})\\
		& \boldsymbol{\Sigma} = 
		  \begin{bmatrix}
        \Sigma_{11} & \dots  & \Sigma_{18}\\
        \vdots & \ddots & \vdots\\
        \Sigma_{81} & \dots  & \Sigma_{88} 
      \end{bmatrix} \\
	\end{split}
\end{equation}
$$

Description:

We assume that the observed responses for a given participant $p$, repetition $r$, condition $c$, and trial $t$, $y_{tcrp}$,
were generated from a Student's t distribution, described by a location parameter $\alpha_{crp}$ and scale parameter $\nu_{crp}$.
Conceptually, these location and scale parameters are analogous to the mean and standard deviation parameters of the 
Gaussian distribution.
The t distribution is used here instead of Gaussian, however, as its fatter tails provide greater robustness to outlying
values in the response variable $y_{tcrp}$.

Next, we decompose the location parameter $\alpha_{crp}$ into group-level fixed effects and subject-level random effects.
We denote the group-level effects with $\mathbf{m}$, a vector that holds the mean response for each condition*repetition 
combination (and is thus of length 4).
We denote the subject-level random effects with $\boldsymbol{\mu}_p$, a set of $P$ vectors that hold the 
deviation of the $p$-th participant's means from the group means in $\mathbf{m}$.

Similar to the location parameter, we decompose the (log of the) scale parameter into group and subject-level effects.
We also denote these effects with vectors, $\boldsymbol{\gamma}$ and $\boldsymbol{\tau}_p$, which are analogous 
to their location-parameter counterparts.

Finally, we assume that each participant --- that is, their eight location and scale parameters, concatenated into a 
vector $(\boldsymbol{\mu}, \boldsymbol{\tau})_p^{\intercal}$ --- were sampled from a single multivariate normal
distribution, with zero mean and covariance matrix $\boldsymbol{\Sigma}$.
This covariance matrix describes the linear relationships among these parameters over subjects.

One of the relationships within $\boldsymbol{\Sigma}$ describes the pivotal quantity of test-retest reliability: 
the correlation in the Stroop contrast between test and retest repetitions, denoted here as $\rho$.
The value of this correlation is obtained through a contrast on $\boldsymbol{\Sigma}$.
Specifically, if $\mathbf{C}$ is a 2-by-8 matrix with columns that implement the Stroop effect contrast (incongruent - congruent) per repetition, $\mathbf{P} = \mathbf{C}\boldsymbol{\Sigma}\mathbf{C}^\intercal$ yields a covariance matrix that,
when normalized, contains the pivotal correlation, $\rho = P_{12}/(\sqrt{P_{11}}\sqrt{P_{22}})$.

Notes:

- 'flat' or dummy-coding with no-intercept model is used to parameterize all factors
- full random effect structure: all covariances estimated
- heterogeneous trial-level variance
- location and scale correlations estimated (hence 8 by 8 and not two 4 by 4s)
- t-distributed error for robustness to outliers


## example code

```r
m2 <- brm(
  bf(
    y ~ 0 + repetition * condition + (0 + repetition * condition | a | subj),
    sigma ~ 0 + repetition * condition + (0 + repetition * condition | a | subj)
    ),
  data = d,
  family = student()
)
```
