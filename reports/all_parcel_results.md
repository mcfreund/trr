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



Here we plot the TRRs estimated by Pearson correlation and hierarchical Bayesian model (no_lscov_symm).

![plot of chunk scatters-uvmv](figure/scatters-uvmv-1.png)

We fit two ellipses for the distribution over control and non-control regions respectively.

![plot of chunk fit-ellipses](figure/fit-ellipses-1.png)

Then we investigate the improvements by replacing Pearson correaltion with Bayesian model:

<img src="figure/scatter-bias-1.png" title="plot of chunk scatter-bias" alt="plot of chunk scatter-bias" width="50%" style="display: block; margin: auto;" />

We seperate this plot for each network to investiagte which network benefits most from hierarchical Bayesian modeling:

![plot of chunk scatters-bias](figure/scatters-bias-1.png)

We then take a look at whether the improvement by HBM correlates with Pearson(TRR).

![plot of chunk scatter-trr-bias](figure/scatter-trr-bias-1.png)

Then we plot the difference between multivariate and univariate methods for each TRR model (Pearson correlation or HBM):

<img src="figure/scatter-uvmv-diff-1.png" title="plot of chunk scatter-uvmv-diff" alt="plot of chunk scatter-uvmv-diff" width="50%" style="display: block; margin: auto;" />

We seperate this plot for each network to investiagte which network benefits most from multivariate methods:

![plot of chunk scatters-uvmv-diff](figure/scatters-uvmv-diff-1.png)
