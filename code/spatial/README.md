# fitting spatial models

These scripts involve summarizing spatial fMRI activity patterns on each trial into a single score or projection.
Two types of spatial models are considered: univariate and multivariate.
Both of these spatial models are forms of linear dimensionality reduction, in which the trial-by-vertex data **X** are summarized over vertices to form scalar projections for each trial, **y**, with weights (vertex-by-1) **w**:

**y** = **Xw**

In the univariate model, the weights **w** are all positive and uniform across vertices, whereas in the multivariate model, the weights are estimated via the discriminant function of LDA.


## resampling_trials.r

For multivariate decoding (LDA), it is easier to interpret results if the classess are balanced.
Firstly, we should ensure that the classess have equal numbers of trials (so that the discriminant function will be based equally on each class).
We handle this by randomly subsampling trials within each class (without replacement) so that each class has the same nubmer of observations.
Secondly, we should also ensure that other experimental factors that likely impact fMRI activity are balanced across the classes.
In particular, we want to ensure that each class is made up of the same number of stimulus colors, stimulus words, sessions, and runs.
For example, we wouldn't want the "congruent" class to be made up of mostly trials with red and blue stimuli when the "incongruent" class is made up of mostly trials with green and black stimuli.
We handle this by performing the random subsampling of trials in a stratified manner, such that the resampled classes are always balanced in these factors.

The resampling script implements this resampling, and saves indices for the resampled trials within an RDS file.
This RDS file is then read by the script that implements the spatial models.
