#!/usr/bin/env bash

path_project=/data/nil-external/ccp/freund/trr

glm_prefix="null_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out
subjects_file=in/subjects_wave12_all.txt
sessions=(baseline proactive reactive)
tasks=(Axcpt Cuedts Stroop Stern)
waves=(1 2)
do_single_subj=false  ## for dev/debugging

source $path_project/code/timeseries/$glm_prefix/deconvolve_detrend.sh
