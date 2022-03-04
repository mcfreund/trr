#!/usr/bin/env bash

glm_prefix="null_2rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out

subjects_file=in/subjects_wave12_all.txt
sessions=(baseline proactive reactive)
tasks=Stroop#(Axcpt Cuedts Stroop Stern)
waves=(1 2)
do_single_subj=false  ## for dev/debugging

source /data/nil-external/ccp/freund/trr/code/timeseries/$glm_prefix/deconvolve_detrend.sh
