Exporting relevant data for Freund, Chen, Chen, & Braver (2024; https://doi.org/10.1101/2024.04.24.591032) to new centralized directory for upload to zenedo.

```bash

path_ruiqi="/data/nil-external/ccp/chenr/trr/out"
path_mike="/data/nil-external/ccp/freund/trr/out"
path_new="/data/nil-external/ccp/freund/trr_data"

## timeseries and timeseries model outputs (vertex by TR):
rsync -avzP --relative ${path_ruiqi}/./timeseries/*/RESULTS/Stroop/baseline* $path_new
rsync -avzP --relative ${path_mike}/./timeseries/*csv $path_new
rsync -avzP --relative ${path_mike}/./timeseries/*h5 $path_new

## spatial model outputs (parcel by trial):
rsync -avzP --relative ${path_mike}/./spatial/*projections__stroop__rda__demean_run*cv_allsess.* $path_new
rsync -avzP --relative ${path_mike}/./spatial/*weights__stroop__rda__demean_run* $path_new
rsync -avzP --relative ${path_mike}/./spatial/*noise_projs__stroop__rda__demean_run* $path_new
rsync -avzP --relative ${path_mike}/./spatial/*noise_projs_mean* $path_new

## inferential (parcel):
rsync -avzP --relative --exclude='*proactive*' --exclude='*reactive*' \
    ${path_ruiqi}/./inferential/schaefer2018_17_400_fsaverage5 $path_new

## zip
zip -rv9 /data/nil-external/ccp/freund/trr/trr_data_export.zip $path_new

## upload to zenodo
#git clone https://github.com/jhpoelen/zenodo-upload.git /data/nil-external/ccp/freund/zenodo-upload
source zenodo.sh  ## export token
/data/nil-external/ccp/freund/zenodo-upload/zenodo_upload.sh 14043319 /data/nil-external/ccp/freund/trr/trr_data_export.zip --verbose

```

to CHPC
NB: ensure directory exists on CHPC

```bash
files="out/spatial/projections__stroop__rda__n_resamples100__demean_run*__cv_allsess.csv"
path_ccplinux1="/data/nil-external/ccp/freund/trr/"
path_chpc="m.freund@login3.chpc.wustl.edu:/home/m.freund/trr/"

## ccplinux1->chpc
rsync -avzP ${path_ccplinux1}${files} ${path_chpc}${files}
```