#!/usr/bin/env bash

## setup ----

conda deactivate
cd $path_project


if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file    
else  ## for dev
    unset subjects
    subjects=DMCC8033964
    subject=DMCC8033964
    wave=1
    task_i=0
    session_i=0
    hemi_i=0
fi



function deconvolve_2rpm {
	
    ## build xmat
    /usr/local/pkg/afni_18/3dDeconvolve \
    -local_times \
    -force_TR 1.2 \
    -input ${name_img_1}" "${name_img_2} \
    -polort A \
    -float \
    -censor ${dir_stimts}/movregs_FD_mask.txt \
    -num_stimts 0 \
    -ortvec ${dir_stimts}/motion_demean_${sessions[$session_i]}.1D movregs \
    -x1D ${dir_out}/X.xmat.1D \
    -xjpeg ${dir_out}/X.jpg \
    -nobucket \
    -errts ${dir_out}/errts_${hemi}

}


hemis=(L R)
wd=$(pwd)


## run ----


for wave in ${waves[@]}; do

    echo "running wave "${wave}

    glm=$glm_prefix"_wave"$wave

    if [ ${wave} = 1 ]
    then
        wave_dir=HCP_SUBJECTS_BACKUPS
    else  ## wave 2 or greater
        wave_dir=DMCC_Phase$((wave+1))  ## add one to get phase
    fi

    # paths:

    stimts=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/  ## JUST NEED MOVEMENT REGS
    out=$path_project/out/timeseries/
    img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/$wave_dir/fMRIPrep_AFNI_ANALYSIS/
    scripts=$path_project/code/timeseries/

    for subject in ${subjects[@]}; do

        echo ${subject}

        for task_i in ${!tasks[@]}; do

            for session_i in ${!sessions[@]}; do

                ## define paths and names
                sess=${sessions[$session_i]:0:3}  ## get short name
                sess=${sess^}  ## Namecase
                dir_stimts=${stimts}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}
                dir_out=${out}${subject}/RESULTS/${tasks[$task_i]}/${sessions[$session_i]}_${glm}

                ## make result dir and cd into it
                mkdir -p ${dir_out}
                cd ${dir_out}

                for hemi in ${hemis[@]}; do

                    name_img_1=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}1_AP_${hemi}.func.gii
                    name_img_2=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}2_PA_${hemi}.func.gii

                    deconvolve_2rpm < /dev/null > ${dir_out}/runtime_3dDeconvolve.log 2>&1 &

                done

                cd ${wd}  ## back to original dir

           	done

        done

        wait  ## run subjs & waves serially, but each task and session in parallel

    done

done
