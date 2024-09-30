#!/bin/bash

#SBATCH -J fit_brms                     ## Name or the job         ##
#SBATCH --array 1-10                    ## Job array               ##
#SBATCH -o code/bash/fit_brms%A_%a.out  ## Output file name        ##
#SBATCH -e code/bash/fit_brms%A_%a.err  ## Error file name         ##
#SBATCH -N 1                            ## Number of nodes         ##
#SBATCH -n 1                            ## Number of tasks         ##
#SBATCH --cpus-per-task 32              ## Number of CPUs per task ##
#SBATCH -t 10:00:00                     ## Walltime                ##

response_var="uv"
test_ses="reactive"
ncore=`expr ${SLURM_JOB_CPUS_PER_NODE} \* ${SLURM_NNODES}`
roi_start=`expr ${SLURM_ARRAY_TASK_ID} \* 40 - 39`
roi_end=`expr ${roi_start} + 39`

# module load gcc/10.2.0
module load python
source activate r_env
Rscript code/inferential/estimate_reliability.R ${response_var} ${ncore} ${roi_start} ${roi_end} ${test_ses}