#!/bin/bash

#SBATCH -J fit_brms                     ## Name or the job         ##
#SBATCH --array 1                       ## Job array               ##
#SBATCH -o .logs/fit_brms%A_%a.out      ## Output file name        ##
#SBATCH -e .logs/fit_brms%A_%a.err      ## Error file name         ##
#SBATCH -N 1                            ## Number of nodes         ##
#SBATCH -n 1                            ## Number of tasks         ##
#SBATCH --cpus-per-task 32              ## Number of CPUs per task ##
#SBATCH -t 10:00:00                     ## Walltime                ##
#SBATCH --account=todd_braver
#SBATCH --partition=tier2_cpu

response_var="uv"
test_ses="baseline"
roi_col="network"
ncore=`expr ${SLURM_JOB_CPUS_PER_NODE} \* ${SLURM_NNODES}`
if [ "$roi_col" = "parcel" ]; then
  roi_start=`expr ${SLURM_ARRAY_TASK_ID} \* 40 - 39`
  roi_end=`expr ${roi_start} + 39`
elif [ "$roi_col" = "network" ]; then
  roi_start=1
  roi_end=17
fi

#module load gcc/10.2.0
#module load python
#source activate r_env
module load R
Rscript code/inferential/estimate_reliability.R ${response_var} ${ncore} ${roi_start} ${roi_end} ${test_ses} ${roi_col}