#!/bin/bash

#SBATCH -J fit_brms                  ## Name or the job         ##
#SBATCH -o code/bash/fit_brms%j.out  ## Output file name        ##
#SBATCH -e code/bash/fit_brms%j.err  ## Error file name         ##
#SBATCH -N 1                         ## Number of nodes         ##
#SBATCH -n 1                         ## Number of tasks         ##
#SBATCH --cpus-per-task 32           ## Number of CPUs per task ##
#SBATCH -t 23:59:59                  ## Walltime                ##

ncore=`expr ${SLURM_JOB_CPUS_PER_NODE} \* ${SLURM_NNODES}`

module load gcc/10.2.0
module load python
source activate r-env
Rscript code/inferential/estimate_reliability.R uv ${ncore} 17 80