#!/bin/bash

#SBATCH -J test_brms                  ## Name or the job         ##
#SBATCH -o code/bash/test_brms%j.out  ## Output file name        ##
#SBATCH -e code/bash/test_brms%j.err  ## Error file name         ##
#SBATCH -N 1                          ## Number of nodes         ##
#SBATCH -n 1                          ## Number of tasks         ##
#SBATCH --cpus-per-task 32            ## Number of CPUs per task ##
#SBATCH -t 4:00:00                    ## Walltime                ##

ncore=`expr ${SLURM_JOB_CPUS_PER_NODE} \* ${SLURM_NNODES}`

module load gcc/10.2.0
module load python
source activate r-env
Rscript code/inferential/estimate_reliability.R uv ${ncore} 1 16