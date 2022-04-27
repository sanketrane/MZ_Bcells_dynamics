#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%A.%a.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%A.%a.err_out
#SBATCH --get-user-env
#SBATCH -J test
#SBATCH -D /opt/mesh/eigg/sanket/MZ_New_dynamics
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr3506@cumc.columbia.edu
#SBATCH --array=1-10
#SBATCH --cpus-per-task=25

srun Rscript --vanilla scripts/ki67_SHM_MZ.R
