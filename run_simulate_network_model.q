#!/bin/bash
#
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=4096
/opt/apps/MATLAB/R2012b/bin/matlab -nojvm -nodisplay -r "simulate_network_model($SLURM_ARRAY_TASK_ID,1,1,1)"



