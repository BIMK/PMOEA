#!/bin/bash
#SBATCH -e slurm-%j.out
#SBATCH -o slurm-%j.out
#SBATCH -N 1 -c 2
#SBATCH -p COMPUTE
#SBATCH -J PMOEA

module load /Share/apps/matlab/R2019b
/Share/apps/matlab/R2019b/bin/matlab -nodesktop -nosplash -nodisplay -r "run_main"
