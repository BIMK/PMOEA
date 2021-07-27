#!/bin/bash
#SBATCH -e slurm-%j.out
#SBATCH -o slurm-%j.out
#SBATCH -p COMPUTE
#SBATCH -J test2

module load /Share/apps/matlab/R2019b
/Share/apps/matlab/R2019b/bin/matlab -nodesktop -nosplash -nodisplay -r "subjob;exit"
