#!/bin/bash
#SBATCH -J E22
#SBATCH -p COMPUTE
#SBATCH -n 28
#SBATCH -N 1

cd /Share/home/E19301122/PMOEA
/Share/apps/matlab/R2016a/bin/matlab -nodisplay -r "run" -logfile lym.log
