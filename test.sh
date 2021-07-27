#!/bin/bash
#SBATCH -J E22
#SBATCH -p COMPUTE -N 2 -n 32

cd /Share/home/E19301122/PMOEA
module load /Share/apps/matlab/R2016a
matlab -nodesktop -nosplash -nodisplay -r "run" -logfile test.log
