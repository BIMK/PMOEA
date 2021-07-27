#!/bin/bash
#SBATCH -J test
#SBATCH -p COMPUTE -N 1 -n 28

cd /Share/home/E19301122/PMOEA_fix
module load /Share/apps/matlab/R2019b
/Share/apps/matlab/R2019b/bin/matlab -nodesktop -nosplash -nodisplay -r "run" -logfile test.log
