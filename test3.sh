#!/bin/bash
#SBATCH -J PMOEA
#SBATCH -p COMPUTE -N 1 -n 28

module load /Share/apps/matlab/R2019b
/Share/apps/matlab/R2019b/bin/matlab -nodesktop -nosplash -nodisplay -r "run" -logfile test.log
