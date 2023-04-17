#!/bin/bash
#SBATCH --mail-user=spmoran@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --job-name=MATLAB1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16g
#SBATCH --time=48:00:00
#SBATCH --account=indikar1
#SBATCH --partition=standard

module load matlab/R2019b
matlab -r "run('/home/spmoran/indikar/temp-sean/CellTracking_Current/Colorassign_Hungarian_NoDisplayCellTracker.m');exit;"
#matlab -nosplash -nodesktop -nodisplay -r "run('/home/spmoran/indikar/temp-sean/scripts/cellTracker.m');exit;"
