#!/bin/bash

#SBATCH --job-name=symbolic_vars
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB 
#SBATCH --time=48:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/%x-%j.log
#SBATCH --array=3-12

module load matlab

matlab -nodisplay -r "symbolicComputations($SLURM_ARRAY_TASK_ID)"


