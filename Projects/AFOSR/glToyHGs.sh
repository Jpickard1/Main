#!/bin/bash

#SBATCH --job-name=symToyHG2
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB 
#SBATCH --time=7-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/toyHG/%x-%j.log
#SBATCH --array=1-3

module load matlab

matlab -nodisplay -r "toy_hypergraphs_sym($SLURM_ARRAY_TASK_ID)"

# matlab -nodisplay -r "symbolicComputations($SLURM_ARRAY_TASK_ID)"


