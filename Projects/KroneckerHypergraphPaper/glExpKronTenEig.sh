#!/bin/bash

#SBATCH --job-name=KTEig
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50GB
#SBATCH --time=3-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --array=0-29
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/%x-%j.log


module load matlab
matlab -nodisplay -r "glEXpKronTenEig(2,2);"

