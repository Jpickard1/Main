#!/bin/bash

#SBATCH --job-name=KTEig9
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=25GB
#SBATCH --time=3-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/%x-%j.log


module load matlab
matlab -nodisplay -r "glExpKronTenEig(3,3);"

