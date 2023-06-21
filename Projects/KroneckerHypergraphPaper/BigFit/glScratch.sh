#!/bin/bash

#SBATCH --job-name=SC-EnronFit
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=25GB
#SBATCH --time=0-10:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/BigFit/%x-%j.log


module load matlab
matlab -nodisplay -r "scratchEnron('GL');"
# matlab -nodisplay -r "driverHipHop2HG('GL', $SLURM_ARRAY_TASK_ID, 0.01, '/nfs/turbo/umms-indikar/Joshua/Main/Code/reproductions/Hip-Hop/HiP-HoP_Pax6_FullConformations/', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/')"

