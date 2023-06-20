#!/bin/bash

#SBATCH --job-name=HH2HG
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=1-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=array
#SBATCH --array=0-29
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/%x-%j.log


module load matlab
matlab -nodisplay -r "driverHipHop2HG('GL', $SLURM_ARRAY_TASK_ID, 0.01, '/nfs/turbo/umms-indikar/Joshua/Main/Code/reproductions/Hip-Hop/HiP-HoP_Pax6_FullConformations/', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/')"

