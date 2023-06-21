#!/bin/bash

#SBATCH --job-name=loadDawn
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=0-01:30:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/BigFit/adjLists/%x-%j.log


module load matlab
matlab -nodisplay -r "loadDAWN('GL',3); loadDAWN('GL',4); loadDAWN('GL',5);"
# matlab -nodisplay -r "driverHipHop2HG('GL', $SLURM_ARRAY_TASK_ID, 0.01, '/nfs/turbo/umms-indikar/Joshua/Main/Code/reproductions/Hip-Hop/HiP-HoP_Pax6_FullConformations/', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/')"

