#!/bin/bash

#SBATCH --job-name=HG2THETA
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=0-01:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
# SBATCH --array=0:39
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HG2theta/%x-%j.log


module load matlab
matlab -nodisplay -r "driverHG2theta('worker',1,'epsilon',0.01,'system','GL', 'n0', 2, 'filePath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/', 'firstPermItrs', 10000, 'gradSamples', 100000, 'maxItrs', 75, 'learningRate', 1e-8', 'outputPath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HG2theta/');"

