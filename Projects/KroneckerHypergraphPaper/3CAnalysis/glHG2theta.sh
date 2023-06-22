#!/bin/bash

#SBATCH --job-name=HG2THETA
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=3-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --array=0-29
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HG2theta/%x-%j.log


module load matlab
matlab -nodisplay -r "driverHG2theta('worker',$SLURM_ARRAY_TASK_ID,'epsilon',0.01,'system','GL', 'n0', 5, 'filePath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HipHop2HG/', 'firstPermItrs', 10000, 'gradSamples', 100000, 'maxItrs', 100, 'learningRate', 1e-11', 'outputPath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/KroneckerHypergraphPaper/3CAnalysis/HG2theta/');"

