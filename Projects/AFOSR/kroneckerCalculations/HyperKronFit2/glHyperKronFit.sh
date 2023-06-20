#!/bin/bash

#SBATCH --job-name=GO1!
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=1-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit2/testResults/%x-%j.log


module load matlab
matlab -nodisplay -r "driverHyperKronFit('system','GL', 'n0', 2, 'filePath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/as20graph.txt', 'firstPermItrs', 10000, 'gradSamples', 100000, 'maxItrs', 500, 'learningRate', 1e-9', 'outputPath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit2/testResults/as20graph_out.m');"

