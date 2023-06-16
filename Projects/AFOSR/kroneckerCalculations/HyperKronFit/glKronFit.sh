#!/bin/bash

#SBATCH --job-name=syntheticData4
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=7-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/testResults/%x-%j.log


module load matlab
matlab -nodisplay -r "driverKronFit('system','DBTM', 'n0', 2, 'filePath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/syntheticTestGraph1.txt', 'firstPermItrs', 1000, 'gradSamples', 2500, 'maxItrs', 5, 'learningRate', 1e-6', 'outputPath', '/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/testResults/syntheticData4.m');"

