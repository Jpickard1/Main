#!/bin/bash

#SBATCH --job-name=syntheticData1
#SBATCH --mail-user=jpic@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH --time=1-00:00:00
#SBATCH --account=indikar0
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/testResults/%x-%j.log


module load matlab
matlab -nodisplay -r "driverKronFit('system','GL', 'n0', 2, 'filePath', 'C:\Joshua\MissingData\Projects\AFOSR\kroneckerCalculations\HyperKronFit\syntheticTestGraph1.txt', 'firstPermItrs', 10000, 'gradSamples', 100000, 'maxItrs', 2000, 'outputPath', "/nfs/turbo/umms-indikar/Joshua/Main/Projects/AFOSR/kroneckerCalculations/HyperKronFit/testResults/syntheticData1.m");"

