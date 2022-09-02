#!/bin/bash
#SBATCH --job-name=testSPP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testSPP_OUTPUT.txt
#SBATCH --error=testSPP_ERROR.txt


for subfolder in ./BOSPP/*; do
    for file in $subfolder/*; do
        echo "$file"
        srun julia vOptBOSPP.jl "$file"
    done
done
