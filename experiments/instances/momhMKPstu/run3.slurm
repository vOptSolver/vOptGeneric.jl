#!/bin/bash
#SBATCH --job-name=largeC

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=largeC_ERROR.txt
#SBATCH --error=largeC_OUTPUT.txt


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp3.jl "$file"
done