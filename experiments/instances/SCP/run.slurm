#!/bin/bash
#SBATCH --job-name=testSCP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testSCP_OUTPUT.txt
#SBATCH --error=testSCP_ERROR.txt


for file in ./BOSCP/*.dat; do
    echo "$file"
    srun julia vOptBOSCP.jl "$file"
done

srun julia latexWriter.jl