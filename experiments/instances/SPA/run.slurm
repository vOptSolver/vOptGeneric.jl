#!/bin/bash
#SBATCH --job-name=testSPA

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testSPA_OUTPUT.txt
#SBATCH --error=testSPA_ERROR.txt


for file in ./BOSPA/*.txt; do
    echo "$file"
    srun julia vOptBOSPA.jl "$file"
done

srun julia latexWriter.jl