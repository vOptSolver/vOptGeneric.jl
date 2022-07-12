#!/bin/bash
#SBATCH --job-name=try_on_smallExamples

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-VERYSHORT
#SBATCH --output=try_on_smallExamples_OUTPUT.txt
#SBATCH --error=try_on_smallExamples_ERROR.txt

#srun ../../../../julia-1.7.3/bin/julia ./latexWriter.jl

for file in ./*.jl; do
    echo "$file" 
    srun julia "$file"
done

srun julia latexWriter.jl
