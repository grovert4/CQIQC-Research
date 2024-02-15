#!/bin/bash
#SBATCH --mail-user=tanmay.grover@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=150
#SBATCH --cpus-per-task=1 
#SBATCH --account=def-aparamek
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=4000MB
#SBATCH --output=/scratch/grovert4/SLURM/slurm-%x-%j.txt
module load julia/1.8.5

srun julia SkX_BiLayer_Run.jl 14.02.2024-Bilayer-constantfield
