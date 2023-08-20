#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=90
#SBATCH --cpus-per-task=1 
#SBATCH --account=def-aparamek
#SBATCH --time=2:45:00
#SBATCH --mem-per-cpu=4000MB
#SBATCH --output=/scratch/andykh/SLURMOutputs/slurm-%x-%j.txt
module load julia/1.8.5

srun julia SkX_MonoLayer_Run.jl 08.17.2023_Monolayer