#!/bin/bash
#SBATCH --mail-user=tanmay.grover@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1 
#SBATCH --account=def-aparamek
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000MB
#SBATCH --output=/scratch/grovert4/SLURM/slurm-%x-%j.txt
module load StdEnv/2020
module load julia/1.8.5

srun julia SkX_Tempering_BiLayer_Run.jl J2=-0.08-05.09.2024-Bilayer-MPI-Aminus

