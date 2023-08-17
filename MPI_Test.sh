#!/bin/bash
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=0:45:00
#SBATCH --output=/scratch/andykh/SLURMOutputs/slurm-%x-%j.txt
module load julia/1.8.5
srun julia MPITest.jl