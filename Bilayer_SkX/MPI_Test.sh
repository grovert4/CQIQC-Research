#!/bin/bash
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=0:45:00
#SBATCH --output=/scratch/andykh/SLURMOutputs/slurm-%x-%j.txt
module load julia/1.9.1
module load NiaEnv/2022a
srun julia MPITest.jl
