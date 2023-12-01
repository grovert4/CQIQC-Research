#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1 
#SBATCH --account=rrg-aparamek
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/andykh/SLURMOutputs/slurm-%x-%j.txt
module load julia/1.9.1
module load NiaEnv/2022a
mpirun julia MPITest.jl
