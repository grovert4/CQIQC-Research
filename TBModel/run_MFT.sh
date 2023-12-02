#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=rrg-aparamek
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/a/aparamek/andykh/SLURMOutputs/slurm-%x-%j.txt
module load julia/1.9.1
module load intelmpi
mpirun julia MFT_wrapper.jl 11.27.2023_Bilayer
