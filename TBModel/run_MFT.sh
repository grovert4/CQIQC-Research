#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account=rrg-aparamek
#SBATCH --time=23:00:00
#SBATCH --output=/scratch/a/aparamek/andykh/SLURMOutputs/slurm-%x-%j.txt
module load CCEnv
module load StdEnv
module load julia/1.9.1
module load intelmpi
mpirun julia --heap-size-hint=3G MFT_wrapper.jl 01.15.2024_Bilayer
