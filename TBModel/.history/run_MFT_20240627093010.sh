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
#SBATCH --time=23:45:00
#SBATCH --output=/scratch/a/aparamek/andykh/SLURMOutputs/slurm-%x-%j.txt
module load CCEnv
module load StdEnv/2023
module load julia/1.10.0
module load intel/2023.2.1
module load intelmpi

mpirun julia --project=.. --heap-size-hint=2G MFT_wrapper_n.jl 07.04.2024_Bilayer ./Bilayer_MFT_n.jl
mpirun julia --project=.. --heap-size-hint=2G MFT_wrapper_n.jl 07.02.2024_Bilayer ./Bilayer_MFT_n.jl

#mpirun julia --project=.. --heap-size-hint=2G MFT_wrapper.jl 06.10-3.2024_Bilayer ./Bilayer_MFT_n.jl
#mpirun julia --project=.. --heap-size-hint=2G MFT_wrapper.jl 06.10-4.2024_Bilayer ./Bilayer_MFT_n.jl
mpirun julia --project=.. --heap-size-hint=2G Extract_LastData.jl "/scratch/a/aparamek/andykh/Data/Bilayer_Data" 07.04- Bilayer