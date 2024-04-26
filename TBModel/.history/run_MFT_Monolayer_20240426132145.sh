#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --account=rrg-aparamek
#SBATCH --time=23:45:00
#SBATCH --output=/scratch/a/aparamek/andykh/SLURMOutputs/slurm-%x-%j.txt
module load CCEnv
module load StdEnv/2023
module load julia/1.10.0
module load intel/2023.2.1
module load intelmpi
mpirun julia --project=.. --heap-size-hint=1G MFT_wrapper.jl 04.26.2024_Monolayer_NN ./Monolayer_MFT.jl
#mpirun julia --project=.. --heap-size-hint=1G Extract_LastData.jl "/scratch/a/aparamek/andykh/Data/Monolayer_Data" 04.19
#mpirun julia --project=.. --heap-size-hint=1G MFT_wrapper.jl 04.20.2024_Monolayer_NN ./Monolayer_MFT.jl
mpirun julia --project=.. --heap-size-hint=1G Extract_LastData.jl "/scratch/a/aparamek/andykh/Data/Monolayer_Data" 04.20

