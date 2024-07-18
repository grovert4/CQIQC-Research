#!/bin/bash
#SBATCH --mail-user=andrew.hardy@mail.utoronto.ca  
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --account=rrg-aparamek
#SBATCH --time=23:45:00
#SBATCH --output=/scratch/a/aparamek/andykh/SLURMOutputs/slurm-%x-%j.txt
module load CCEnv
module load StdEnv/2023
module load gcc/12 openmpi
module load python/3.11 imkl cmake clang mpi4py hdf5 boost fftw
source ~/triqsvenv/bin/activate
julia --project=../../TBModel ../examples/models/SkX_0.jl
julia --project=../../TBModel ../examples/interactions/SkX_Coulomb_0.jl
julia --project=../../TBModel run_RPA.jl --input ../Inputs/07.17-0.2024_SkX.yml --run_bare true
# julia --project=../../TBModel ../examples/models/SkX_2.jl
# julia --project=../../TBModel ../examples/interactions/SkX_Coulomb_2.jl
# julia --project=../../TBModel run_RPA.jl --input ../Inputs/07.17-2.2024_SkX.yml --run_bare true
# julia --project=../../TBModel ../examples/models/SkX_4.jl
# julia --project=../../TBModel ../examples/interactions/SkX_Coulomb_4.jl
# julia --project=../../TBModel run_RPA.jl --input ../Inputs/07.17-4.2024_SkX.yml --run_bare true