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
#SBATCH --output=/scratch/andykh/SLURMOutputs/slurm-%x-%j.txt
MODULES="gcc/12.3 flexiblas openmpi cmake fftw hdf5 boost python/3.10.13 llvm/16 eigen clang"
module purge
module load ${MODULES}


# ===== RAY Configuration =====
# ===== Call your code below =====
python suscep_skyrmion.py 0
python suscep_skyrmion.py 1
python suscep_skyrmion.py 2
python suscep_skyrmion.py 3
python suscep_skyrmion.py 4



