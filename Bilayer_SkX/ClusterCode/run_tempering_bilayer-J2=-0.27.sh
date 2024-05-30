#!/bin/bash\n#SBATCH --mail-user=tanmay.grover@mail.utoronto.ca  \n#SBATCH --mail-type=ALL\n#SBATCH --mail-type=BEGIN\n#SBATCH --mail-type=END\n#SBATCH --mail-type=FAIL\n#SBATCH --ntasks=20\n#SBATCH --cpus-per-task=1 \n#SBATCH 
--account=def-aparamek\n#SBATCH --time=12:00:00\n#SBATCH --mem-per-cpu=8000MB\n#SBATCH --output=/scratch/grovert4/SLURM/slurm-%x-%j.txt\nmodule load StdEnv/2020\nmodule load julia/1.8.5\n\nsrun julia SkX_Tempering_BiLayer_Run.jl J2=-0.27-30.05.2024-Bilayer-tempering