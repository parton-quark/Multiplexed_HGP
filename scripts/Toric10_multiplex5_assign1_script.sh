#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=5_1000_1_0_101_1_Toric10
#SBATCH --output=output_5_1000_1_0_101_1_Toric10.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 5 1000 1 0 101 1 ../input_matrices/Toric10
