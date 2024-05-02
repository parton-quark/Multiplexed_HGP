#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=8_100000_1_0_101_5_N128_K2_symmetric
#SBATCH --output=output_8_100000_1_0_101_5_N128_K2_symmetric.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 8 100000 1 0 101 5 ../input_matrices/N128_K2_symmetric
