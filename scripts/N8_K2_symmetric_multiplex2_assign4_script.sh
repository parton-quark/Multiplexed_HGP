#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_100000_1_0_101_4_N8_K2_symmetric
#SBATCH --output=output_2_100000_1_0_101_4_N8_K2_symmetric.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 2 100000 1 0 101 4 ../input_matrices/N8_K2_symmetric
