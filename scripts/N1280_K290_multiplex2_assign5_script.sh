#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_10000_1_0_101_5_N1280_K290
#SBATCH --output=output_2_10000_1_0_101_5_N1280_K290.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 2 10000 1 0 101 5 ../input_matrices/N1280_K290