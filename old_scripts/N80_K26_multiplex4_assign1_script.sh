#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=4_100000_1_0_101_1_N80_K26
#SBATCH --output=output_4_100000_1_0_101_1_N80_K26.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 4 100000 1 0 101 1 ../input_matrices/N80_K26
