#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=4_100000_0.6_0_61_5_N512_K8_symmetric
#SBATCH --output=output_4_100000_0.6_0_61_5_N512_K8_symmetric.txt

source /apps/unit/NemotoU/NicholasSoftware/bin/activate
module load python/3.11.4
python ../multiplexing_VH_decoder.py 4 100000 0.6 0 61 5 ../input_matrices/N512_K8_symmetric
