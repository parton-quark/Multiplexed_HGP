#!/bin/bash
#SBATCH -p compute
#SBATCH --time=99:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_100000_1_0_101_3_N512_K8_symmetric
#SBATCH --output=output_2_100000_1_0_101_3_N512_K8_symmetric.txt
python ../multiplexing_VH_decoder.py 2 100000 1 0 101 3 N512_K8_symmetric
