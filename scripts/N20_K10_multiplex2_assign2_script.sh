#!/bin/bash
#SBATCH -p compute
#SBATCH --time=99:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_100000_1_0_101_2_N20_K10
#SBATCH --output=output_2_100000_1_0_101_2_N20_K10.txt
python ../multiplexing_VH_decoder.py 2 100000 1 0 101 2 N20_K10
