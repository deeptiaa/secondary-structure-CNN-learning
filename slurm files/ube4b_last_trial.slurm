#!/bin/bash
#SBATCH --job-name=pab1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_72824_test_18207_t1.txt 72824 > ube4b_72824_t1.txt
