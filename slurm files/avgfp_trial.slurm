#!/bin/bash
#SBATCH --job-name=avgfp
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=05:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py avgfp_all/avgfp_trainbatch_train_41372_test_10342_t1.txt 41372 > avgfp_41372_t1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all/avgfp_trainbatch_train_41372_test_10342_t2.txt 41372 > avgfp_41372_t2.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all/avgfp_trainbatch_train_41372_test_10342_t3.txt 41372 > avgfp_41372_t3.txt
