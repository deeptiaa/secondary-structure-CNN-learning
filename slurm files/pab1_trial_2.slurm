#!/bin/bash
#SBATCH --job-name=pab1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_6033_test_7542_t1.txt 6033 > pab1_6033_v1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_6033_test_7542_t2.txt 6033 > pab1_6033_v2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_6033_test_7542_t3.txt 6033 > pab1_6033_v3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_1206_test_7542_t1.txt 1206 > pab1_1206_v1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_1206_test_7542_t2.txt 1206 > pab1_1206_v2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_1206_test_7542_t3.txt 1206 > pab1_1206_v3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_240_test_7542_t1.txt 240 > pab1_240_v1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_240_test_7542_t2.txt 240 > pab1_240_v2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_240_test_7542_t3.txt 240 > pab1_240_v3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_48_test_7542_t1.txt 48 > pab1_48_v1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_48_test_7542_t2.txt 48 > pab1_48_v2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all_3/pab1_trainbatch_train_48_test_7542_t3.txt 48 > pab1_48_v3.txt

