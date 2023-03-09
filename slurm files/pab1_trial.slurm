#!/bin/bash
#SBATCH --job-name=pab1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_48_test_12_t1.txt 48 > pab1_48_t1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_48_test_12_t2.txt 48 > pab1_48_t2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_48_test_12_t3.txt 48 > pab1_48_t3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_240_test_61_t1.txt 240 > pab1_240_t1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_240_test_61_t2.txt 240 > pab1_240_t2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_240_test_61_t3.txt 240 > pab1_240_t3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_1206_test_302_t1.txt 1206 > pab1_1206_t1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_1206_test_302_t2.txt 1206 > pab1_1206_t2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_1206_test_302_t3.txt 1206 > pab1_1206_t3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_6033_test_1509_t1.txt 6033 > pab1_6033_t1.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_6033_test_1509_t2.txt 6033 > pab1_6033_t2.txt
python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_6033_test_1509_t3.txt 6033 > pab1_6033_t3.txt

python -u  pred_clean_sbatch_bgl3.py pab1_all/pab1_trainbatch_train_30168_test_7542_t1.txt 30168 > pab1_30168_t1.txt
