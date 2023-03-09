#!/bin/bash
#SBATCH --job-name=bgl3
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_32_test_9_t1.txt 32 > bgl3_32_t1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_32_test_9_t2.txt 32 > bgl3_32_t2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_32_test_9_t3.txt 32 > bgl3_32_t3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_164_test_41_t1.txt 164 > bgl3_164_t1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_164_test_41_t2.txt 164 > bgl3_164_t2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_164_test_41_t3.txt 164 > bgl3_164_t3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_823_test_206_t1.txt 823 > bgl3_823_t1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_823_test_206_t2.txt 823 > bgl3_823_t2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_823_test_206_t3.txt 823 > bgl3_823_t3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_4117_test_1030_t1.txt 4117 > bgl3_4117_t1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_4117_test_1030_t2.txt 4117 > bgl3_4117_t2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_4117_test_1030_t3.txt 4117 > bgl3_4117_t3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all/bgl3_trainbatch_train_20590_test_5147_t1.txt 20590 > bgl3_20590_t1.txt
