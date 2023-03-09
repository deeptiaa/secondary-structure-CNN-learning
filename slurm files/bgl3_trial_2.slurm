#!/bin/bash
#SBATCH --job-name=bgl3
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_4117_test_5147_t1.txt 4117 > bgl3_5147_v1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_4117_test_5147_t2.txt 4117 > bgl3_5147_v2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_4117_test_5147_t3.txt 4117 > bgl3_5147_v3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_823_test_5147_t1.txt 823 > bgl3_823_v1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_823_test_5147_t2.txt 823 > bgl3_823_v2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_823_test_5147_t3.txt 823 > bgl3_823_v3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_164_test_5147_t1.txt 164 > bgl3_164_v1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_164_test_5147_t2.txt 164 > bgl3_164_v2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_164_test_5147_t3.txt 164 > bgl3_164_v3.txt

python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_32_test_5147_t1.txt 32 > bgl3_32_v1.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_32_test_5147_t2.txt 32 > bgl3_32_v2.txt
python -u  pred_clean_sbatch_bgl3.py bgl3_all_2/bgl3_trainbatch_train_32_test_5147_t3.txt 32 > bgl3_32_v3.txt
