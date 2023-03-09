#!/bin/bash
#SBATCH --job-name=1be9
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_265_test_332_t1.txt 265 > 1be9_265_v1.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_265_test_332_t2.txt 265 > 1be9_265_v2.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_265_test_332_t3.txt 265 > 1be9_265_v3.txt

python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_52_test_332_t1.txt 52 > 1be9_52_v1.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_52_test_332_t2.txt 52 > 1be9_52_v2.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all_2/pro_1be9_trainbatch_train_52_test_332_t3.txt 52 > 1be9_52_v3.txt

