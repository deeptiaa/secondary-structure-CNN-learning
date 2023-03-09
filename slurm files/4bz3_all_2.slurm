#!/bin/bash
#SBATCH --job-name=4bz3
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_728_test_911_t1.txt 728 > 4bz3_728_v1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_728_test_911_t2.txt 728 > 4bz3_728_v2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_728_test_911_t3.txt 728 > 4bz3_728_v3.txt

python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_145_test_911_t1.txt 145 > 4bz3_145_v1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_145_test_911_t2.txt 145 > 4bz3_145_v2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_145_test_911_t3.txt 145 > 4bz3_145_v3.txt

python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_28_test_911_t1.txt 28 > 4bz3_28_v1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_28_test_911_t2.txt 28 > 4bz3_28_v2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all_2/pro_4bz3_trainbatch_train_28_test_911_t3.txt 28 > 4bz3_28_v3.txt

