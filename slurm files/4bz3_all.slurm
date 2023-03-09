#!/bin/bash
#SBATCH --job-name=4bz3
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_28_test_8_t1.txt 28 > 4bz3_28_t1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_28_test_8_t2.txt 28 > 4bz3_28_t2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_28_test_8_t3.txt 28 > 4bz3_28_t3.txt

python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_145_test_37_t1.txt 145 > 4bz3_145_t1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_145_test_37_t2.txt 145 > 4bz3_145_t2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_145_test_37_t3.txt 145 > 4bz3_145_t3.txt

python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_728_test_182_t1.txt 728 > 4bz3_728_t1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_728_test_182_t2.txt 728 > 4bz3_728_t2.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_728_test_182_t3.txt 728 > 4bz3_728_t3.txt

python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_3643_test_911_t1.txt 3643 > 4bz3_3643_t1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_3643_test_911_t1.txt 3643 > 4bz3_3643_t1.txt
python -u  pred_clean_sbatch_bgl3.py 4bz3_all/4bz3_trainbatch_train_3643_test_911_t1.txt 3643 > 4bz3_3643_t1.txt

