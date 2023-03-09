#!/bin/bash
#SBATCH --job-name=1d5r
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_36_test_10_t1.txt 36 > 1d5r_36_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_36_test_10_t2.txt 36 > 1d5r_36_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_36_test_10_t3.txt 36 > 1d5r_36_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_186_test_47_t1.txt 186 > 1d5r_186_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_186_test_47_t2.txt 186 > 1d5r_186_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_186_test_47_t3.txt 186 > 1d5r_186_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_932_test_234_t1.txt 932 > 1d5r_932_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_932_test_234_t2.txt 932 > 1d5r_932_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_932_test_234_t3.txt 932 > 1d5r_932_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_4666_test_1167_t1.txt 4666 > 1d5r_4666_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_4666_test_1167_t2.txt 4666 > 1d5r_4666_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1d5r_all/1d5r_trainbatch_train_4666_test_1167_t3.txt 4666 > 1d5r_4666_t3.txt

