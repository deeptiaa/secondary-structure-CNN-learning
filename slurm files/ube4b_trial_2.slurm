#!/bin/bash
#SBATCH --job-name=ube4b
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_14564_test_3642_t1.txt 14564 > ube4b_14564_v1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_14564_test_3642_t2.txt 14564 > ube4b_14564_v2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_14564_test_3642_t3.txt 14564 > ube4b_14564_v3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_2912_test_3642_t1.txt 2912 > ube4b_2912_v1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_2912_test_3642_t2.txt 2912 > ube4b_2912_v2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_2912_test_3642_t3.txt 2912 > ube4b_2912_v3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_582_test_3642_t1.txt 582 > ube4b_582_v1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_582_test_3642_t2.txt 582 > ube4b_582_v2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_582_test_3642_t3.txt 582 > ube4b_582_v3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_116_test_3642_t1.txt 116 > ube4b_116_v1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_116_test_3642_t2.txt 116 > ube4b_116_v2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_116_test_3642_t3.txt 116 > ube4b_116_v3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_23_test_3642_t1.txt 23 > ube4b_23_v1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_23_test_3642_t2.txt 23 > ube4b_23_v2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all_2/ube4b_trainbatch_train_23_test_3642_t3.txt 23 > ube4b_23_v3.txt
