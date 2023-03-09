#!/bin/bash
#SBATCH --job-name=pab1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_23_test_6_t1.txt 23 > ube4b_23_t1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_23_test_6_t2.txt 23 > ube4b_23_t2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_23_test_6_t3.txt 23 > ube4b_23_t3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_116_test_29_t1.txt 116 > ube4b_116_t1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_116_test_29_t2.txt 116 > ube4b_116_t2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_116_test_29_t3.txt 116 > ube4b_116_t3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_582_test_146_t1.txt 582 > ube4b_582_t1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_582_test_146_t2.txt 582 > ube4b_582_t2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_582_test_146_t3.txt 582 > ube4b_582_t3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_2912_test_729_t1.txt 2912 > ube4b_2912_t1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_2912_test_729_t2.txt 2912 > ube4b_2912_t2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_2912_test_729_t3.txt 2912 > ube4b_2912_t3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_14564_test_3642_t1.txt 14564 > ube4b_14564_t1.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_14564_test_3642_t2.txt 14564 > ube4b_14564_t2.txt
python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_14564_test_3642_t3.txt 14564 > ube4b_14564_t3.txt

python -u  pred_clean_sbatch_bgl3.py ube4b_all/ube4b_trainbatch_train_72824_test_18207_t1.txt 72824 > ube4b_72824_t1.txt
