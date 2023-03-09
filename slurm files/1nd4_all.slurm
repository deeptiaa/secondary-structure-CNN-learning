#!/bin/bash
#SBATCH --job-name=1nd4
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_32_test_8_t1.txt 32 > 1nd4_32_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_32_test_8_t2.txt 32 > 1nd4_32_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_32_test_8_t3.txt 32 > 1nd4_32_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_162_test_41_t1.txt 162 > 1nd4_162_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_162_test_41_t2.txt 162 > 1nd4_162_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_162_test_41_t3.txt 162 > 1nd4_162_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_815_test_204_t1.txt 815 > 1nd4_815_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_815_test_204_t2.txt 815 > 1nd4_815_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_815_test_204_t3.txt 815 > 1nd4_815_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_4076_test_1019_t1.txt 4076 > 1nd4_4076_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_4076_test_1019_t2.txt 4076 > 1nd4_4076_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all/1nd4_trainbatch_train_4076_test_1019_t3.txt 4076 > 1nd4_4076_t3.txt
