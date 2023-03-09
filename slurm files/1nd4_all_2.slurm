#!/bin/bash
#SBATCH --job-name=1nd4
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_815_test_1019_t1.txt 815 > 1nd4_815_v1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_815_test_1019_t2.txt 815 > 1nd4_815_v2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_815_test_1019_t3.txt 815 > 1nd4_815_v3.txt

python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_162_test_1019_t1.txt 162 > 1nd4_162_v1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_162_test_1019_t2.txt 162 > 1nd4_162_v2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_162_test_1019_t3.txt 162 > 1nd4_162_v3.txt

python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_32_test_1019_t1.txt 32 > 1nd4_32_v1.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_32_test_1019_t2.txt 32 > 1nd4_32_v2.txt
python -u  pred_clean_sbatch_bgl3.py 1nd4_all_2/pro_1nd4_trainbatch_train_32_test_1019_t3.txt 32 > 1nd4_32_v3.txt
