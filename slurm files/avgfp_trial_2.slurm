#!/bin/bash
#SBATCH --job-name=avgfp
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_8264_test_10342_t1.txt 8264 > avgfp_8264_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_8264_test_10342_t2.txt 8264 > avgfp_8264_v2.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_8264_test_10342_t3.txt 8264 > avgfp_8264_v3.txt

python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_331_test_10342_t1.txt 331 > avgfp_331_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_331_test_10342_t2.txt 331 > avgfp_331_v2.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_331_test_10342_t3.txt 331 > avgfp_331_v3.txt

python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_66_test_10342_t1.txt 66 > avgfp_66_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_66_test_10342_t2.txt 66 > avgfp_66_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_66_test_10342_t3.txt 66 > avgfp_66_v1.txt

python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_13_test_10342_t1.txt 13 > avgfp_13_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_13_test_10342_t2.txt 13 > avgfp_13_v1.txt
python -u  pred_clean_sbatch_bgl3.py avgfp_all_3/avgfp_trainbatch_train_13_test_10342_t3.txt 13 > avgfp_13_v1.txt

