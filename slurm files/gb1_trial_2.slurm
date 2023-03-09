#!/bin/bash
#SBATCH --job-name=gb1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_85772_test_21444_t1.txt 85772 > gb1_85772_v1.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_17154_test_21444_t1.txt 17154 > gb1_17154_v1.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_3430_test_21444_t1.txt 3430 > gb1_3430_v1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_3430_test_21444_t2.txt 3430 > gb1_3430_v2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_3430_test_21444_t3.txt 3430 > gb1_3430_v3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_685_test_21444_t1.txt 685 > gb1_685_v1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_685_test_21444_t2.txt 685 > gb1_685_v2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_685_test_21444_t3.txt 685 > gb1_685_v3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_136_test_21444_t1.txt 136 > gb1_136_v1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_136_test_21444_t2.txt 136 > gb1_136_v2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_136_test_21444_t3.txt 136 > gb1_136_v3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_27_test_21444_t1.txt 27 > gb1_27_v1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_27_test_21444_t2.txt 27 > gb1_27_v2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all_2/gb1_trainbatch_train_27_test_21444_t3.txt 27 > gb1_27_v3.txt
