#!/bin/bash
#SBATCH --job-name=gb1
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_27_test_7_t1.txt 27 > gb1_27_t1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_27_test_7_t2.txt 27 > gb1_27_t2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_27_test_7_t3.txt 27 > gb1_27_t3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_136_test_35_t1.txt 136 > gb1_136_t1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_136_test_35_t2.txt 136 > gb1_136_t2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_136_test_35_t3.txt 136 > gb1_136_t3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_685_test_172_t1.txt 685 > gb1_685_t1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_685_test_172_t1.txt 685 > gb1_685_t2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_685_test_172_t1.txt 685 > gb1_685_t3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_3430_test_858_t1.txt 3430 > gb1_3430_t1.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_3430_test_858_t1.txt 3430 > gb1_3430_t2.txt
python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_3430_test_858_t1.txt 3430 > gb1_3430_t3.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_17154_test_4289_t1.txt 17154 > gb1_17154_t1.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_85772_test_21444_t1.txt 85772 > gb1_85772_t1.txt

python -u  pred_clean_sbatch_bgl3.py gb1_all/gb1_trainbatch_train_428867_test_107217_t1.txt 428867 > gb1_428867_t1.txt
