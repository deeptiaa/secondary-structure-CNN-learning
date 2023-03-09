#!/bin/bash
#SBATCH --job-name=3dqw
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_20_test_6_t1.txt 20 > 3dqw_20_t1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_20_test_6_t2.txt 20 > 3dqw_20_t2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_20_test_6_t3.txt 20 > 3dqw_20_t3.txt

python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_105_test_27_t1.txt 105 > 3dqw_105_t1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_105_test_27_t2.txt 105 > 3dqw_105_t2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_105_test_27_t3.txt 105 > 3dqw_105_t3.txt

python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_530_test_133_t1.txt 530 > 3dqw_530_t1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_530_test_133_t2.txt 530 > 3dqw_530_t2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_530_test_133_t3.txt 530 > 3dqw_530_t3.txt

python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_2652_test_663_t1.txt 2652 > 3dqw_2652_t1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_2652_test_663_t1.txt 2652 > 3dqw_2652_t2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all/3dqw_trainbatch_train_2652_test_663_t1.txt 2652 > 3dqw_2652_t3.txt
