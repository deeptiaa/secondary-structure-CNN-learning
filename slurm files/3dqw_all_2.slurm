#!/bin/bash
#SBATCH --job-name=3dqw
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_530_test_663_t1.txt 530 > 3dqw_530_v1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_530_test_663_t2.txt 530 > 3dqw_530_v2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_530_test_663_t3.txt 530 > 3dqw_530_v3.txt

python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_105_test_663_t1.txt 105 > 3dqw_105_v1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_105_test_663_t2.txt 105 > 3dqw_105_v2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_105_test_663_t3.txt 105 > 3dqw_105_v3.txt

python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_20_test_663_t1.txt 20 > 3dqw_20_v1.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_20_test_663_t2.txt 20 > 3dqw_20_v2.txt
python -u  pred_clean_sbatch_bgl3.py 3dqw_all_2/pro_3dqw_trainbatch_train_20_test_663_t3.txt 20 > 3dqw_20_v3.txt
