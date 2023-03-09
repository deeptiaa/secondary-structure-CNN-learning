#!/bin/bash
#SBATCH --job-name=1be9
#SBATCH --partition=compute
#SBATCH --account=pfaendtner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --chdir=/gscratch/pfaendtner/danoop/

python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_52_test_14_t1.txt 52 > 1be9_52_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_52_test_14_t2.txt 52 > 1be9_52_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_52_test_14_t3.txt 52 > 1be9_52_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_265_test_67_t1.txt 265 > 1be9_265_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_265_test_67_t2.txt 265 > 1be9_265_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_265_test_67_t3.txt 265 > 1be9_265_t3.txt

python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_1328_test_332_t1.txt 1328 > 1be9_1328_t1.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_1328_test_332_t2.txt 1328 > 1be9_1328_t2.txt
python -u  pred_clean_sbatch_bgl3.py 1be9_all/1be9_trainbatch_train_1328_test_332_t3.txt 1328 > 1be9_1328_t3.txt

