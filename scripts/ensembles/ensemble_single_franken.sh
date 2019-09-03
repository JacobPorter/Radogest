#!/bin/sh
# $1 estimator
srun -A fungcat --partition pegasus_q -t 10-1 --exclusive ~/Applications/Simple_Estimator/simple_estimator.py train -i 8 -f 4 -p 4 -v 8 --estimator $1 /scratch/fungcat/jsporter/Models/Ensembles_CD/franken_ensemble/$1/ /scratch/fungcat/jsporter/FrankenEnsemble/ensemble_train/ensemble_training.cd.10mil.plinko_ova_ovo.opal_ova_ovo.prob  /scratch/fungcat/jsporter/Data/Data_Ensemble_CDs/1.10mil.taxid &> /scratch/fungcat/jsporter/Train_Out/Ensembles_CD/ensemble_training.franken.$1.i8.f4.p4.train.out &
