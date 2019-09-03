#!/bin/sh
# $1 estimator
# $2 type: ova, ovo, ova_ovo
srun -A fungcat --partition pegasus_q -t 7-1 --exclusive ~/Applications/Simple_Estimator/simple_estimator.py train -i 8 -f 4 -p 4 -v 8 --estimator $1 /scratch/fungcat/jsporter/Models/Ensembles_CD/opal_ensemble/$2/$1/ /scratch/fungcat/jsporter/Opal/coding_domain/ensemble_train/opal.ensemble_train.$2.prob /scratch/fungcat/jsporter/Data/Data_Ensemble_CDs/1.10mil.taxid &> /scratch/fungcat/jsporter/Train_Out/Ensembles_CD/ensemble_training.opal.$1.$2.i8.f4.p4.train.out &
