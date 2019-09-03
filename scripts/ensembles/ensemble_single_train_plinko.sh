#!/bin/sh
# $1 estimator
# $2 type: ova, ovo, ova_ovo
srun -A fungcat --partition pegasus_q -t 7-1 --exclusive ~/Applications/Simple_Estimator/simple_estimator.py train -i 8 -f 4 -p 4 -v 8 --estimator $1 /scratch/fungcat/jsporter/Models/Ensembles_CD/train_plinko/$2/$1/ /scratch/fungcat/jsporter/Data/Data_Ensemble_CDs/train_plinko/ensemble_training_data.train_plinko.cd.10mil.$2.prob  /scratch/fungcat/jsporter/Data/Data_Ensemble_CDs/train_plinko/1.train_plinko.10mil.taxid &> /scratch/fungcat/jsporter/Train_Out/Ensembles_CD/ensemble_training.train_plinko.$1.$2.i8.f4.p4.train.out &
