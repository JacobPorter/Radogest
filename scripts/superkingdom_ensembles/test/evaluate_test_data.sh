#!/bin/sh
# $1 Type
# $2 Classifier
# $3 Features
# $4 Response
# 
nohup srun -t 6-23 --partition standard /home/jsp4cu/Applications/Simple_Estimator/simple_estimator.py evaluate /scratch/jsp4cu/Models/$1_ensembles/$2/ $3 $4 &> /scratch/jsp4cu/Prediction_Out/$1.$2.evaluate.out &
