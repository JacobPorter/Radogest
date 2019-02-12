#!/bin/sh
nohup srun --partition standard -t 6-23 --exclusive /home/jsp4cu/Applications/Simple_Estimator/simple_estimator.py train -e $1 -i 16 -f $2 -p $2 -v 10 /scratch/jsp4cu/Models/$3/$1/ $4 $5 &> /scratch/jsp4cu/Train_Out/ensembles/ensemble.$3.$1.f_$2.train.out &
