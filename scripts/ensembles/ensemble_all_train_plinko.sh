#!/bin/sh
# $1 type
nohup ./ensemble_single_train_plinko.sh AdaBoostTree $1 &
nohup ./ensemble_single_train_plinko.sh GradBoost $1 &
nohup ./ensemble_single_train_plinko.sh LogisticRegression $1 &
nohup ./ensemble_single_train_plinko.sh RandomForestClassifier $1 &



