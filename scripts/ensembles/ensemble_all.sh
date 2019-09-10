#!/bin/sh
# $1 type
nohup ./ensemble_single.sh AdaBoostTree $1 &
nohup ./ensemble_single.sh GradBoost $1 &
nohup ./ensemble_single.sh LogisticRegression $1 &
nohup ./ensemble_single.sh RandomForestClassifier $1 &



