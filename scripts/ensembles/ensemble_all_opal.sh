#!/bin/sh
# $1 type
nohup ./ensemble_single_opal.sh AdaBoostTree $1 &
nohup ./ensemble_single_opal.sh GradBoost $1 &
nohup ./ensemble_single_opal.sh LogisticRegression $1 &
nohup ./ensemble_single_opal.sh RandomForestClassifier $1 &



