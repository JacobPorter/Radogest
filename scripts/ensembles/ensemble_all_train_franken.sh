#!/bin/sh
nohup ./ensemble_single_franken.sh AdaBoostTree &
nohup ./ensemble_single_franken.sh GradBoost &
nohup ./ensemble_single_franken.sh LogisticRegression &
nohup ./ensemble_single_franken.sh RandomForestClassifier &



