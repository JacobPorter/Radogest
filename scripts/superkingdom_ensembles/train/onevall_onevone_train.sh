#!/bin/sh
FEATURES=/scratch/jsp4cu/Data/Superkingdom/full_predictions/onevall_onevone.train.prob
RESPONSE=/scratch/jsp4cu/Data/Superkingdom/All/train/1.all.train.taxid
# Types OneVAll_ensembles  OneVAll_OneVOne_ensembles  OneVOne_ensembles
TYPE=OneVAll_OneVOne_ensembles
FOLDS=4
./train_ensemble.sh LogisticRegression $FOLDS $TYPE $FEATURES $RESPONSE
./train_ensemble.sh RandomForestClassifier $FOLDS $TYPE $FEATURES $RESPONSE
./train_ensemble.sh AdaBoostTree $FOLDS $TYPE $FEATURES $RESPONSE
./train_ensemble.sh DecisionTreeClassifier $FOLDS $TYPE $FEATURES $RESPONSE
