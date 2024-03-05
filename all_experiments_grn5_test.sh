#!/bin/bash

BUDGET=1050000
INPUT_FILE="GRN5_DATA.txt"
MODEL="grn5new"
EXE_PATH="./build/ES"
DIFF="norm_v4"
EXP_NAME="$MODEL-$BUDGET-$DIFF"
SET_DIV="-sd 0 34 35 49"
TEST_SET=""

mkdir $EXP_NAME
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE $TEST_SET -m $MODEL -e lsoda -a cmaes -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-cmaes.csv"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE $TEST_SET -m $MODEL -e lsoda -a de -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-de.csv"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE $TEST_SET -m $MODEL -e lsoda -a sade -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-sade.csv"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE $TEST_SET -m $MODEL -e lsoda -a es-i -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-es-i.csv"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE $TEST_SET -m $MODEL -e lsoda -a es-ni -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-es-ni.csv"