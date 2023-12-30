#!/bin/bash

BUDGET=1050000
INPUT_FILE="GRN10.txt"
MODEL="grn10"
EXE_PATH="./build/ES"
EXP_NAME="$MODEL-$BUDGET"

mkdir $EXP_NAME
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a cmaes -n $BUDGET -s \$i >> $EXP_NAME/output-$MODEL-cmaes.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a de -n $BUDGET -s \$i >> $EXP_NAME/output-$MODEL-de.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a es-i -n $BUDGET -s \$i >> $EXP_NAME/output-$MODEL-es-i.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a es-ni -n $BUDGET -s \$i >> $EXP_NAME/output-$MODEL-es-ni.txt"