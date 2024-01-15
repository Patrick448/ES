#!/bin/bash

BUDGET=1050000
INPUT_FILE="DadosZeEduardo_DATA.txt"
MODEL="grn5ncyc"
EXE_PATH="./build/ES"
DIFF="norm"
EXP_NAME="$MODEL-$BUDGET-$DIFF"
SET_DIV="-sd 0 113 114 161"

mkdir $EXP_NAME
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a cmaes -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-cmaes.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a de -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-de.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a sade -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-sade.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a es-i -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-es-i.txt"
./batch_run.sh 4 30 "$EXE_PATH -i $INPUT_FILE -m $MODEL -e lsoda -a es-ni -n $BUDGET -s \$i $SET_DIV >> $EXP_NAME/output-$MODEL-es-ni.txt"