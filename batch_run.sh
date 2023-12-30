#!/bin/bash

CONT=0
BATCH_SIZE=$1
MAX_RUNS=$2
echo "Command received: $3"

start=$(date +%s)

while [ $CONT -lt $MAX_RUNS ]
do
  END=$(($CONT+$BATCH_SIZE-1))
  if [ "$END" -eq "$MAX_RUNS" ] || [ "$END" -gt "$MAX_RUNS" ]; then
    END=$(($MAX_RUNS-1))
  fi

  echo "running batch $CONT to $END"
  for ((i=CONT;i<=END;i++))
  do
   #./build/ES -i GRN5.txt -m grn5 -e lsoda -a cmaes -n 10000 -s $i >> bash-test.txt &
   eval $3 &
   echo $i
  done

  echo "waiting for processes..."
  ps -U $USER -o pid,psr,comm | grep  "ES"
  wait
  CONT=$(($CONT+$BATCH_SIZE))
done

echo "finished"
end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"

exit

