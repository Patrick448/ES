#!/bin/bash
for i in {0..29}
do
   ./build/ES -i GRN5.txt -m grn5 -e lsoda -a cmaes -n 1050000 -s $i >> testing-grn5-lsoda.txt &
done

exit

