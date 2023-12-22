#!/bin/bash
for i in {0..5}
do
   ./build/ES grn5 lsoda cmaes 1050000 $i >> testing-grn5-lsoda.txt &
done

exit

