#!/bin/bash
for i in {0..29}
do
   ./build/ES grn10 lsoda sade 1050000 $i >> sade-output5-new-limits.txt &
done

exit

