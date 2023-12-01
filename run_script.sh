#!/bin/bash
for i in {0..29}
do
   ./build/ES grn5 lsoda sade 1050000 $i >> sade-output5-old-limits.txt &
done

exit

