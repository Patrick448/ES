#!/bin/bash
for i in {0..29}
do
   ./build/ES grn5 lsoda de 1050000 $i >> de-output5-old-limits.txt &
done

exit

