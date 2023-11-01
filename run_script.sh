#!/bin/bash
for i in {0..29}
do
   ./build/ES grn5 lsoda cmaes 105000 $i >> output.txt &
done

exit
