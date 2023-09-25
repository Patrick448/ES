#!/bin/bash
for i in {0..30}
do
   ./build/ES grn5 lsoda cmaes 10500 $i >> output.txt &
done

exit
