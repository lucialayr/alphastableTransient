#!/bin/bash


for a in 2 1.5 1 0.5

do


for k in  '0' #'-1' '-0.39'  

do


input="[${a},${k}]"

echo $input

python code/cluster/transient_trajectories_splitting.py $a $k &
python code/cluster/transient_final_state_splitting.py $a $k &

done
done

wait

echo "All jobs finished."



