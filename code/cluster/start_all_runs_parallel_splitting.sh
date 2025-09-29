#!/bin/bash


for a in 2 1.5 1 0.5

do


for k in  '-1' '-0.39' '-0.2' 

do


input="[${a},${k},'test_data/']"

echo $input

python transient_trajectories_splitting.py $a $k 'test_data'&
python transient_final_state_splitting.py $a $k 'test_data'&

done
done

wait


