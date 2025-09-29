#!/bin/bash


for a in 2 1.5 1 0.5

do


for k in  '-1' '-0.39' '-0.2' 

do


input="[${a},${k}]"

echo $input

srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python transient_trajectories_splitting.py $a $k &
srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python transient_final_state_splitting.py $a $k &

done
done

wait


