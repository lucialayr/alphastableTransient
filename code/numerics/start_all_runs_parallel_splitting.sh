#!/bin/bash


for a in 2 1.5 1 0.5

do

for k in '-1' '-0.39' '0'

do

input="[${a},${k}]"

echo $input

python code/numerics/transient_trajectories_splitting.py $a $k &

python code/numerics/transient_final_state_splitting.py $a $k 0 &
python code/numerics/transient_final_state_splitting.py $a $k 1 &
python code/numerics/transient_final_state_splitting.py $a $k 2 &
python code/numerics/transient_final_state_splitting.py $a $k 3 &

done
done

wait

echo "All jobs finished."

for a in 2 1.5 1 0.5

do

for k in '-1' '-0.39' '0'

do

python code/numerics/merge_batches_data.py $a $k

rm "data/final_states_a${a}_k${k}_batch0.csv"
rm "data/final_states_a${a}_k${k}_batch1.csv"
rm "data/final_states_a${a}_k${k}_batch2.csv"
rm "data/final_states_a${a}_k${k}_batch3.csv"

done
done

echo "All data merged."



