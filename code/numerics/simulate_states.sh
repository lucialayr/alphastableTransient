#!/bin/bash

# Default number of simulations per batch (can be overridden by command line argument)
NUM_SIMS=${1:-10000}

echo "========================================"
echo "Starting parallel simulations"
echo "Number of simulations per batch: ${NUM_SIMS}"
echo "Number of batches per parameter combination: 10"
echo "Total simulations per parameter combination: $((NUM_SIMS * 10))"
echo "Time: $(date)"
echo "========================================"

for a in 2 1.5 1 0.5
do
    for k in '-1' '-0.39' '0'
    do
        echo ""
        echo "Launching parameter combination: alpha=${a}, k=${k}"
        echo "Starting 10 parallel batches (${NUM_SIMS} simulations each) at $(date)"
        
        # Launch 10 batches in parallel for this parameter combination
        python code/numerics/transient_final_state_splitting.py $a $k 0 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 1 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 2 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 3 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 4 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 5 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 6 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 7 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 8 ${NUM_SIMS} &
        python code/numerics/transient_final_state_splitting.py $a $k 9 ${NUM_SIMS} &
        
        # Wait for all 10 batches to complete before starting next parameter combination
        wait
        
        echo "Completed parameter combination: alpha=${a}, k=${k} at $(date)"
    done
done

echo ""
echo "========================================"
echo "All simulation jobs finished at $(date)"
echo "========================================"

echo ""
echo "========================================"
echo "Starting data merging at $(date)"
echo "========================================"

for a in 2 1.5 1 0.5
do
    for k in '-1' '-0.39' '0'
    do
        echo "Merging data for alpha=${a}, k=${k}"
        
        # Verify all batch files exist before merging
        all_files_exist=true
        for batch in {0..9}; do
            if [ ! -f "data/final_states_a${a}_k${k}_batch${batch}.csv" ]; then
                all_files_exist=false
                echo "WARNING: Missing batch file ${batch} for alpha=${a}, k=${k}"
                break
            fi
        done
        
        if [ "$all_files_exist" = true ]; then
            python code/numerics/merge_batches_data.py $a $k
            
            # Clean up batch files after successful merge
            for batch in {0..9}; do
                rm "data/final_states_a${a}_k${k}_batch${batch}.csv"
            done
            
            echo "Successfully merged and cleaned up batches for alpha=${a}, k=${k}"
        else
            echo "WARNING: Missing batch files for alpha=${a}, k=${k}. Skipping merge."
        fi
    done
done

echo ""
echo "========================================"
echo "All data merged at $(date)"
echo "========================================"



