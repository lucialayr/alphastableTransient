import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from scipy.stats import levy_stable, moment
from scipy.special import gamma
from scipy import linalg
from scipy.integrate import quad
import pandas as pd
from ast import literal_eval
import sys
import psutil
import os

from splitting_functions import stable_rv, Phi

def get_memory_usage():
    """Get current memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # Convert to MB

def simulate(alpha, k, batch, num_simulations):

    # simulation setting
    # num_simulations is now passed as a parameter

    print(f"Starting simulation: alpha={alpha}, k={k}, batch={batch}, num_simulations={num_simulations}")

    Xzero = 1.0  # initial value

    if k == 0:
        T = 40 #bi-stable systems need longer to stabilize
    else:
         T = 10 # we integrate solutions on the time interval [0,T]
    N = 10000     # number of steps pro time unit
    dt = 1.0/N   # time mesh

    sigma = .6

    noise_amplitude = sigma*np.power(dt,1/alpha)

    X = [Xzero]*num_simulations # SDE values
    runs = [j + num_simulations*batch for j in range(num_simulations)] # Set run IDs once

    total_iterations = T*N
    # Print progress at 25%, 50%, 75%, and 100%
    print_milestones = [int(total_iterations * 0.25), int(total_iterations * 0.5), 
                        int(total_iterations * 0.75), total_iterations]
    
    for n in range(1, total_iterations + 1):
            for j in range(num_simulations):             
                xi = stable_rv(alpha, 1)[0]  # Extract scalar from array
                y =  X[j] - k*dt + noise_amplitude*xi
                X[j] = Phi(dt, y)
            
            # Print progress at 25% milestones
            if n in print_milestones:
                progress = (n / total_iterations) * 100
                print(f"[alpha={alpha}, k={k}, batch={batch}] Progress: {progress:.0f}% (iteration {n}/{total_iterations})")

    # Convert the results to a DataFrame
    results_df = pd.DataFrame({'run': runs, 'final_state': X})

    # Save to CSV
    results_df.to_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{batch}.csv", index=False)
    
    print(f"[alpha={alpha}, k={k}, batch={batch}] Simulation complete!")

var1 = literal_eval(sys.argv[1])  # alpha
var2 = literal_eval(sys.argv[2])  # k
var3 = literal_eval(sys.argv[3])  # batch
var4 = literal_eval(sys.argv[4])  # num_simulations

simulate(var1, var2, var3, var4)
