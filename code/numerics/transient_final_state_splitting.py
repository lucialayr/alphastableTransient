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

from splitting_functions import stable_rv, Phi

def simulate(alpha, k, batch):

    # simulation setting
    num_simulations = 25000  # number of samples in Monte Carlo

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
    runs = [Xzero]*num_simulations # SDE values

    for n in range(1,T*N+1):
            for j in range(num_simulations):             
                xi = stable_rv(alpha)
                y =  X[j] - k*dt + noise_amplitude*xi
                X[j] = Phi(dt, y)
                runs[j] = j + 25000*batch #make sure that run ID is unique across batches

    # Convert the results to a DataFrame
    results_df = pd.DataFrame({'run': runs, 'final_state': X})

    # Save to CSV
    results_df.to_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{batch}.csv", index=False)

var1 = literal_eval(sys.argv[1])
var2 = literal_eval(sys.argv[2])
var3 = literal_eval(sys.argv[3])

simulate(var1, var2, var3)
