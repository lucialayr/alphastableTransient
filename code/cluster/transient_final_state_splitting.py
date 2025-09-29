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

def Phi(t,x):
    return x/np.sqrt( x*x*(1-np.exp(-2*t))  +  np.exp(-2*t)  )

def stable_rv(alpha):
    if alpha == 2:
        return rng.standard_normal()
    elif alpha == 1: 
        return rng.standard_cauchy()
    else:
        U = np.pi*(rng.random()-0.5)  # uniform r.v. on [-Pi/2,Pi/2]
        E = -np.log(rng.random())     # standard exponential r.v.
        a=np.sin(alpha*U)/np.power(np.cos(U),1/alpha)
        a=a*np.power(np.cos(U*(1-alpha))/E,(1-alpha)/alpha)
        return(a)

def simulate(alpha, k, path):

    # simulation setting
    num_simulations = 10000  # number of samples in Monte Carlo

    Xzero = 1.0  # initial value
    T = 20       # we integrate solutions on the time interval [0,T]
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
                run[j] = j

    # Convert the results to a DataFrame
    results_df = pd.DataFrame({'run': runs, 'final_state': X})

    # Save to CSV
    results_df.to_csv(f"{path}final_states_a{alpha}_k{round(k, 2)}.csv", index=False)

var1 = literal_eval(sys.argv[1])
var2 = literal_eval(sys.argv[2])
var3 = literal_eval(sys.argv[3])
simulate(var1, var2, var3)
