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

def simulate(alpha, k):

	Xzero = 1.0  # initial value
	T = 50       # we integrate solutions on the time interval [0,T]
	N = 1000     # number of steps pro time unit
	dt = 1.0/N   # time mesh

	sigma = .6

	noise_amplitude = sigma*np.power(dt,1/alpha)

	num_simulations = 1

	# Store all raw simulation paths
	sim_paths = np.zeros((num_simulations, T * N + 1))

	for sim_idx in range(num_simulations):
		X = np.zeros(T * N + 1)
		X[0] = Xzero
		for n in range(1, T * N + 1):
			xi = stable_rv(alpha)
			y = X[n - 1] - k*dt + noise_amplitude*xi
			X[n] = Phi(dt, y)
		sim_paths[sim_idx, :] = X
	
	# Convert the dictionary to a DataFrame		
	traj_df = pd.DataFrame(sim_paths)
	
	# Save to CSV
	traj_df.to_csv(f"test_data/trajectories_a{alpha}_k{round(k, 2)}.csv", index=False)
	
var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])

simulate(var1, var2)

				   
