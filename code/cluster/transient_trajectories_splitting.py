import scipy.stats as stats
import time
import datetime
from multiprocessing import Pool
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math

from multiprocessing import Pool
import os
from ast import literal_eval
import sys


def simulate(alpha, k):

	Xzero = 0.5


	dt = 0.005
	N = 10**3
	T = int(N/dt)
	
	epsilon = dt
	
	t = np.arange(dt, T+dt, dt)

	i = 0

	samplesize = 10
		
	traj_data = {}
	
	for s in range(0, samplesize):
		dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = s, size = T)
		dtdL = 3*dL*(dt)**(1/alpha)
		
		X = np.zeros(T) 
		
		Xtemp = Xzero
		
		for j in range(1, T):
			Linc = 0.1*dtdL[j]
			Xtemp = Xtemp - (dt*(Xtemp**3 - Xtemp + k))/(1 + epsilon*abs(Xtemp**3 - Xtemp + k)) + Linc
			X[j-1] = Xtemp
	
		traj_data[f'{s}'] = X# Convert the dictionary to a DataFrame
		
	# Convert the dictionary to a DataFrame		
	traj_df = pd.DataFrame(traj_data)
	
	# Save to CSV
	traj_df.to_csv(f"/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/trajectories_a{alpha}_k{round(k, 2)}.csv", index=False)
	
		
var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])

simulate(var1, var2)
				   
