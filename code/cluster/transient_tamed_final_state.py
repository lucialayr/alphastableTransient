import scipy.stats as stats
import numpy as np
import pandas as pd
from ast import literal_eval
import sys

def simulate(alpha, k):
    Xzero = 0.5

    dt = 0.005
    N = 10**2
    T = int(N/dt)

    epsilon = dt
    samplesize = 20000

    # Initialize lists to store the results
    runs = []
    final_states = []

    for s in range(samplesize):
        dL = stats.levy_stable.rvs(alpha=alpha, beta=0, loc=0, scale=1, random_state=s, size=T)
        dtdL = 3*dL * (dt)**(1/alpha)

        X = np.zeros(T)
        Xtemp = Xzero

        for j in range(1, T):
            Linc = 0.1*dtdL[j]
            Xtemp = Xtemp - (dt*(Xtemp**3 - Xtemp + k))/(1 + epsilon*abs(Xtemp**3 - Xtemp + k)) + Linc
            X[j-1] = Xtemp

        # Store the run index and the final state
        runs.append(s)
        final_states.append(X[j-1])

    # Convert the results to a DataFrame
    results_df = pd.DataFrame({'run': runs, 'final_state': final_states})

    # Save to CSV
    results_df.to_csv(f"/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/final_states_a{alpha}_k{round(k, 2)}.csv", index=False)

var1 = literal_eval(sys.argv[1])
var2 = literal_eval(sys.argv[2])

simulate(var1, var2)
