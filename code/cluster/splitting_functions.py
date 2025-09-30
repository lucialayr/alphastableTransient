import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from scipy.stats import levy_stable, moment
from scipy.special import gamma
from scipy import linalg
from scipy.integrate import quad

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