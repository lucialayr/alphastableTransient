import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from scipy.stats import levy_stable, moment
from scipy.special import gamma
from scipy import linalg
from scipy.integrate import quad

def Phi(t,x):
    ''' Approximation Phi for Splitting Method'''
    return x/np.sqrt( x*x*(1-np.exp(-2*t))  +  np.exp(-2*t)  )

def stable_rv(alpha, size):
    """
    Generates an array of alpha-stable random variables.
    The seed must be set immediately before calling this function.
    """
    if alpha == 2:
        return rng.standard_normal(size=size)
    elif alpha == 1:
        return rng.standard_cauchy(size=size)
    else:
        U = np.pi * (rng.random(size=size) - 0.5)  # uniform r.v. on [-Pi/2,Pi/2]
        E = -np.log(rng.random(size=size))         # standard exponential r.v.
        a = np.sin(alpha * U) / np.power(np.cos(U), 1 / alpha)
        a = a * np.power(np.cos(U * (1 - alpha)) / E, (1 - alpha) / alpha)
        return a

    
