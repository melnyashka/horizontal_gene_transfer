import numpy as np
import pandas as pd

T_max = 1500 # maximal time 
dT = 0.1 # Discretization time 

K = 1000     # Maximal capacity of the system
N0 = 1000    # Initial number of population
sigma0=0.1  #Initial standard variation of the population
C = 0.5    # competition
p = 0.3      # Probability of mutation
b_r = 1     # birth rate
d_r = 1      # death rate
beta = 0 
mu = 1
sigma = 0.01
tau = 0.17    # transfer rate



X0 = np.random.normal(1, sigma0, N0) # Initial population



X = [None]*int(T_max/dT)  # history of all populations up to time T_max
X[0] = np.sort(X0)

def horizontal_transfer(x):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return list(map(lambda i: ht_rate*(n_tot-i), range(n_tot)))
    


def Next_Generation(x):
    n_tot = x.size
    if n_tot==0:
        return x
    else:
        beta_birth = np.divide(1,np.repeat(b_r, n_tot))
        beta_death = np.divide(1,d_r*x**2 + n_tot*C/K)
        beta_transfer = np.divide(1,horizontal_transfer(x))
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(beta_transfer)])
        b_mat = (times < dT)

        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))


