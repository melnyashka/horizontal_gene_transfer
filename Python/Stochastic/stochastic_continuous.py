import numpy as np
import pandas as pd

T_max = 10 # maximal time 
dT = 0.1 # Discretization time 

K = 1000     # Maximal capacity of the system
N0 = 1000    # Initial number of population
C = 0.5/K    # competition
p = 0.3      # Probability of mutation
b_r = 1     # birth rate
d_r = 1      # death rate
beta = 0 
mu = 1
tau = 0.6    # transfer rate
sigma = 0.1

X0 = np.random.normal(1, 0.1, N0) # Initial population



X = [None]*int(T_max/dT)  # history of all populations up to time T_max
X[0] = np.sort(X0)

def horizontal_transfer(x):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return list(map(lambda i: ht_rate*(n_tot-i), range(n_tot)))
    


def Next_Generation(x):
    if x.size==0:
        return x
    else:
        n_tot = len(x)
        print(n_tot)
        lambda_birth = b_r
        lambda_death = d_r*x**2 + n_tot*C
        lambda_transfer = horizontal_transfer(x)
        lambda_total =  lambda_birth + lambda_death + lambda_transfer
        times = np.array([np.random.exponential(lambda_total/lambda_birth),np.random.exponential(lambda_total/lambda_death), np.random.exponential(lambda_total/lambda_transfer)])
        b_mat = (times < dT)

        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))


