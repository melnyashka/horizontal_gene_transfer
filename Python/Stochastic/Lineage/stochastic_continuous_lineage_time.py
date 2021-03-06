import numpy as np
import matplotlib.pyplot as plt
import gc
from datetime import datetime
from itertools import compress

def horizontal_transfer(x, tau, beta, mu):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return list(map(lambda i: ht_rate*(n_tot-i), range(n_tot)))
    


def Next_Generation_lineage(x,l,t, parameters):
    b_r, d_r, C, K, sigma, d_e = parameters['b_r'], parameters['d_r'], parameters['C'], parameters['K'], parameters['sigma'],parameters['d_e']
    n_tot = x.size
    f=np.vectorize(lambda e: np.random.choice(x[x>=e]),otypes=[np.float64])
    if n_tot==0:
        return x
    else:
        beta_birth = np.divide(1,np.repeat(b_r, n_tot))
        beta_death = np.divide(1,d_r*np.power(np.absolute(x),d_e) + n_tot*C/K)
        beta_transfer = np.divide(1,horizontal_transfer(x, tau = parameters['tau'], beta = parameters['beta'], mu = parameters['mu']))
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(beta_transfer)])
        b_mat = (times < parameters['dT'])
        x_birth=x[b_mat[0]]
        l_add_birth=l[b_mat[0]]
        for i in range(x_birth.size):
            l_add_birth[i]=l_add_birth[i]+[(t,x_birth[i])]
        x_HT=x[b_mat[2]]
        if x_HT.size==0:
            x_HT_add=x_HT
        else:
            x_HT_add=f(x_HT)
        l_add_HT=l[b_mat[2]]
        for i in range(x_HT.size):
            l_add_HT[i]=l_add_HT[i]+[(t,x_HT[i])]
        x_new=np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       np.random.normal(loc=x_birth, scale=sigma, size=None),
                                       x_HT_add))
        l_new= np.concatenate((l[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       l_add_birth,
                                       l_add_HT))
        return x_new,l_new

