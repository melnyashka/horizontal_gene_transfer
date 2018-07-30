import matplotlib.pyplot as plt
import numpy as np

parameters = dict(T_max = 1000, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                  sigma0 = 0.1,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.01,
                  tau = 0.17,  # transfer rate
                  L = 10, #length of the numerical interval of traits (for PDE!)
                  dX = 0.001, #discretization of the space of traits
                  eps = 1
                  )
X_min, X_max= -L, L
nX = int((X_max-X_min)/dX) # number of traits
X = np.arange(0,dX*nX,dX) # list of possible traits
F0 = np.exp(-np.power(np.absolute(X),2)/parameters['sigma0']*parameters['eps'])) # initial density 

def x_to_i(x, parameters): # return the index by the trait
    return int(np.divide((x + L), dX))

def f_to_x(f, parameters): # return the vector of traits by given vector of densities
    L, dX = parameters['L'], parameters['dX']
    return np.vectorize(map(lambda i: -L+i*dX, range(len(f))))

def tau_ht(x,y, n_tot, parameters): # horizontal transfer function
    tau, beta, mu = parameters['tau'], parameters['beta'], parameters['sigma']
    return tau*np.divide(np.heaviside(y-x, 1), (beta+n_tot*mu))

def Next_Generation_PDE(f,parameters)
    dT, d_r, d_e, K, C  = parameters['dT'], parameters['d_r'], parameters['d_e'], parameters['K'], parameters['C']
    n_tot = np.sum(f)
    x = f_to_x(f, parameters)
    death = d_r*np.power(np.absolute(x),d_e) + n_tot*C/K

    dX, eps = parameters['dX'], parameters['eps']
    birth_part = np.vectorize(map(lambda y: b_r*np.sum(f*np.exp(-np.abs(x-y)/eps))*dX, x)) # birth part of the integro-differential equation
    transfer_part = f*np.vectorize(map(lambda y: np.divide(tau_ht(x,y,n_tot,parameters)*f, np.sum(f)), x))*dX 
    new_f = f - dT*(death+n_tot)*f + birth_part + transfer_part
    return new_f
    


