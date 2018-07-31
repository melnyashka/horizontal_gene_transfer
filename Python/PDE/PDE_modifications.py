import matplotlib.pyplot as plt
import numpy as np
import gc 

def i_to_x(i, parameters): # return the vector of traits by given vector of densities
    L, dX = parameters['L'], parameters['dX']
    return  -L+i*dX

def I_to_x(parameters):
    return np.vectorize(lambda i: i_to_x(i, parameters),otypes=[float])




def f_to_x(f, parameters): # return the vector of traits by given vector of densities
    L, dX = parameters['L'], parameters['dX']
    return list(map(lambda i: -L+i*dX, range(int(2*L/dX))))



def tau_ht(x,y, n_tot, parameters): # horizontal transfer function
    tau, beta, mu = parameters['tau'], parameters['beta'], parameters['sigma']
    return tau*np.divide(np.heaviside(y-x, 1), (beta+n_tot*mu))



def Next_Generation_PDE(f,parameters):
    dT, b_r, d_r, d_e, K, C  = parameters['dT'], parameters['b_r'], parameters['d_r'], parameters['d_e'], parameters['K'], parameters['C']
    L, dX, T_max, dT = parameters['L'], parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= -L, L
    nX = int((X_max-X_min)/dX) # number of traits
    x = np.arange(X_min,X_max,dX)
    n_tot = np.sum(f)
    death = d_r*np.power(np.absolute(x),d_e) + n_tot*C/K
    dX, sigma = parameters['dX'], parameters['sigma']
    birth_part = np.array(list(map(lambda y: b_r*np.sum(f*np.exp(-np.abs(x-y)/sigma))*dX, x))) # birth part of the integro-differential equation
    transfer_part = (f*np.array(list(map(lambda y: np.sum(np.divide(tau_ht(x,y,n_tot,parameters)*f, n_tot)), x))))*dX 
    new_f = f - dT*(death+n_tot)*f + birth_part + transfer_part
    gc.collect()
    return new_f



def birth_kernel(x,parameters):
    return parameters['b_r']*np.exp(-(x-X)**2/parameters['sigma'])
    
def birth_term(x,f,parameters):
    return np.dot(birth_kernel(x,parameters),f)

def Birth(f,parameters):
    return np.vectorize(lambda x: birth_term(x,f, parameters),otypes=[float])(X)

