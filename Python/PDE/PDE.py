import matplotlib.pyplot as plt
import numpy as np
import gc 

def x_to_i(x, parameters): # return the index by the trait
    return int(np.divide((x + L), dX))

def f_to_x(f, parameters): # return the vector of traits by given vector of densities
    L, dX = parameters['L'], parameters['dX']
    return list(map(lambda i: -L+i*dX, range(len(f))))

def tau_ht(x,y, n_tot, parameters): # horizontal transfer function
    tau, beta, mu = parameters['tau'], parameters['beta'], parameters['sigma']
    return tau*np.divide(np.heaviside(y-x, 1), (beta+n_tot*mu))

def Next_Generation_PDE(f,parameters):
    dT, b_r, d_r, d_e, K, C  = parameters['dT'], parameters['b_r'], parameters['d_r'], parameters['d_e'], parameters['K'], parameters['C']
    n_tot = np.sum(f)
    x = np.array(f_to_x(f, parameters))
    death = d_r*np.power(np.absolute(x),d_e) + n_tot*C/K
    dX, eps = parameters['dX'], parameters['eps']
    birth_part = np.array(list(map(lambda y: b_r*np.sum(f*np.exp(-np.abs(x-y)/eps))*dX, x))) # birth part of the integro-differential equation
    transfer_part = (f*np.array(list(map(lambda y: np.sum(np.divide(tau_ht(x,y,n_tot,parameters)*f, n_tot)), x))))*dX 
    new_f = f - dT*(death+n_tot)*f + birth_part + transfer_part
    gc.collect()
    return new_f

