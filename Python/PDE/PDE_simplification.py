import matplotlib.pyplot as plt
import numpy as np
import gc 
def compute_things(parameters):
    L, dX = parameters['L'], parameters['dX']
    X_min, X_max= -L, L
    nX = int((X_max-X_min)/dX)
    # Create a new dictionary with: computed vector of traits, computed birth/death rate
    X = np.vectorize(lambda i: -parameters['L']+i*parameters['dX'],otypes=[float])(range(nX)) # reconstruct traits by index
    death_kernel = np.vectorize(lambda x: parameters['d_r']*np.power(np.absolute(x),parameters['d_e']))(X) 
    # ht_kernel = np.vectorize(lambda x: np.dot(parameters['tau'],np.heaviside(x-X, 1)))(X)
    
    init_density = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])
    rho0=np.sum(init_density)
    init_density=init_density/(parameters['C']*rho0)*(parameters['b_r'])
    new_dict = {'X' : X,
                'death_kernel': death_kernel, 
                'init_density': init_density,
                'X_min': X_min}
    return new_dict    

def tau_ht(x,y, n_tot, parameters): # horizontal transfer function
    tau, beta, mu = parameters['tau'], parameters['beta'], parameters['sigma']
    return tau*np.divide(np.heaviside(y-x, 1), (beta+n_tot*mu))

def Next_Generation_PDE(f,parameters, pre_values):
    dT, b_r, dX, sigma = parameters['dT'], parameters['b_r'], parameters['dX'], parameters['sigma']
    x = pre_values['X']
    n_tot = np.sum(f)
    death = pre_values['death_kernel'] + n_tot*parameters['C']
    birth_part = b_r/(sigma*parameters['eps'])*dX*np.vectorize(lambda y: np.sum(np.dot(f,np.exp(-(np.abs(y-x)/(parameters['eps']*sigma)**2)))), otypes=[float])(x)
    #sigma* laplacian of f
    #Laplacian_f= sigma*(np.append(f[1:],f[-1])+np.insert(f[:-1],0,f[-1])-2*f)/(parameters['dX']**2)
    transfer_part = np.vectorize(lambda x: parameters['tau']*(np.sum(f[:int((x-pre_values['X_min'])/parameters['dX'])])-np.sum(f[int((x-pre_values['X_min'])/parameters['dX']):]))/n_tot,otypes=[float])(x)
    new_f = np.maximum(0,f + dT/parameters['eps']*(-death*f + birth_part+ transfer_part*f))
    gc.collect()
    return new_f



