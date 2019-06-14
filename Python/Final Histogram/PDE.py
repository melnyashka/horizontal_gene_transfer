
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime

def Pre_Initialization_PDE(parameters_PDE):
    dX, T_max, dT = parameters_PDE['dX'], parameters_PDE['T_max'], parameters_PDE['dT']
    X_min, X_max= parameters_PDE['X_min'], parameters_PDE['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
       
    f0 = np.exp(-np.power(np.absolute(X-parameters_PDE['x_mean0']),2)/parameters_PDE['sigma0']*parameters_PDE['eps'])# initial density 
    rho0=(parameters_PDE['b_r'])/(parameters_PDE['C'])
    f0=f0/(np.sum(f0)*parameters_PDE['dX'])*rho0

    f = np.empty([nT, nX])        #densities for all times
    f[0]=f0
    
    # Computing constant death and mutation kernels
    Death = parameters_PDE['d_r']*np.power(np.absolute(X),parameters_PDE['d_e'])
    Mutation_kernel=np.exp(-(X/(parameters_PDE['sigma']*parameters_PDE['eps']))**2*1/2)*parameters_PDE['dX']
    Mutation_kernel=parameters_PDE['b_r']*Mutation_kernel/np.sum(Mutation_kernel)
    
    pre_init_values = dict(
        f = f,
        T = T,
        nX = nX,
        nT = nT,
        X = X, 
        Death = Death, 
        Mutation_kernel = Mutation_kernel
        )
    return pre_init_values

def Next_Generation_PDE(f,parameters_PDE, pre_init_values):
    dX, T_max, dT = parameters_PDE['dX'], parameters_PDE['T_max'], parameters_PDE['dT']
    X_min, X_max= parameters_PDE['X_min'], parameters_PDE['X_max']
    Death, Mutation_kernel = pre_init_values['Death'], pre_init_values['Mutation_kernel']
    X_min_new=int(-X_min/dX)          # Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-2*X_min+X_max)/dX)
    
    rho = np.sum(f)*parameters_PDE['dX'] # I don't like the idea of comparing float with an integer, but okay
    
    death_term = Death + rho*parameters_PDE['C']
    birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
    if rho>np.power(10.,-7):

        X_large=np.arange(2*X_min,2*X_max,dX)
        transfer_kernel=  parameters_PDE['tau']*1/np.pi*(np.arctan(np.divide(X_large,parameters_PDE['delta']))-np.arctan(-np.divide(X_large,parameters_PDE['delta'])))*dX

        X_min_new=int(-2*X_min/dX)#Bounds for the new indexes that must be kept after the convolution
        X_max_new=int((-3*X_min+X_max)/dX)
        T = np.convolve(transfer_kernel,f/rho)[X_min_new:X_max_new]# that's the transfer term    
    else :
        T = 0
     #new_f = np.maximum(0,f + dT/parameters_PDE['eps']*(np.multiply(f,(-death_term+ transfer_part))+birth_part))
    new_f = f + dT/parameters_PDE['eps']*(np.multiply(f,(-death_term + T))+birth_part)
    return new_f



#######################################    
####### EXECUTABLE PART ###############
#######################################

if __name__ == "__main__":

    parameters_PDE = dict(T_max = 100, # maximal time 
                      dT = 0.01, # Discretization time 
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
                      tau = 0.4,  # transfer rate
                      X_min = -1, #length of the numerical interval of traits (for PDE!)
                      X_max=3,
                      dX = 0.1, #discretization of the space of traits
                      eps = 1,
                      delta=0.005
                      )
        
    # Loop over given parameters_PDE

pre_init_values = Pre_Initialization_PDE(parameters_PDE) 
U_P, nT = pre_init_values['f'], pre_init_values['nT']

for i in range(nT-1):
    U_P[i+1] = Next_Generation_PDE(U_P[i],parameters_PDE, pre_init_values)
    if i%1000==0:
        print('T= '+str(i*parameters_PDE['dT']), flush = True)
        plt.plot(U_P[i+1])
        plt.show()

