
# coding: utf-8

# In[12]:


# %load 'horizontal_gene_transfer/Mesocentre/PDE_Full.py'
import matplotlib 
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime

from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import cm

def Pre_Initialization_PDE(parameters):
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= parameters['X_min'], parameters['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
       
    f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
    rho0=(parameters['b_r'])/(parameters['C'])
    f0=f0/(np.sum(f0)*parameters['dX'])*rho0

    f = np.empty([nT, nX])        #densities for all times
    f[0]=f0
    
    # Computing constant death and mutation kernels
    Death = parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
    Mutation_kernel=np.exp(-(X/(parameters['sigma']*parameters['eps']))**2*1/2)*parameters['dX']
    Mutation_kernel=parameters['b_r']*Mutation_kernel/np.sum(Mutation_kernel)
    
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

def Next_Generation_PDE(f,parameters, pre_init_values):
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= parameters['X_min'], parameters['X_max']
    Death, Mutation_kernel = pre_init_values['Death'], pre_init_values['Mutation_kernel']

    X_min_new=int(-X_min/dX)          # Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-2*X_min+X_max)/dX)
    I=range(pre_init_values['nX'])    # all the indexes
    
    rho = np.sum(f)*parameters['dX'] # I don't like the idea of comparing float with an integer, but okay
    if rho==0:
        return 0
    
    death_term = Death + rho*parameters['C']
    birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
    transfer_part = parameters['tau']/rho*parameters['dX']*np.fromiter((np.sum(f[:i])-np.sum(f[i:]) for i in I),float)
    #new_f = np.maximum(0,f + dT/parameters['eps']*(np.multiply(f,(-death_term+ transfer_part))+birth_part))
    new_f = f + dT/parameters['eps']*(np.multiply(f,(-death_term + transfer_part))+birth_part)
    return new_f

#######################################    
####### EXECUTABLE PART ###############
#######################################

parameters = dict(T_max = 100, # maximal time 
                  dT = 0.001, # Discretization time 
                  sigma0 = 0.01,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.01,
                  tau = 0.3,  # transfer rate
                  X_min = -0.2, #length of the numerical interval of traits (for PDE!)
                  X_max=1.5,
                  dX = 0.01, #discretization of the space of traits
                  eps = 1
                  )
    
# Loop over given parameters
param='tau' # here we check the value
param_m=0.2
param_M=0.55
param_step=25
J=np.arange(param_m,param_M,param_step)
for j in J:
    pre_init_values = Pre_Initialization_PDE(parameters) 
    f, nT = pre_init_values['f'], pre_init_values['nT']
    parameters[param]=j
    print(' !!!!!!!!!!!!!! '+str(param)+' = '+str(j), flush = True)    
    gc.collect()
    #SIMULATION !!!!!!!!!!!!!!!!!!
    for i in range(nT-1):
        f[i+1] = Next_Generation_PDE(f[i],parameters, pre_init_values)
        if i%1000==0:
            print('T= '+str(i*parameters['dT']), flush = True)

            

