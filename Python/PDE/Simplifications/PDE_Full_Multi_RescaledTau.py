import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime

import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import cm

#PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
parameters = dict(T_max = 0.3, # maximal time 
                  dT = 0.00001, # Discretization time 
                  sigma0 = 0.01,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 1,
                  tau = 1.,  # transfer rate
                  X_min = -0.2, #length of the numerical interval of traits (for PDE!)
                  X_max=1.5,
                  dX = 0.01, #discretization of the space of traits
                  eps = 0.001
                  )
dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
X_min, X_max= parameters['X_min'], parameters['X_max']


#GRID !!!!!!!!!!!!!!!!!!

#TIME
nT = int(T_max/dT)#number of times
T = np.arange(0,T_max,dT)#space of time

nX = int((X_max-X_min)/dX) # number of traits
I=range(nX)#all the indexes
X = np.arange(X_min,X_max,dX)#space of traits

#X_large=np.arange(2*X_min,2*X_max,dX)



#INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!               
f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
rho0=(parameters['b_r'])/(parameters['C']*500)
f0=f0/np.sum(f0)*rho0

f = np.empty([nT, nX])#densities for all times
f[0]=f0




#LOOP OVER PARAMETERS !!!!!!!!!!!!!!!!!!!
param='tau'
param_m=0
param_M=2
param_step=0.1
J=np.arange(param_m,param_M,param_step)
for j in J:
    f[0]=f0
    parameters[param]=j
    print(' !!!!!!!!!!!!!! '+str(param)+' = '+str(j))
    
    
    
    #AUXILIARY FUNCTIONS
    
    #DEATH
    Death= parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
    # ht_kernel = np.vectorize(lambda x: np.dot(parameters['tau'],np.heaviside(x-X, 1)))(X)
    
    #MUTATION
    sigma_eps= parameters['sigma']*parameters['eps']
    constant=np.sqrt(2*np.pi)
    Mutation_kernel=  parameters['b_r']/(sigma_eps*constant)*np.exp(-(X/sigma_eps)**2*1/2)*parameters['dX']
    X_min_new=int(-X_min/dX)#Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-2*X_min+X_max)/dX)
    
    
    def Next_Generation_PDE(f,parameters):
        rho = np.sum(f)
        if rho==0:
            return 0
        death_term = Death + rho*parameters['C']
        birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
        transfer_part = parameters['tau']/(rho*parameters['eps'])*parameters['dX']*np.fromiter((np.sum(f[:i])-np.sum(f[i:]) for i in I),float)
        new_f = np.maximum(0,f + dT/parameters['eps']*(np.multiply(f,(-death_term+ transfer_part))+birth_part))
        return new_f
    
    
    #SIMULATION !!!!!!!!!!!!!!!!!!


    for i in range(nT-1):
        f[i+1] = Next_Generation_PDE(f[i],parameters)
        if i%1000==0:
            print('T= '+str(i*parameters['dT']))












#PLOT !!!!!!!!!!!!!!!!!


    

    figure = plt.figure()
    im = imshow(f[100:].transpose()[::-1],cmap=cm.coolwarm,aspect='auto',extent=(0,parameters['T_max'],X_min,X_max),vmin=0)
    colorbar(im)
    
    plt.xlabel('time')
    plt.ylabel('trait')
    
    par_str = '' # create the title, with values of parameters
    for k, v in parameters.items():
        if k == 'N0' or k == 'd_r' or k=='X_min': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    par_str=par_str+' rho0='+str(rho0)
    plt.title(par_str)
    
    plt.tight_layout()
    #levellines=np.array([-0.01])
    #cset = contour(u.transpose(),levellines,linewidths=2,cmap=cm.Set2)
    #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
    
    current_time = datetime.now().time()
    figure.savefig(str("Figures/Full/RescaledTau_plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!
    



    gc.collect()








