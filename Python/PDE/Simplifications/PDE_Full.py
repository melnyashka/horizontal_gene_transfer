import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime


#PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
parameters = dict(T_max = 100, # maximal time 
                  dT = 0.01, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                  sigma0 = 0.1,  #Initial standard variation of the population
                  x_mean0 = 1.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.1,
                  tau = 0.2,  # transfer rate
                  X_min = -1, #length of the numerical interval of traits (for PDE!)
                  X_max=2.5,
                  dX = 0.01, #discretization of the space of traits
                  eps = 1
                  )
dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
X_min, X_max= parameters['X_min'], parameters['X_max']


#GRID !!!!!!!!!!!!!!!!!!

#TIME
nT = int(T_max/dT)#number of times
T = np.arange(0,T_max,dT)#space of time

nX = int((X_max-X_min)/dX) # number of traits
I=range(nX)
X = np.arange(X_min,X_max,dX)#space of traits

#X_large=np.arange(2*X_min,2*X_max,dX)



#INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!               
f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
rho0=(parameters['b_r'])/(parameters['C'])
f0=f0/np.sum(f0)*rho0

f = np.empty([nT, nX])#densities for all times
f[0]=f0





#AUXILIARY FUNCTIONS

#DEATH
Death= parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
# ht_kernel = np.vectorize(lambda x: np.dot(parameters['tau'],np.heaviside(x-X, 1)))(X)

#MUTATION
sigma_eps= parameters['sigma']*parameters['eps']
constant=np.sqrt(2*np.pi)
Mutation_kernel=  parameters['b_r']/(sigma_eps*constant)*np.exp(-(X/sigma_eps)**2)*parameters['dX']
X_min_new=int(-X_min/dX)
X_max_new=int((-2*X_min+X_max)/dX)


def Next_Generation_PDE(f,parameters):
    dT, b_r, dX, sigma = parameters['dT'], parameters['b_r'], parameters['dX'], parameters['sigma']
    rho = np.sum(f)
    death_term = Death + rho*parameters['C']
    birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
    transfer_part = parameters['tau']/rho*np.fromiter((np.sum(f[:i])-np.sum(f[i:]) for i in I),float)
    new_f = np.maximum(0,f + dT/parameters['eps']*(f*(-death_term+ transfer_part)+birth_part))
    return new_f


#SIMULATION !!!!!!!!!!!!!!!!!!
for i in range(nT-1):
    f[i+1] = Next_Generation_PDE(f[i],parameters)
    if i%100==0:
        print('T= '+str(i*parameters['dT']))












#PLOT !!!!!!!!!!!!!!!!!

import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import cm

figure = plt.figure()
im = imshow(f.transpose()[::-1],cmap=cm.coolwarm,aspect='auto',extent=(0,parameters['T_max'],X_min,X_max),vmin=0)
colorbar(im)

plt.xlabel('time')
plt.ylabel('trait')

par_str = '' # create the title, with values of parameters
for k, v in parameters.items():
    if k == 'N0' or k == 'd_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth
plt.title(par_str)

plt.tight_layout()
#levellines=np.array([-0.01])
#cset = contour(u.transpose(),levellines,linewidths=2,cmap=cm.Set2)
#clabel(cset,inline=True,fmt='%1.1f',fontsize=10)

current_time = datetime.now().time()
figure.savefig(str("Figures/Full/plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!



gc.collect()








