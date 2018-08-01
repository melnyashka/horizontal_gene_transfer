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
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.1,
                  tau = 0.,  # transfer rate
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
T = np.fromiter((i*parameters['dT'] for i in range(nT)),float)#space of time

nX = int((X_max-X_min)/dX) # number of traits
X = np.fromiter((X_min+i*parameters['dX'] for i in range(nX)),float)#space of traits





#INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!               
f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
rho0=(parameters['b_r'])/(parameters['C'])
f0=f0/np.sum(f0)*rho0

f = np.empty([nT, nX])#densities for all times
f[0]=f0





#AUXILIARY FUNCTIONS

Death= parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
# ht_kernel = np.vectorize(lambda x: np.dot(parameters['tau'],np.heaviside(x-X, 1)))(X)


def Next_Generation_PDE(f,parameters):
    dT, b_r, dX, sigma = parameters['dT'], parameters['b_r'], parameters['dX'], parameters['sigma']
    rho = np.sum(f)
    death_term = Death + rho*parameters['C']
    birth_part = b_r/(sigma*parameters['eps'])*dX*np.vectorize(lambda y: np.sum(np.dot(f,np.exp(-(np.abs(y-X)/(parameters['eps']*sigma)**2)))), otypes=[float])(X)
    Laplacian_f= sigma*(np.append(f[1:],f[-1])+np.insert(f[:-1],0,f[-1])-2*f)/(parameters['dX']**2)
    transfer_part = np.vectorize(lambda x: parameters['tau']*(np.sum(f[:int((x-X_min)/dX)])-np.sum(f[int((x-X_min)/dX):]))/rho,otypes=[float])(X)
    #transfer_part=0
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








