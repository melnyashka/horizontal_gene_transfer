import matplotlib.pyplot as plt
import numpy as np
import gc 



#PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
parameters = dict(T_max = 500, # maximal time 
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
                  tau = 0.8,  # transfer rate
                  X_min = -2, #length of the numerical interval of traits (for PDE!)
                  X_max=4,
                  dX = 0.01, #discretization of the space of traits
                  eps = 1
                  )
dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
X_min, X_max= parameters['X_min'], parameters['X_max']


#GRID !!!!!!!!!!!!!!!!!!
#TIME
nT = int(T_max/dT)#number of times
T = np.fromiter((i*parameters_HJ['dT'] for i in range(nT)),float)#space of time

nX = int((X_max-X_min)/dX) # number of traits
X = np.fromiter((X_min+i*parameters_HJ['dX'] for i in range(nX)),float)#space of traits





#INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!               
f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
rho0=(parameters['b_r'])/(parameters['C']*100)
f0=f0/np.sum(f0)*rho0

f = np.empty([nT, nX])#densities for all times
f[0]=f0





#AUXILIARY FUNCTIONS

Death= parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
# ht_kernel = np.vectorize(lambda x: np.dot(parameters['tau'],np.heaviside(x-X, 1)))(X)



new_dict = {'X' : X,
            'death_kernel': death_kernel, 
            'init_density': init_density,
            'X_min': X_min}









for i in range(nT-1):
    XT[i+1] = Next_Generation_PDE(XT[i],parameters, pre_values)
    if i%10==0:
        print('T= '+str(i*parameters['dT']))










#PLOT !!!!!!!!!!!!!!!!!

import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import cm

figure = plt.figure()
im = imshow(XT,cmap=cm.coolwarm,aspect='auto')
colorbar(im)









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



