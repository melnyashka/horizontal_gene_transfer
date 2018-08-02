#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 17:00:05 2018

@author: samuelnordmann
"""
import numpy as np
import matplotlib.pyplot as plt
import gc
from datetime import datetime
from itertools import compress
#matplotlib.use("Agg")

#Parameters
parameters_HJ = dict(T_max = 3, # maximal time 
                dT = 1/10000, # Discretization time 
                C=0.5, # carrying capacity of the system
                rho0 = 1000,    # Initial number of population
                sigma0=1,  #Initial standard variation of the population
                x0=0.,
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                d_e=2,      #exponent for the death function
                beta = 0, 
                mu = 1,
                sigma = 1,
                tau = 0, # transfer rate
                X_min=-10,
                X_max=30,
                dX=1/100, #discretization of space trait
                u_inf=-50
                )



#Grid of time and trait
nT=int(parameters_HJ['T_max']/parameters_HJ['dT']) #number of times
T = np.fromiter((i*parameters_HJ['dT'] for i in range(nT)),float)

X_min = parameters_HJ['X_min']   # Minimal trait
X_max = parameters_HJ['X_max']  # Maximum amount of traits
nX =int((X_max-X_min)/parameters_HJ['dX']) # number of traits
X = np.fromiter((X_min+i*parameters_HJ['dX'] for i in range(nX)),float)#space of traits


def x_to_i_HJ (x):
    return int(x-X_min)/parameters_HJ['dX']
def i_to_x_HJ (i):
    return X_min+i*parameters_HJ['dX']





#INITIAL TIME
u0=np.fromfunction(lambda j,i:  -(i_to_x_HJ(i)-parameters_HJ['x0'])**2/parameters_HJ['sigma0']**2, (1,nX),dtype=float)
u= np.empty((nT,nX),dtype=float)





#!!!!!!!!!!!!!!!! LOOP !!!!!!!!!!!!!!!!!!!!!


param='tau'
param_m=0
param_M=101
param_step=10
J=np.arange(param_m,param_M,param_step)
for j in J:
    u[0]=u0[0]
    parameters_HJ[param]=j
    print(' !!!!!!!!!!!!!! '+str(param)+' = '+str(j))
        
    #DEATH:
    def death (x):
        return parameters_HJ['d_r']*x**parameters_HJ['d_e']
    Death= (np.vectorize(death,otypes=[np.float64]))(X)
    
    
    
    
    #BIRTH:
    Birth=parameters_HJ['b_r']
    #def birth_kernel(p,z):
    #    return np.exp(-(z/parameters_HJ['sigma'])**2+p*z)*parameters_HJ['b_r']
    #def birth_term(p):
    #    return np.sum(np.vectorize(lambda x: birth_kernel(p,x),otypes=[np.float64])(X))*parameters_HJ['dX']
    #
    #Birth= np.vectorize(birth_term,otypes=[np.float64])
    #
    
    
    #HORIZONTAL TRANSFER:
    def tau_kernel(x,y):
        return (np.arctan(x-y)-np.arctan(y-x))*parameters_HJ['tau']
    def HT_term(y):
        return np.vectorize(lambda x: tau_kernel(x,y), otypes=[np.float64])(X)
    
    
        
    
    
    
    
    #Evolution Operator
    def Next_time_HJ(u,parameters_HJ):
        grad_u= (u[1:]-u[:-1])/parameters_HJ['dX']
        grad_u=np.insert(grad_u,0,grad_u[0])
        x_M=i_to_x_HJ(np.argmax(u))
        u_add=Birth-Death+parameters_HJ['sigma']**2/2*np.dot(grad_u,grad_u)+HT_term(x_M)
        u_new= u+u_add*parameters_HJ['dT']
        M=np.max(u_new)
        return np.maximum(u_new-max(0,M), parameters_HJ['u_inf']),M
        
    
    
    
    
    
    #Simulation
    M=[] #maximum of u
    m=0.
    t_extinction=[]
    for i in range(int(parameters_HJ['T_max']/parameters_HJ['dT']-1)):
        if i%100==0:
            print('T= '+str(i*parameters_HJ['dT']))
            #print('m = '+str(m))
        u[i+1],m=Next_time_HJ(u[i], parameters_HJ)
        if m<0:
            M=M+[m]
            t_extinction=t_extinction+[i*parameters_HJ['dT']]
            print('coucou!') 
    
    
        
        
        
        
        
        
    #Plot
    
    from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    
    
    figure = plt.figure()
    #im = imshow(u.transpose(),cmap=cm.coolwarm)
    im = imshow(u.transpose()[::-1],cmap=cm.coolwarm,aspect='auto',extent=(0,parameters_HJ['T_max'],X_min,X_max),vmin=-30)
    colorbar(im)
    par_str = '' # create a string of parameters to pass into plots
    for k, v in parameters_HJ.items():
        if k == 'N0' or k == 'd_r' or k=='tau': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    
    plt.xlabel('time')
    plt.ylabel('trait');
    
    string_e='\nNumber of extinctions: '+str(len(t_extinction))
    if len(t_extinction)>0:
        string_e=string_e+'. First extinction at time t= '+str(t_extinction[0])
    plt.title(par_str+string_e)
    plt.tight_layout()
    #levellines=np.array([-0.01])
    #cset = contour(u.transpose(),levellines,linewidths=2,cmap=cm.Set2)
    #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
    
    
    current_time = datetime.now().time()
    figure.savefig(str("Figures/Simplified_PositiveCorrection/plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!
#   figure.savefig(str("scratch/gene/Figures" +"plot_" + str(current_time)[0:8]+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!






#plt.hist2d(u,bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)

    

    