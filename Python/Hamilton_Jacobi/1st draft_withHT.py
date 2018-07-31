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


#Parameters
parameters_HJ = dict(T_max = 10, # maximal time 
                dT = 0.1, # Discretization time 
                C=0.5, # carrying capacity of the system
                rho0 = 1000,    # Initial number of population
                sigma0=1,  #Initial standard variation of the population
                x0=0.,
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                d_e=2,      #exponent for the death function
                beta = 0, 
                mu = 1,
                sigma = 0.01,
                tau = 2, # transfer rate
                L=4, #space trait X=[-L,L]
                dX=1/100, #discretization of space trait
                u_inf=-50
                )



#Grid of time and trait
nT=int(parameters_HJ['T_max']/parameters_HJ['dT']) #number of times
T=[t*parameters_HJ['dT'] for t in range(nT)] #list of all times
T= np.array(T)

X_min = -parameters_HJ['L']   # Minimal trait
X_max = parameters_HJ['L']  # Maximum amount of traits
nX =int((X_max-X_min)/parameters_HJ['dX']) # number of traits
X = [X_min+x*parameters_HJ['dX'] for x in range(nX)] # list of possible traits
X=np.array(X)


def x_to_i_HJ (x):
    return int(x-X_min)/parameters_HJ['dX']
def i_to_x_HJ (i):
    return X_min+i*parameters_HJ['dX']





#INITIAL TIME
u0=np.fromfunction(lambda j,i:  -(i_to_x_HJ(i)-parameters_HJ['x0'])**2/parameters_HJ['sigma0']**2, (1,nX),dtype=float)
u= np.empty((nT,nX),dtype=float)
u[0]=u0[0]



#DEATH:
def death (x):
    return parameters_HJ['d_r']*x**parameters_HJ['d_e']
Death= (np.vectorize(death,otypes=[np.float64]))(X)




#BIRTH:
def birth_kernel(p,z):
    return np.exp(-(z/parameters_HJ['sigma'])**2+p*z)*parameters_HJ['b_r']
def birth_term(p):
    return np.sum(np.vectorize(lambda x: birth_kernel(p,x),otypes=[np.float64])(X))*parameters_HJ['dX']

Birth= np.vectorize(birth_term,otypes=[np.float64])



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
    u_add=-Death+Birth(grad_u)+HT_term(x_M)
    u_new= u+u_add*parameters_HJ['dT']
    
    return np.maximum(u_new-np.max(u_new), parameters_HJ['u_inf'])
    





#Simulation
for i in range(int(parameters_HJ['T_max']/parameters_HJ['dT']-1)):
    if i%10==0:
        print('T= '+str(i*parameters_HJ['dT']))
    u[i+1]=Next_time_HJ(u[i], parameters_HJ)
    
    
    
    
    
    
#Plot

from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

figure = plt.figure()
im = imshow(u,cmap=cm.coolwarm)
colorbar(im)
par_str = '' # create a string of parameters to pass into plots
for k, v in parameters_HJ.items():
    if k == 'N0' or k == 'b_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth

plt.ylabel('time')
plt.xlabel('trait');
plt.title(par_str)
levellines=np.array([-10,-2,-0.01])
cset = contour(u,levellines,linewidths=2,cmap=cm.Set2)
clabel(cset,inline=True,fmt='%1.1f',fontsize=10)


current_time = datetime.now().time()
figure.savefig(str("Figures/plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!







#plt.hist2d(u,bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)

    

    