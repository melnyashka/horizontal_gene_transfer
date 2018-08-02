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
parameters_HJ = dict(T_max = 5, # maximal time 
                dT = 1/100, # Discretization time 
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
                tau = 0., # transfer rate
                X_min=-3,
                X_max=3,
                dX=1/100, #discretization of space trait
                u_inf=-5000
                )



#Grid of time and trait
nT=int(parameters_HJ['T_max']/parameters_HJ['dT']) #number of times
T = np.fromiter((i*parameters_HJ['dT'] for i in range(nT)),float)

dX,X_min,X_max= parameters_HJ['dX'],parameters_HJ['X_min'],parameters_HJ['X_max']
nX =int((X_max-X_min)/parameters_HJ['dX']) # number of traits
X = np.fromiter((X_min+i*parameters_HJ['dX'] for i in range(nX)),float)#space of traits

Iplus= np.arange(int(-X_min/dX),nX,1)#indexes for positive values of x
Iminus=np.arange(0,int(-X_min/dX),1)#indexes for negative values of x
Xplus= X[Iplus]
Xminus=X[Iminus]

def x_to_i_HJ (x):
    return int(x-X_min)/parameters_HJ['dX']
def i_to_x_HJ (i):
    return X_min+i*parameters_HJ['dX']





#INITIAL TIME
u0=np.fromiter((-1/2*((x-parameters_HJ['x0'])/parameters_HJ['sigma0'])**2 for x in X),dtype=float)
u= np.empty((nT,nX),dtype=float)
u[0]=u0



#DEATH:

Death= parameters_HJ['d_r']*X**parameters_HJ['d_e']




#BIRTH:
    
Birth=parameters_HJ['b_r']    
    
#def birth_kernel(p,z):
#    return np.exp(-(z/parameters_HJ['sigma'])**2+p*z)*parameters_HJ['b_r']
#def birth_term(p):
#    return np.sum(np.vectorize(lambda x: birth_kernel(p,x),otypes=[np.float64])(X))*parameters_HJ['dX']
#
#Birth= np.vectorize(birth_term,otypes=[np.float64])
#

#def birth_plus (p):#p is the gradient for positive z 
#    return np.sum(parameters_HJ['b_r']*parameters['dX']*np.exp(np.multiply(Xplus,p-Xplus/(2*parameters['sigma']**2))))
#Birth_plus=np.vectorize(birth_plus, otypes=[np.float64])
#
#def birth_minus (p):#and here p is for negative z
#    return np.sum(parameters_HJ['b_r']*parameters['dX']*np.exp(np.multiply(Xminus,p-Xminus/(2*parameters['sigma']**2))))
#Birth_minus=np.vectorize(birth_plus, otypes=[np.float64])



#HORIZONTAL TRANSFER:
def HT(y):
    return (np.arctan(X-y)-np.arctan(y-X))*parameters_HJ['tau']


    




#Evolution Operator
def Next_time_HJ(u,parameters_HJ):
    grad_u= (u[1:]-u[:-1])/parameters_HJ['dX']
    grad2_u=np.power(grad_u,2)
    grad2_u=(np.insert(grad2_u,0,2*grad2_u[0]-grad2_u[1])+np.append(grad2_u,2*grad2_u[-1]-grad2_u[-2]))/2
    x_M=i_to_x_HJ(np.argmax(u))
    u_add=Birth-Death+(parameters_HJ['sigma']**2/2)+HT(x_M)
    u_new= u+u_add*parameters_HJ['dT']
    M=np.max(u_new)
    return u_new-max(0,M),M
    




#Simulation
M=[] #maximum of u
m=0.
t_extinction=[]
for i in range(int(parameters_HJ['T_max']/parameters_HJ['dT']-1)):
    if i%500==0:
        print('T= '+str(i*parameters_HJ['dT']))
        print('m = '+str(m))
    u[i+1],m=Next_time_HJ(u[i], parameters_HJ)
    if m<0:
        M=M+[m]
        t_extinction=t_extinction+[i*parameters_HJ['dT']]

    
    
    
    
    
    
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
    if k == 'N0' or k == 'd_r': 
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
figure.savefig(str("Figures/Simplified/plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!







#plt.hist2d(u,bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)

    

    