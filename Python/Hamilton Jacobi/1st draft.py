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




parameters_HJ = dict(T_max = 1000, # maximal time 
                dT = 0.1, # Discretization time 
                C=0.5, # carrying capacity of the system
                rho0 = 1000,    # Initial number of population
                sigma0=0.1,  #Initial standard variation of the population
                x0=0.,
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                d_e=2,      #exponent for the death function
                beta = 0, 
                mu = 1,
                sigma = 0.01,
                tau = 0.1, # transfer rate
                L=4, #space trait X=[-L,L]
                dX=1/100 #discretization of space trait
                )

# Auxiliary part
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
plt.plot(u0)
u= np.empty((nT,nX),dtype=float)
u[0]=u0[0]


def death (x):
    return parameters_HJ['d_r']*x**parameters_HJ['d_e']
Death= (np.vectorize(death,otypes=[np.float64]))(X)


def birth(x):
    return parameters_HJ['b_r']

def mutation_kernel (x):
    return np.exp(-(x/sigma)**2)

def birth_kernel(x,p):
    return mutation_kernel (x)*birth(x)*np.exp(p*x)

Birth= (np.vectorize(birth_kernel,otypes=[np.float64,np.float64]))

def Next_time_HJ(u,parameters_HJ):
    
    grad_u= (u[1:]-u[:-1])/parameters_HJ['dX']
    grad_u=np.insert(grad_u,0,grad_u[0])
    
    
