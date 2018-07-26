#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 16:14:55 2018

@author: samuelnordmann
"""

# import stochastic_continuous # Check if it works on your machine! 
parameters = dict(T_max = 1000, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                 sigma0=0.1,  #Initial standard variation of the population
                 x_mean0=0.,
                C = 0.5,    # competition
#                p = 0.03,      # Probability of mutation
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                beta = 0, 
                mu = 1,
                sigma = 0.01,
                tau = 0.17  # transfer rate
                )

# Change some parameters if needed!
# parameters['tau'] = 0.17
# idea: create a grid of parameters we want to check, and then run the experiments inside the loop! 

X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max
X = np.sort(X0)

Abs=[]   
Ord=[]


tau_M=3
tau_m=0.001
tau_step=0.05

param='b_r'
param_M=10
param_m=1
param_step=1


for tau in np.flip(np.arange(tau_m,tau_M,tau_step),0):
    parameters['tau']=tau
    for i in range(int(parameters['T_max']/parameters['dT']-1)):
        X=Next_Generation(X, parameters)
        if X.size<10:#If the population dies, we go out of the loop
            break
    if X.size<10:
        continue
    Ord.append(parameters['tau'])
    Abs.append(parameters[param])
    break





