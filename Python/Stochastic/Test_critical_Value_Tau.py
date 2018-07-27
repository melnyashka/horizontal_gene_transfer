#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 16:14:55 2018

@author: samuelnordmann
"""
import gc
from datetime import datetime


# import stochastic_continuous # Check if it works on your machine! 
parameters = dict(T_max = 1000, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                 sigma0=0.1,  #Initial standard variation of the population
                 x_mean0=0.,
                C = 0.5,    # competition
#                p = 1,      # Probability of mutation
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                d_e=2,      #exponent for the death function
                beta = 0, 
                mu = 1,
                sigma = 0.01,
                tau = 0.17  # transfer rate
                )

# Change some parameters if needed!
# parameters['tau'] = 0.17
# idea: create a grid of parameters we want to check, and then run the experiments inside the loop! 



















X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
X0.sort()
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max


Abs=[]   
Ord=[]


tau_M=1.5
tau_m=0.001
tau_step=0.02

param='d_e'
param_M=4
param_m=1.5
param_step=0.2

for param_value in np.arange(param_m,param_M,param_step):
    parameters[param]=param_value
    print(param+' = '+str(param_value))
    for tau in np.flip(np.arange(tau_m,tau_M,tau_step),0):
        parameters['tau']=tau
        X = X0
        for i in range(int(parameters['T_max']/parameters['dT']-1)):
            X=Next_Generation(X, parameters)
            if X.size<10:#If the population dies, we go out of the loop
                break
        if X.size<10:
            continue
        Ord.append(parameters['tau'])
        Abs.append(parameters[param])
        print(Ord)
        print(Abs)
        break



par_str = '' # create a string of parameters to pass into plots
for k, v in parameters.items():
    if k == 'N0' or k == 'b_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth





figure = plt.figure()
plt.plot(Abs,Ord)
plt.xlabel(param)
plt.ylabel('Tau_critical')
plt.title(par_str)
plt.show()
current_time = datetime.now().time()
figure.savefig(str("Figures/CriticalTau_param="+param+"__"+ str(current_time)[0:8]+".pdf"), bbox_inches='tight')
plt.clf()
plt.gcf().clear()






























X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
X0.sort()
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max


Abs=[]   
Ord=[]


tau_M=1
tau_m=0.001
tau_step=0.01

param='b_r'
param_M=1
param_m=0.1
param_step=0.1

for param_value in np.arange(param_m,param_M,param_step):
    parameters[param]=param_value
    print(param+' = '+str(param_value))
    for tau in np.flip(np.arange(tau_m,tau_M,tau_step),0):
        parameters['tau']=tau
        X = X0
        for i in range(int(parameters['T_max']/parameters['dT']-1)):
            X=Next_Generation(X, parameters)
            if X.size<10:#If the population dies, we go out of the loop
                break
        if X.size<10:
            continue
        Ord.append(parameters['tau'])
        Abs.append(parameters[param])
        print(Ord)
        print(Abs)
        break



par_str = '' # create a string of parameters to pass into plots
for k, v in parameters.items():
    if k == 'N0' or k == 'b_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth





figure = plt.figure()
plt.plot(Abs,Ord)
plt.xlabel(param)
plt.ylabel('Tau_critical')
plt.title(par_str)
plt.show()
current_time = datetime.now().time()
figure.savefig(str("Figures/CriticalTau_param="+param+"__"+ str(current_time)[0:8]+".pdf"), bbox_inches='tight')
plt.clf()
plt.gcf().clear()























X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
X0.sort()
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max


Abs=[]   
Ord=[]


tau_M=1.5
tau_m=0.001
tau_step=0.02

param='sigma'
param_M=0.5
param_m=0.01
param_step=0.05

for param_value in np.arange(param_m,param_M,param_step):
    parameters[param]=param_value
    print(param+' = '+str(param_value))
    for tau in np.flip(np.arange(tau_m,tau_M,tau_step),0):
        parameters['tau']=tau
        X = X0
        for i in range(int(parameters['T_max']/parameters['dT']-1)):
            X=Next_Generation(X, parameters)
            if X.size<10:#If the population dies, we go out of the loop
                break
        if X.size<10:
            continue
        Ord.append(parameters['tau'])
        Abs.append(parameters[param])
        print(Ord)
        print(Abs)
        break



par_str = '' # create a string of parameters to pass into plots
for k, v in parameters.items():
    if k == 'N0' or k == 'b_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth





figure = plt.figure()
plt.plot(Abs,Ord)
plt.xlabel(param)
plt.ylabel('Tau_critical')
plt.title(par_str)
plt.show()
current_time = datetime.now().time()
figure.savefig(str("Figures/CriticalTau_param="+param+"__"+ str(current_time)[0:8]+".pdf"), bbox_inches='tight')
plt.clf()
plt.gcf().clear()









































gc.collect()



