# Import packages here

import numpy as np
import scipy
from scipy import stats
import pandas as pd
import random


# Auxiliary part
T_max= 10 #Maximal time
dT=0.1 #Step of discretization for time
nT=int(T_max/dT) #number of times
T=[t*dT for t in range(nT)] #list of all times

dT2=dT #interval of time for which events are taken into account

X_min = 0   # Minimal trait
X_max = 4   # Maximum amount of traits
dX = 0.1      # Step of discretization for time
nX = int((X_max-X_min)/dX) # number of traits
X = [x*dX for x in range(nX)] # list of possible traits

K = 1000    # Maximal capacity of the system
C = 0.05     # Competition
p = 0.03    # Probability of mutation
beta = 0    # 
mu = 1      # 

sigma = 0.01     # variance of the mutation kernel


Xtot0=1000  #initial size of population
X0 = np.random.normal(loc=0.0, scale=1.0, size=Xtot0)
X0.sort()



X=X0
Bb=np.random.choice([True,False], size=Xtot0) 
Bd=np.random.choice([True,False], size=Xtot0) 
Bht=np.random.choice([True,False], size=Xtot0) 



#X_add_birth=np.random.normal(loc=X[Bb], scale=sigma, size=None)
#X_add_ht= np.vectorize(lambda i: np.random.choice(X[(i+1):]))(np.arange(X.size)[Bht][:-1])
#X_keep=X[np.logical_not(np.logical_or(Bd,Bht))]
#X_new=np.concatenate((X_keep,X_add_birth,X_add_ht))


np.sort(np.concatenate((X[np.logical_not(np.logical_or(Bd,Bht))],np.random.normal(loc=X[Bb], scale=sigma, size=None),np.vectorize(lambda i: np.random.choice(X[(i+1):]))(np.arange(X.size)[Bht][:-1]))))

