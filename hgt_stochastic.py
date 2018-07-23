# Import packages here

import numpy as np
import random as rnd

# Auxiliary part
T_max= 10 #Maximal time
dT=0.1 #Step of discretization for time
nT=int(T_max/dT) #number of times
T=[t*dT for t in range(nT)] #list of all times

X_min = 0   # Minimal trait
X_max = 4   # Maximum amount of traits
dX = 0.1      # Step of discretization for time
nX = int((X_max-X_min)/dX) # number of traits
X = [x*dX for x in range(nX)] # list of possible traits

K = 1000    # Maximal capacity of the system
C = 0.5     # Competition
p = 0.03    # Probability of mutation
beta = 0    # 
mu = 1      # 

sigma = 0.1     # variance of the mutation kernel

N_0 = 1000*rnd.gauss(1, K) # Initial number of populations 

N = [[0] * nX for _ in range(nT)] #Intialization of the population. Matrix T*X

def mut_kern(x):
    return rnd.gauss(0, sigma^2)

def b(x):   # Birth rate function
    return 4-x

def d(x):   # Death rate function
    return 1

def death_prob(x, Ntot):
    # nx = number of individuals in a population with a trait x
    # C = Constant competition rate
    return d(x) + Ntot*C

# Compute the probability of birth, death and horizontal transfer events

def tau(x,y):   # function tau
    if (x<y): 
        return 1 
    else:
        return 0
    
def horizontal_transfer(x, y, Ntot):    # HT rate, i.e tau/(beta+mu*Ntot)
    return tau(x,y)/(beta+mu*Ntot)



