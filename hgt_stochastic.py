# Import packages here

import numpy as np


# Auxiliary part

X_min = 0   # Minimal trait
X_max = 4   # Maximum amount of traits

K = 1000    # Maximal capacity of the system
C = 0.5     # Competition
p = 0.03    # Probability of mutation

sigma = 0.1     # variance of the mutation kernel

N_0 = 1000 # Initial number of populations 

def b(x):   # Birth rate function
    return 4-x


def d(x):   # Death rate function
    return 1

# Compute the probability of birth, death and horizontal transfer events

def death_prob(x, nx):
    # nx = number of individuals in a population with a trait x
    # C = Constant competition rate
    return d(x) + nx*C

def horizontal_transfer(x, y, nu):
    # x = "donor" trait
    # y = "recepient" trait
    # nu = 
    return 0 
