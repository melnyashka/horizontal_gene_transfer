# Import packages here

import numpy as np

# Auxiliary part
T_max= 10 #Maximal time
dT=0.001 #Step of discretization for time
nT=int(T_max/dT) #number of times
T=[t*dT for t in range(nT)] #list of all times

dT2=dT #interval of time for which events are taken into account

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


#INITIAL TIME
N0=[0 for x in X]
Ntot0=1000#initial size of population

m=1#Initial law of repartition. Gaussian centered at m with variance sigma0^2
sigma0=0.1
X_weighted= list(map(lambda z: np.exp(-((m-z)/sigma0)**2), X))
S=sum(X_weighted)
Initial_repartition=list(map(lambda z: z/S, X_weighted))#law of repartition
    
x=np.random.choice(X, size=Ntot0, replace=True, p=Initial_repartition)#choice of Ntot random variable with respect to the law
for y in x:
    N0[resc_x(y)]+=1#creation of the vector of the initial population






#################################################

#FUNCTION/PARAMETERS

def b(x):   # Birth rate function
    return 4-x

def d(x):   # Death rate function
    return 1

def death_prob(x, Ntot):
    # Ntot = total size of the population
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

























#TESTS
N[0]=N0


