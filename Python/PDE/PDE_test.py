parameters = dict(T_max = 100, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                  sigma0 = 0.1,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.01,
                  tau = 0.17,  # transfer rate
                  L = 3, #length of the numerical interval of traits (for PDE!)
                  dX = 0.001, #discretization of the space of traits
                  eps = 1
                  )

L, dX, T_max, dT = parameters['L'], parameters['dX'], parameters['T_max'], parameters['dT']
X_min, X_max= -L, L
nX = int((X_max-X_min)/dX) # number of traits
X = np.arange(X_min,X_max,dX) # list of possible traits
F0 = np.exp(-np.power(X-parameters['x_mean0'],2)/parameters['sigma0']*parameters['eps']) # initial density 
XT = np.empty([int(T_max/dT), nX])
XT[0] = F0

for i in range(int(T_max/dT)-1):
    XT[i + 1] = Next_Generation_PDE(XT[i], parameters)
