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

for i in range(int(parameters['T_max']/parameters['dT']-1)):
    for x in X:
        Abs.extend([i*parameters['dT']])
        Ord.extend([x])
    X=Next_Generation(X, parameters)
    
build_and_save(Abs, Ord, parameters) # build and save a plot in folder Figures in home directory (you must create it first!)
gc.collect()
