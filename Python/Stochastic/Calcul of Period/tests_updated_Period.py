# import stochastic_continuous # Check if it works on your machine! 
parameters = dict(T_max = 2000, # maximal time 
                dT = 0.1, # Discretization time 
                K = 1000, # Maximal capacity of the system
                N0 = 1000,    # Initial number of population
                sigma0=0.1,  #Initial standard variation of the population
                x_mean0=0.,
                C = 0.5,    # competition
#                p = 0.03,      # Probability of mutation
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                d_e=2,      #exponent for the death function
                beta = 0, 
                mu = 1,
                sigma = 0.01,
                tau = 0.22 # transfer rate
                )

# Change some parameters if needed!
# parameters['tau'] = 0.17
# idea: create a grid of parameters we want to check, and then run the experiments inside the loop! 

X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max
X = np.sort(X0)

Abs=[]#Save all the individuals
Ord=[]
Abs_mean=[]#Save the mean trait
Ord_mean=[]
Abs_tot=[]#Save the total size of the population
Ord_tot=[]
for i in range(int(parameters['T_max']/parameters['dT']-1)):
    if i%1000==0:
        print('T= '+str(i*parameters['dT']))
    for x in X:
        Abs.append(i*parameters['dT'])
        Ord.append(x)
    Abs_mean.append(i*parameters['dT'])
    Ord_mean.append(np.mean(X))
    Abs_tot.append(i*parameters['dT'])
    Ord_tot.append(X.size)
    X=Next_Generation(X, parameters)
    
Abs=np.array(Abs)
Ord=np.array(Ord)
Abs_mean=np.array(Abs_mean)
Ord_mean=np.array(Ord_mean)
Abs_tot=np.array(Abs_tot)
Ord_tot=np.array(Ord_tot)
#build_and_save(Abs, Ord, parameters) # build and save a plot in folder Figures in home directory (you must create it first!)
gc.collect()
