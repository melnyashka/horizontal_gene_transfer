import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gc
from datetime import datetime

def horizontal_transfer(x, tau, beta, mu):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return list(map(lambda i: ht_rate*(n_tot-i), range(n_tot)))
    
def Next_Generation(x, parameters):
    b_r, d_r, C, K, sigma, d_e = parameters['b_r'], parameters['d_r'], parameters['C'], parameters['K'], parameters['sigma'],parameters['d_e']
    n_tot = x.size
    if n_tot==0:
        return x
    else:
        beta_birth = np.divide(1,np.repeat(b_r, n_tot))
        beta_death = np.divide(1,d_r*np.power(np.absolute(x),d_e) + n_tot*C/K)
        beta_transfer = np.divide(1,horizontal_transfer(x, tau = parameters['tau'], beta = parameters['beta'], mu = parameters['mu']))
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(beta_transfer)])
        b_mat = (times < parameters['dT'])

        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),
                                       np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))

def build_and_save(Abs, Ord, parameters, path): # function for creating and saving plots
    par_str = '' # create a string of parameters to pass into plots
    for k, v in parameters.items():
        if k == 'N0' or k == 'b_r': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    figure = plt.figure()
    # plt.hist2d(Abs,Ord,bins=3/2*parameters['K'],cmap=plt.cm.jet,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)
    plt.hist2d(Abs, Ord, bins=parameters['K']/2, cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('time')
    plt.ylabel('trait');
    plt.title(par_str)
    # plt.show()
    current_time = datetime.now().time()
    # figure.savefig(str(path +"plot_" + str(current_time)[0:8].replace(':','_')+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    
    figure.savefig(str(path +"plot_" + str(current_time)[0:8]+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    

########################################################################################
###                  Executable part 
########################################################################################

parameters = dict(T_max = 1000, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 10000, # Maximal capacity of the system
                  N0 = 10000,    # Initial number of population
                 sigma0=0.1,  #Initial standard variation of the population
                 beta = 0, 
                 d_e = 2,
                mu = 1,
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                C = 0.4,    # competition
                sigma = 0.01, # mutation variance
                tau = 0.1  # transfer rate
                )

# Let us check the next taus:
tau_i = np.arange(0.05,1.,0.05)

for i in range(len(tau_i)):
    print(str(i) + " and " + str(tau_i[i]))
    parameters['tau'] = tau_i[i]
    X0 = np.random.normal(1, parameters['sigma0'], parameters['N0']) # Initial population
    # X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max
    X = np.sort(X0)

    Abs=[]   
    Ord=[]

    for i in range(int(parameters['T_max']/parameters['dT']-1)):
        for x in X:
            Abs.extend([i*parameters['dT']])
            Ord.extend([x])
        X=Next_Generation(X, parameters)
        if i%100==0: 
            gc.collect()
    build_and_save(Abs, Ord, parameters, path = "/scratch/gene/Figures/") 
 
    
