import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gc
from datetime import datetime
from time import time 

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


# Now we have to create a grid:
def create_grid(parameters):
    T_max, dT, x_mean, sigma0 = parameters["T_max"], parameters["dT"], parameters["x_mean"], parameters["sigma0"]
    Abs = np.arange(0, T_max, dT)
    Ord = np.arange(-1 - x_mean - 5*sigma0, x_mean + 5*sigma0, 0.0001) # here we make a static grid within given limits with the given accuracy (the last number of the function). We need it in order to count the number of elements fitting in the grid, and plot the obtained "discretized" population
    return Abs, Ord

def discretize(x, Ord):
    x_discret = np.vectorize(lambda i: ((Ord[i]<x)&(x<Ord[i+1])).sum())(range(len(Ord) - 1))
    return x_discret

def build_and_save(Abs, Ord, XT, mean_x, len_x, parameters, path): # function for creating and saving plots
    par_str = '' # create a string of parameters to pass into plots
    for k, v in parameters.items():
        if k == 'N0' or k == 'b_r': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    ord_ax = Ord[:-1]
    X_min, X_max = ord_ax[np.all(XT == 0, axis = 0)][0], ord_ax[np.all(XT == 0, axis = 0)][-1]
    xt = XT[:,~np.all(XT == 0, axis=0)]
    fig1.imshow(xt.transpose()[::-1],cmap=plt.cm.jet,aspect='auto',extent=(0,parameters['T_max'],X_min,X_max),  vmin=0, vmax = int(np.max(XT)/2 +1))
    # plt.hist2d(Abs, Ord, bins=parameters['K']/2, cmap=plt.cm.jet)
    fig1.set_xlabel('time')
    fig1.set_ylabel('trait');
    fig1.set_title('Population dynamics')
    # plt.show()
    fig2 = plt.subplot(grid[0,3:])
    plt.plot(Abs[:-1], len_x[:-1])
    fig2.set_title("Population size")
    fig3 = plt.subplot(grid[1,3:])
    plt.plot(Abs[:-1], mean_x[:-1])
    fig3.set_title("Mean trait")
    current_time = datetime.now().time()
    # figure.savefig(str(path +"plot_" + str(current_time)[0:8].replace(':','_')+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    figure.savefig(str(path +"plot_" + str(current_time)[0:8]+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    plt.close() 

########################################################################################
###                  Executable part 
########################################################################################

print("Hello")

parameters = dict(T_max = 2000, # maximal time 
                  dT = 0.1, # Discretization time 
                  K = 10000, # Maximal capacity of the system
                  N0 = 10000,    # Initial number of population
                 sigma0=0.1,  #Initial standard variation of the population
                 beta = 0, 
                 d_e = 2,
                mu = 1,
                b_r = 1,     # birth rate
                d_r = 1,      # death rate
                C = 0.5,    # competition
                sigma = 0.01, # mutation variance
                tau = 0.1,  # transfer rate
                x_mean = 1
                )

# Let us check the next taus:
tau_i = np.arange(0.2,1.,0.05)

for i in range(len(tau_i)):
    print(str(i) + " and " + str(tau_i[i]), flush = True)
    parameters['tau'] = tau_i[i]
    X0 = np.random.normal(parameters["x_mean"], parameters['sigma0'], parameters['N0']) # Initial population
    # X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max
    X = np.sort(X0)
    len_x, mean_x = np.empty([len(Abs)]), np.empty([len(Abs)])
    Abs, Ord = create_grid(parameters)
    XT = np.empty([len(Abs), len(Ord)-1])
    for i in range(int(parameters['T_max']/parameters['dT']-1)):
        start_time = time()
        XT[i] = discretize(X, Ord)
        len_x[i], mean_x[i] = len(X), np.mean(X)
        X=Next_Generation(X, parameters)
        if i%100 == 0: print("That is our "+str(i)+ "-th iteration, it lasted "+str(round(time()-start_time, 3))+"s", flush = True)
    gc.collect()
    build_and_save(Abs, Ord, XT, len_x, mean_x, parameters, path = "/scratch/gene/horizontal_gene_transfer/Mesocentre/Figures/") 
 
    

    
