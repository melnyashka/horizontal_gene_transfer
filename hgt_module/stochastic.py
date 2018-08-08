import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gc
from datetime import datetime
from time import time 

def horizontal_transfer(x, tau, beta, mu):
    """ 
    Function, which applies a horizontal transfer kernel to an ALREADY sorted list (for convention). Basically, we just spare time by applying a heaviside function to every pair of individuals in a population. 
    returns an array of elements
    """
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return np.vectorize(lambda i: ht_rate*(n_tot-i))(range(n_tot))
    
def Next_Generation(x, parameters):
    """
    Next_Generation(x, parameters) returns a vector of population for a next time step. Input is the vector of population on the previous step and the parameters of the model, passed as a dictionary.
    Important note: resulting vector is NOT of the same size as x, since we consider continuous space of traitsand the population is described only by a vector of traits (not by number of individuals in each trait, as in discrete model)
    """
    b_r, d_r, C, K, sigma, d_e = parameters['b_r'], parameters['d_r'], parameters['C'], parameters['K'], parameters['sigma'],parameters['d_e']
    n_tot = x.size # Compute size of population
    if n_tot==0:
        return x # We just spare all the numerical procedure if the population goes extinct
    else:
        beta_birth = np.divide(1,np.repeat(b_r, n_tot)) 
        beta_death = np.divide(1,d_r*np.power(np.absolute(x),d_e) + n_tot*C/K)
        beta_transfer = np.divide(1,horizontal_transfer(x, tau = parameters['tau'], beta = parameters['beta'], mu = parameters['mu']))
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(beta_transfer)]) # compute arrival times of each events
        b_mat = (times < parameters['dT']) # create boolean matrix by excluding all the arrival times which are bigger than the time step
        # We return the sorted array of population by adding all the "newborns", deleting all the dead bodies, changing traits of randomly selected individuals according to a horizontal transfer rate
        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),
                                       np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))

def create_grid(parameters):
    """
    here we make a static grid within given limits with the given accuracy (the last number of the function). We need it in order to count the number of elements fitting in the grid, and plot the obtained "discretized" population. We need exclusively for plotting reason, all the computations are done for continuous model!
    """
    T_max, dT, x_mean, sigma0 = parameters["T_max"], parameters["dT"], parameters["x_mean"], parameters["sigma0"]
    Abs = np.arange(0, T_max, dT)
    Ord = np.arange(-1 - x_mean - 5*sigma0, x_mean + 5*sigma0, 0.001 
    return Abs, Ord

def discretize(x, Ord):
    """
    Input: population, described by continuous traits, and a pre-created grid
    Output: vector of size len(Ord)-1, counting number of individuals in x, fitting in each cell of grid
    """
    x_discret = np.vectorize(lambda i: ((Ord[i]<x)&(x<Ord[i+1])).sum())(range(len(Ord) - 1))
    return x_discret

def build_and_save(Abs, Ord, XT, len_x, mean_x, parameters, path): 
    """
    function for creating and saving plots. Saves plots in pdf format, where one big plot on the left describes the evolution of population over time, and two small plots depicture the evolution of the population (normalized by K), and the mean trait. 
    """
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
    figure.savefig(str(path +"plot_" + str(current_time)[0:8]+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    plt.close() 
