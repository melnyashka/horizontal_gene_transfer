import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime
from time import time

def horizontal_transfer(x, tau, beta, mu):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/(beta+mu*n_tot)
    return np.vectorize(lambda i: ht_rate*(n_tot-i))(range(n_tot))
    
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
    T_max, dT, x_mean, sigma0 = parameters["T_max"], parameters["dT"], parameters["x_mean0"], parameters["sigma0"]
    Abs = np.arange(0, T_max, dT)
    Ord = np.arange(x_mean - 1.5, x_mean + 1.5, 0.001) # here we make a static grid within given limits with the given accuracy (the last number of the function). We need it in order to count the number of elements fitting in the grid, and plot the obtained "discretized" population
    return Abs, Ord

def discretize(x, Ord):
    x_discret = np.vectorize(lambda i: ((Ord[i]<x)&(x<Ord[i+1])).sum())(range(len(Ord) - 1))
    return x_discret

def build_and_save(Abs, Ord, XT, len_x, mean_x, parameters, path): # function for creating and saving plots
    par_str = '' # create a string of parameters to pass into plots
    for k, v in parameters.items():
        if k == 'N0' or k == 'b_r': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    par_str = 'T_max='+str(parameters['T_max'])+', dT='+str(parameters['dT'])+', tau='+str(parameters['tau'])+', N0='+str(parameters['N0'])
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    ord_ax = Ord[:-1]
    X_min, X_max = ord_ax[~np.all(XT == 0, axis = 0)][0], ord_ax[~np.all(XT == 0, axis = 0)][-1]
    xt = XT[:,~np.all(XT == 0, axis=0)]
    fig1.imshow(xt.transpose()[::-1],cmap=plt.cm.jet,aspect='auto',extent=(0,parameters['T_max'],X_min,X_max),  vmin=0, vmax = int(np.max(XT)/2 +1))
    # plt.hist2d(Abs, Ord, bins=parameters['K']/2, cmap=plt.cm.jet)
    fig1.set_xlabel('time')
    fig1.set_ylabel('trait');
    fig1.set_title('Population dynamics')
    # plt.show()
    fig2 = plt.subplot(grid[0,3:])
    plt.plot(Abs[:-1], len_x[:-1])
    fig2.set_title("N/K (norm. popul. size)")
    fig3 = plt.subplot(grid[1,3:])
    plt.plot(Abs[:-1], mean_x[:-1])
    fig3.set_title("Mean trait")
    current_time = datetime.now().time()
    # figure.savefig(str(path +"plot_" + str(current_time)[0:8].replace(':','_')+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    figure.savefig(str(path +"plot_" + str(current_time)[0:8]+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    np.savetxt(str(path+"mean_trait"+str(current_time)[0:8]+".txt"), mean_x, delimiter=',')
    np.savetxt(str(path+"population_size_"+str(current_time)[0:8]+".txt"), len_x, delimiter=",") 
    plt.close() 

def Pre_Initialization_PDE(parameters):
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= parameters['X_min'], parameters['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
       
    f0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/parameters['sigma0']*parameters['eps'])# initial density 
    rho0=(parameters['b_r'])/(parameters['C'])
    f0=f0/(np.sum(f0)*parameters['dX'])*rho0

    f = np.empty([nT, nX])        #densities for all times
    f[0]=f0
    
    # Computing constant death and mutation kernels
    Death = parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
    Mutation_kernel=np.exp(-(X/(parameters['sigma']*parameters['eps']))**2*1/2)*parameters['dX']
    Mutation_kernel=parameters['b_r']*Mutation_kernel/np.sum(Mutation_kernel)
    X0 = np.random.normal(parameters["x_mean0"], parameters['sigma0'], parameters['N0']) 
    
    pre_init_values = dict(
        f = f,
        T = T,
        nX = nX,
        nT = nT,
        X = X, # vector of traits
        X0 = X0, # initial population
        Death = Death, 
        Mutation_kernel = Mutation_kernel
        )
    return pre_init_values

def Next_Generation_PDE(f,parameters, pre_init_values):
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= parameters['X_min'], parameters['X_max']
    Death, Mutation_kernel = pre_init_values['Death'], pre_init_values['Mutation_kernel']

    X_min_new=int(-X_min/dX)          # Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-2*X_min+X_max)/dX)
    I=range(pre_init_values['nX'])    # all the indexes
    
    rho = np.sum(f)*parameters['dX'] # I don't like the idea of comparing float with an integer, but okay
    
    death_term = Death + rho*parameters['C']
    birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
    transfer_part = parameters['tau']/rho*parameters['dX']*np.fromiter((np.sum(f[:i])-np.sum(f[i:]) for i in I),float)
    #new_f = np.maximum(0,f + dT/parameters['eps']*(np.multiply(f,(-death_term+ transfer_part))+birth_part))
    new_f = f + dT/parameters['eps']*(np.multiply(f,(-death_term + transfer_part))+birth_part)
    return new_f

def build_and_save_PDE(f, parameters, pre_init_values, path):
    par_str = '' # create a string of parameters to pass into pots
    for k,v in parameters.items():
        if k == 'N0' or k=='b_r':
            smth =",\n"
        else:
            smth=", "
        par_str+=k+"="+str(v)+smth
    # Now we have to compute the mean trait and the population size at each time! 
    X, nT = pre_init_values['X'], pre_init_values['nT']
    X_min, X_max, dX = parameters['X_min'], parameters['X_max'], parameters['dX']
    sum_f = np.sum(f, axis = 1)*dX
    computed_mean = np.divide(np.sum(X*f, axis = 1), sum_f)
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    fig1.imshow(f.transpose()[::-1],cmap=plt.cm.jet, aspect = 'auto', extent = (0,parameters['T_max'], X_min, X_max), 
            vmin = 0, vmax = 4)
    fig1.set_xlabel('time')
    fig1.set_ylabel('trait')
    fig1.set_title('Population dynamics')
    fig2 = plt.subplot(grid[0,3:])
    plt.plot(np.arange(nT)*parameters['dT'], sum_f)
    fig2.set_title("Rho")
    fig3 = plt.subplot(grid[1,3:])
    plt.plot(np.arange(nT)*parameters['dT'], computed_mean)
    fig3.set_title("Mean trait")
    current_time = datetime.now().time()
    figure.savefig(str(path +"plot_"+str(current_time)[0:8]+".pdf"), bbox_inches ='tight')
    np.savetxt(str(path+"mean_trait_"+str(current_time)[0:8]+".txt"), computed_mean, delimiter = ',')
    np.savetxt(str(path+"population_size_"+str(current_time)[0:8]+".txt"), sum_f, delimiter = ',')
	# plt.show()


########################################################################################
###                  Executable part 
########################################################################################

if __name__ == "__main__":

    parameters = dict(T_max = 600, # maximal time 
                      dT = 0.1, # Discretization time 
                      sigma0 = 0.01,  #Initial standard variation of the population
                      x_mean0 = 0.,
                      K = 10000, # Maximal capacity of the system
                      N0 = 10000,    # Initial number of population
                      C = 0.5,    # competition
    #                 p = 0.03,      # Probability of mutation
                      b_r = 1,     # birth rate
                      d_r = 1,      # death rate
                      d_e = 2,   #exponetial power
                      beta = 0, 
                      mu = 1,
                      sigma = 0.01,
                      tau = 0.3,  # transfer rate
                      X_min = -0.5, #length of the numerical interval of traits (for PDE!)
                      X_max=1,
                      dX = 0.5, #discretization of the space of traits
                      eps = 1
                      )
    # Let us check the next taus:
    #tau_i = np.arange(0.2,0.55,0.025)
    tau_i = [0.46,0.46,0.46,0.46]
    for i in range(len(tau_i)):
        print(str(i) + " and " + str(round(tau_i[i],3)), flush = True)
        parameters['tau'] = tau_i[i] 
        pre_init_values = Pre_Initialization_PDE(parameters) 
        f, nT, X0 = pre_init_values['f'], pre_init_values['nT'], pre_init_values['X0']
        X = np.sort(X0)
        Abs, Ord = create_grid(parameters)
        len_x, mean_x = np.empty([len(Abs)]), np.empty([len(Abs)])
        XT = np.empty([len(Abs), len(Ord)-1])
        for i in range(nT-1):
            start_time = time()
            XT[i], len_x[i], mean_x[i] = discretize(X, Ord), len(X)/parameters['K'], np.mean(X)
            X = Next_Generation(X, parameters)
            #f[i+1] = Next_Generation_PDE(f[i],parameters, pre_init_values)
            if i%100 == 0: print("That is our "+str(i)+ "-th iteration, it lasted "+str(round(time()-start_time, 3))+"s", flush = True)
        #gc.collect()
        build_and_save(Abs, Ord, XT, len_x, mean_x, parameters, path = "NewFigures/") 
        #build_and_save_PDE(f, parameters, pre_init_values, path = "/scratch/gene/horizontal_gene_transfer/Mesocentre/Figures/PDE/")
    #path = "/scratch/gene/horizontal_gene_transfer/Mesocentre/Figures/Stoch/"
