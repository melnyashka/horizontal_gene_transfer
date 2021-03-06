
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime

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
    
    pre_init_values = dict(
        f = f,
        T = T,
        nX = nX,
        nT = nT,
        X = X, 
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
    if rho>np.power(10.,-7):
        transfer_part = parameters['tau']/rho*parameters['dX']*np.fromiter((np.sum(f[:i])-np.sum(f[i:]) for i in I),float)
    else :
        transfer_part = 0
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
    par_str='T_max='+str(parameters['T_max'])+', dT='+str(parameters['dT'])+', tau='+str(parameters['tau'])
    # Now we have to compute the mean trait and the population size at each time! 
    X, nT = pre_init_values['X'], pre_init_values['nT']
    sum_f = np.sum(f, axis = 1)
    computed_mean = np.divide(np.sum(X*f, axis = 1), sum_f)
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    fig1.imshow(f.transpose()[::-1],cmap=plt.cm.jet, aspect = 'auto', extent = (0,parameters['T_max'], parameters['X_min'],parameters['X_max']), 
            vmin = 0, vmax = 4)
    fig1.set_xlabel('time')
    fig1.set_ylabel('trait')
    fig1.set_title('Population dynamics')
    fig2 = plt.subplot(grid[0,3:])
    plt.plot(np.arange(nT)/parameters['dT'], sum_f)
    fig2.set_title("Population size")
    fig3 = plt.subplot(grid[1,3:])
    plt.plot(np.arange(nT)/parameters['dT'], computed_mean)
    fig3.set_title("Mean trait")
    current_time = datetime.now().time()
    figure.savefig(str(path +"plot_"+str(current_time)[0:8]+".pdf"), bbox_inches ='tight')
    plt.show()
    # plt.close()

#######################################    
####### EXECUTABLE PART ###############
#######################################

if __name__ == "__main__":

    parameters = dict(T_max = 1000, # maximal time 
                      dT = 0.001, # Discretization time 
                      sigma0 = 0.01,  #Initial standard variation of the population
                      x_mean0 = 0.,
                      C = 0.5,    # competition
    #                 p = 0.03,      # Probability of mutation
                      b_r = 1,     # birth rate
                      d_r = 1,      # death rate
                      d_e = 2,   #exponetial power
                      beta = 0, 
                      mu = 1,
                      sigma = 0.01,
                      tau = 0.1,  # transfer rate
                      X_min = -0.5, #length of the numerical interval of traits (for PDE!)
                      X_max=1.5,
                      dX = 0.01, #discretization of the space of traits
                      eps = 1
                      )
        
    # Loop over given parameters
    param='tau' # here we check the value
    param_m=0.95
    param_M=1
    param_step=0.01
    J=[0.4,0.9,1.2,0.02]
    for j in J:
        pre_init_values = Pre_Initialization_PDE(parameters) 
        f, nT = pre_init_values['f'], pre_init_values['nT']
        parameters[param]=j
        print(' !!!!!!!!!!!!!! '+str(param)+' = '+str(j), flush = True)    
        gc.collect()
        #SIMULATION !!!!!!!!!!!!!!!!!!
        for i in range(nT-1):
            f[i+1] = Next_Generation_PDE(f[i],parameters, pre_init_values)
            if i%1000==0:
                print('T= '+str(i*parameters['dT']), flush = True)
        build_and_save_PDE(f, parameters, pre_init_values, "Figures/")
