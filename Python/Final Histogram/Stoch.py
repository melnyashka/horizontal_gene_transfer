import numpy as np
import matplotlib.pyplot as plt
import gc
from datetime import datetime
"""
def horizontal_transfer(X,x,n_tot):
    return np.sum(np.exp(-parameters_S['tau']/n_tot*np.tanh((x-X)/parameters_S['delta'])))/n_tot


def HT(X):
    n_tot = len(X)
    return np.vectorize((lambda x: horizontal_transfer(X,x,n_tot)), otypes=[float])(X)
"""

def horizontal_transfer(x, tau):
    # Do transfer in an already sorted list!!!
    # x = x.sort()
    n_tot = len(x)
    ht_rate = tau/n_tot
    return list(map(lambda i: ht_rate*(n_tot-i), range(n_tot)))



def Next_Generation(x, parameters):
    b_r, d_r, C, K, sigma, d_e = parameters['b_r'], parameters['d_r'], parameters['C'], parameters['K'], parameters['sigma'],parameters['d_e']
    n_tot = x.size
    if n_tot==0:
        return x
    else:
        beta_birth = np.divide(1,np.repeat(b_r, n_tot))
        beta_death = np.divide(1,d_r*np.power(np.absolute(x),d_e) + n_tot*C/K)
        #beta_transfer = HT(X)
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(HT(X))])
        b_mat = (times < parameters['dT'])

        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),
                                       np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))



def build_and_save(Abs, Ord, parameters_S, path): # function for creating and saving plots
    par_str = '' # create a string of parameters_S to pass into plots
    for k, v in parameters_S.items():
        if k == 'N0' or k == 'b_r': 
            smth = ",\n" 
        else: 
            smth = ", "
        par_str += k + "=" + str(v) + smth
    figure = plt.figure()
    plt.hist2d(Abs,Ord,bins=3/2*parameters_S['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters_S['K'],cmin=parameters_S['K']/100)
    #plt.hist2d(Abs, Ord, bins=parameters_S['K']/2, cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('time')
    plt.ylabel('trait');
    plt.title(par_str)
    plt.show()
    current_time = datetime.now().time()
    figure.savefig(str(path +"plot_" + str(current_time)[0:8].replace(':','_')+".pdf"), bbox_inches='tight') # Possibly different delimeter on Linux and Windows!
    


# import stochastic_continuous # Check if it works on your machine! 
parameters_S = dict(T_max = 100, # maximal time 
                dT = 0.01, # Discretization time 
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
                tau = 0.4, # transfer rate
                delta=0.01
                )

#Initialization
X0 = np.random.normal(parameters_S['x_mean0'], parameters_S['sigma0'], parameters_S['N0']) # Initial population
#X = [None]*int(parameters_S['T_max']/parameters_S['dT'])  # history of all populations up to time T_max
X = np.sort(X0)
Size_S=[X0.size/parameters_S['K']]
Mean_S=[np.mean(X0)]

#Initialization of histogram
X_max=parameters_HJ['X_max']
X_min=parameters_HJ['X_min']
bins=int((X_max-X_min)/parameters_HJ['dX'])
rg=(X_min,X_max)

Abs_Stoch=np.histogram(X0,bins,rg,True)[1]
U_S_0=np.histogram(X0,bins,rg,True)[0]*parameters_S['N0']/parameters_S['K']

U_S=np.zeros((int(parameters_S['T_max']/parameters_S['dT']),bins))
U_S[0]=U_S_0


nbr_iteration= int(parameters_S['T_max']/parameters_S['dT'])
for i in range(nbr_iteration-1):
    if i%500==0:
        print(str(i)+'th iteration out of '+str(nbr_iteration)+'.')
        print(X.size/parameters_S['K'])
        plt.plot(U_S[i])
        plt.show()
    X=Next_Generation(X, parameters_S)
    U_S[i+1]=np.histogram(X,bins,rg,False)[0]/parameters_S['K']
    Size_S=Size_S+[X.size/parameters_S['K']]
    Mean_S=Mean_S+[np.mean(X)]
#build_and_save(Abs, Ord, parameters_S, path = "Figures/") # build and save a plot in folder Figures in home directory (you must create it first!)
#gc.collect()






