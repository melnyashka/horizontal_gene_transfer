import matplotlib 
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gc 
from datetime import datetime
from time import time

def Pre_Initialization_HJ(parameters):
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max= parameters['X_min'], parameters['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
       
    # u0 = np.exp(-np.power(np.absolute(X-parameters['x_mean0']),2)/((1+np.power(X,2))/2*(parameters['sigma0']*parameters['eps'])**2) # initial density 
    u0 = -((np.abs(X)<=1)*(np.power(X,2)/2)+(np.abs(X)>1)*(np.abs(X)-1/2))
    # Helene's trick:
    # u=(abs(x)<=1).*x.^2/2+(abs(x)>1).*(abs(x)-1/2) # first we initialize u
    rho0=np.sum(np.exp(u0/parameters['eps']))*dX  # then we initialize rho
    f0 = np.exp(u0/parameters['eps']) # then we initialize f

    f = np.empty([nT, nX])        #densities for all times
    U = np.empty([nT, nX])        #
    U[0]=u0 # matrix for u
    f[0]=f0 # matrix for f
    Rho=np.empty([nT])
    Rho[0] = rho0
    
    # Computing constant death and mutation kernels
    Death = parameters['d_r']*np.power(np.absolute(X),parameters['d_e'])
    
    pre_init_values = dict(
        f = f, # maybe we don't even need this one!
        U = U,
        T = T,
        nX = nX,
        nT = nT,
        X = X, # vector of traits
        Death = Death, 
        Rho = Rho,
        # Mutation_kernel = Mutation_kernel
        )
    return pre_init_values

def Next_Generation_AP(u, rho, parameters, pre_init_values):
    ### That's our draft for the AP scheme
    dX, T_max, dT = parameters['dX'], parameters['T_max'], parameters['dT']
    X_min, X_max, eps = parameters['X_min'], parameters['X_max'],  parameters['eps']
    Death, X = pre_init_values['Death'], pre_init_values['X']

    # now we define a grid to evaluate an integral, cutting all the unnecessary values
    ymax = 8
    Ny=int(2*ymax/dX)
    dy=2*ymax/Ny
    y=np.arange(-ymax, ymax, dy)
    Xy = np.ones([len(y),1])*X # to be checked! 
    Yx = (np.ones([len(X), 1])*y).T

    interp_grid = Xy+eps*Yx
    interp_grid_b = interp_grid*(interp_grid>=X_min)*(interp_grid<=X_max)+X_min*(interp_grid<X_min)+X_max*(interp_grid>X_max)
    u_interp = np.interp(interp_grid_b, X, u)

    # Neumann conditions check
    u_sq=u*np.ones([len(y),1])
    u_sqp1, u_sqm1 = np.empty(np.shape(u_sq)), np.empty(np.shape(u_sq))
    u_sqp1[-1], u_sqm1[0] = 2*u_sq[-1]-u_sq[-2], 2*u_sq[0]-u_sq[1]
    u_sqp1[:-1], u_sqm1[1:] = u_sq[1:], u_sq[:-1]

    # Writing right and left slope
    left_s = (u[1]-u[0])/dX
    right_s = (u[-1]-u[-2])/dX
    u_left = left_s*(interp_grid-X_min)+u[0]
    u_right = right_s*(interp_grid - X_max)+u[-1]
    
    test=eps*np.abs(Yx)
    flux = (test>=dX)*(interp_grid>=X_min)*(interp_grid<=X_max)*(u_sq-u_interp)/eps     +(interp_grid<X_min)*(u_sq-u_left)/eps +(interp_grid>X_max)*(u_sq-u_right)/eps     + (test<dX)*(test>0)*((interp_grid>=X_min)*(interp_grid<=X_max)*                          ((Yx>0)*(u_sq-u_sqm1)*Yx/dX +(Yx<0)*(u_sqp1-u_sq)*Yx/dX)+                          (interp_grid<X_min)*left_s*Yx + (interp_grid>X_max)*right_s*Yx)+(test==0)*0
    
    # Write birth and transfer kernel
    min_u = np.min(u)
    transfer = parameters['tau']*np.vectorize(lambda y: (np.heaviside(X-y,1)-np.heaviside(y-X,1)))(y in X) # can be replaced by arctan
    rho_u=np.sum(np.exp(-(u-min_u)/eps))*dX
    fsurrhou=np.exp(-(u-min_u)/eps)/rho_u
    T = transfer*np.ones([len(X),1])*fsurrhou
    T = np.sum(T,axis=0)*dX # that's the transfer term

    mut_kern = np.exp(-np.power(Yx,2)/2)
    cste=np.sum(mut_kern)*dy
    m_sq = mut_kern/cste*np.ones([len(X)])

    H = parameters['b_r']*np.exp(flux - np.power(np.ones([len(X)])*Yx,2)/2)/cste*m_sq #*m(x) add a gaussian mutation kernel!
    H = np.sum(H,axis=0)*dy # that's the birth term

    A = u + dT*(pre_init_values['Death'] - H - T) # that's the birth-death-transfer term
    min_A = np.min(A)
    C = eps*np.log(dX)- min_A+eps*np.log(np.sum(np.exp(-(A-min_A)/eps)))

    func = lambda x: C - eps*np.log(x)-dT*x
    invd_func = lambda x: -x/(eps+dT*x)

    # Compute rho:
    tol=np.power(10.,-13)
    err = 10.
    i = 0
    r = np.sum(np.exp(u/eps))*dX # "dummy" rho

    while err>tol:
        i+=1
        newr=r-invd_func(r)*func(r)
        err=np.abs(newr-r)
        r=newr
        if i>1000:
            break
    rho = r
    u = dT*rho + A
    return u, rho

def build_and_save_HJ(u, parameters, pre_init_values, path):
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
    f = np.exp(u/parameters['eps'])
    sum_f = np.sum(f, axis = 1)*dX
    computed_mean = np.divide(np.sum(X*f, axis = 1), sum_f)
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    fig1.imshow(u.transpose()[::-1],cmap=plt.cm.jet, aspect = 'auto', extent = (0,parameters['T_max'], X_min, X_max), 
            vmin = 0, vmax = np.max(u))
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
    plt.show()


#######################################    
####### EXECUTABLE PART ###############
#######################################

parameters = dict(T_max = 10, # maximal time 
                  dT = 0.0001, # Discretization time 
                  sigma0 = 0.01,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.,
                  tau = 0.4,  # transfer rate
                  X_min = -1, #length of the numerical interval of traits (for PDE!)
                  X_max=1,
                  dX = 0.01, #discretization of the space of traits
                  eps = 0.001
                  )

pre_init_values = Pre_Initialization_HJ(parameters) 
U, Rho, nT = pre_init_values['U'], pre_init_values['Rho'], pre_init_values['nT']
for i in range(nT-1):
    U[i+1], Rho[i+1] = Next_Generation_AP(U[i], Rho[i], parameters, pre_init_values)
    if i%100==0: 
        print(str(i)+"-th iteration")

