# %load 'horizontal_gene_transfer/Python/Hamilton_Jacobi/Hamilton_Jacobi_APscheme.py'
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
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
    u0 = -((np.abs(X-parameters['x_mean0'])<=1)*(np.power(X-parameters['x_mean0'],2)/(2*parameters['sigma0']**2))
           +(np.abs(X-parameters['x_mean0'])>1)*(np.abs(X-parameters['x_mean0'])/parameters['sigma0']**2-1/(2*parameters['sigma0']**2)))
    #u0= -np.divide(np.power(X-parameters['x_mean0'],2),1+2*np.power(X-parameters['x_mean0'],2))/(2*parameters['sigma0']**2)
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
    ymax = 8*parameters['sigma']
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
    u_right = right_s*(interp_grid-X_max)+u[-1]
    
    test=eps*np.abs(Yx)
    flux = (test>=dX)*((interp_grid>=X_min)*(interp_grid<=X_max)*(u_interp-u_sq)/eps     +(interp_grid<X_min)*(u_left-u_sq)/eps +(interp_grid>X_max)*(u_right-u_sq)/eps)     + (test<dX)*(test>0)*((interp_grid>=X_min)*(interp_grid<=X_max)*                          ((Yx>0)*(u_sqm1-u_sq)*Yx/dX +(Yx<0)*(u_sq-u_sqp1)*Yx/dX)+                          (interp_grid<X_min)*left_s*Yx + (interp_grid>X_max)*right_s*Yx)+(test==0)*0
    
    # Write birth and transfer kernel
    max_u = np.max(u)
    rho_u=np.sum(np.exp((u-max_u)/eps))*dX
    fsurrhou=np.exp((u-max_u)/eps)/rho_u
    
    X_large=np.arange(2*X_min,2*X_max,dX)
    #transfer_kernel=  parameters['tau']*(np.arctan(np.divide(X_large,parameters['delta']))-np.arctan(-np.divide(X_large,parameters['delta'])))*dX
    transfer_kernel = parameters['tau']*(np.heaviside(X_large,1)-np.heaviside(-X_large,1))*dX
    
    X_min_new=int(-2*X_min/dX)#Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-3*X_min+X_max)/dX)
    T = np.convolve(transfer_kernel,fsurrhou)[X_min_new:X_max_new]# that's the transfer term


    mut_kern = np.exp(-np.power(Yx,2)/(2*parameters['sigma']**2))
    cste=np.sum(mut_kern[:,0])*dy
    m_sq = mut_kern/cste*np.ones([len(X)])

    H = parameters['b_r']*np.exp(flux - np.power(np.ones([len(X)])*Yx,2)/2)/cste #*m_sq #*m(x) add a gaussian mutation kernel!
    H = np.sum(H,axis=0)*dy # that's the birth term

    A = u + dT*(-pre_init_values['Death'] + H + T) # that's the birth-death-transfer term
    max_A = np.max(A)
    E = eps*np.log(dX) + max_A + eps*np.log(np.sum(np.exp((A-max_A)/eps)))

    func = lambda x: E - eps*np.log(x)-parameters['C']*dT*x#the unknown rho is the root of this function
    invd_func = lambda x: -x/(eps+dT*x)#inverse of the derivative of the above function
    fprime= lambda x: -eps/x-dT
    fprime2= lambda x:eps/(x**2)

    # Compute rho:
    tol=np.power(10.,-13)
    err = 10.
    i = 0
    r0 = np.sum(np.exp(u/eps))*dX # "dummy" rho

    # rho=scipy.optimize.newton(func, r0, fprime=fprime, args=(), tol=tol, maxiter=1000, fprime2=fprime2)
    r = r0
    while err>tol:
        i+=1
        newr=r-invd_func(r)*func(r)
        err=np.abs(newr-r)
        r=newr
        if i>1000:
            break
    rho = r
    u = -parameters['C']*dT*rho + A
    return u, rho

def build_and_save_HJ(u, rho, parameters, pre_init_values, path):
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
    sum_f = rho
    computed_mean = np.divide(np.sum(X*f, axis = 1), sum_f)*dX
    figure = plt.figure()
    plt.suptitle(par_str, y = 1.1)
    grid = plt.GridSpec(2,5, wspace = 0.9, hspace = 0.5)
    fig1 = plt.subplot(grid[:,:-2])
    fig1.imshow(u.transpose()[::-1],cmap=plt.cm.jet, aspect = 'auto', extent = (0,parameters['T_max'], X_min, X_max), 
            vmin = np.min(u)*dX, vmax = np.max(u)*dX)
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

parameters = dict(T_max = 50, # maximal time 
                  dT = 1e-3, # Discretization time 
                  sigma0 = 1,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  sigma = 0.01,
                  tau = 0.35,  # transfer rate
                  X_min = -4, #length of the numerical interval of traits (for PDE!)
                  X_max=4,
                  dX = 0.01, #discretization of the space of trait6
                  eps = 1e-6,
                  delta=1e-4
                  )

for tau_i in np.arange(0.,1.,0.05):
    for eps_i in np.geomspace(1e-6, 1, num = 7):
        parameters['tau'] = tau_i
        parameters['eps'] = eps_i
        print("Currently, tau="+str(tau_i)+" and eps="+str(eps_i),flush = True)
        pre_init_values = Pre_Initialization_HJ(parameters) 
        U, Rho, nT = pre_init_values['U'], pre_init_values['Rho'], pre_init_values['nT']
        for i in range(nT-1):
            U[i+1], Rho[i+1] = Next_Generation_AP(U[i], Rho[i], parameters, pre_init_values)
            if i%1000==0: 
                print(str(i)+"-th iteration", flush = True)
        build_and_save_HJ(U, Rho, parameters, pre_init_values, "/scratch/gene/horizontal_gene_transfer/Mesocentre/Figures/HJ/")
