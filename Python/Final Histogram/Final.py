#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 21:56:40 2018

@author: samuelnordmann
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import gc




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################  Parameters  ##############################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



parameters = dict(T_max = 100, # maximal time 
                  dT = 0.01, # Discretization time 
                  sigma0 = 1,  #Initial standard variation of the population
                  x_mean0 = 0.,
                  rho0=2.,
                  K=300,
                  C = 0.5,    # competition
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  sigma = 1,
                  tau = 0.4,  # transfer rate
                  X_min = -1, #length of the numerical interval of traits (for PDE!)
                  X_max=3,
                  dX = 0.1, #discretization of the space of traits
                  eps = 0.01,
                  delta=0.01
                  )

parameters_S =parameters.copy()
parameters_PDE =parameters.copy()

parameters_HJ =parameters.copy()
parameters_HJ['T_max']=int(parameters_HJ['T_max']*parameters_HJ['eps'])
parameters_HJ['dT']=parameters_HJ['dT']*parameters_HJ['eps']
#parameters_HJ['sigma']=parameters_HJ['sigma']/parameters_HJ['eps']





#Length and speed of video setting
frameNumber = 100
n=int(parameters['T_max']/parameters['dT'])
c=int(n/frameNumber)












##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################     FUNCTIONS        ######################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
##############################    FUNCTIONS  STOCH        ################################
##########################################################################################
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
        beta_transfer = np.divide(1,horizontal_transfer(x, tau = parameters['tau']))
        times = np.array([np.random.exponential(beta_birth),np.random.exponential(beta_death), np.random.exponential(beta_transfer)])
        b_mat = (times < parameters['dT'])

        return np.sort(np.concatenate((x[np.logical_not(np.logical_or(b_mat[1],b_mat[2]))],
                                       np.random.normal(loc=x[b_mat[0]], scale=sigma, size=None),
                                       np.vectorize(lambda i: np.random.choice(x[(i+1):]),otypes=[np.float64])(np.arange(n_tot)[b_mat[2]][:-1]))))


##########################################################################################
##############################    FUNCTIONS  PDE        ######################################
##########################################################################################

def Pre_Initialization_PDE(parameters_PDE):
    dX, T_max, dT = parameters_PDE['dX'], parameters_PDE['T_max'], parameters_PDE['dT']
    X_min, X_max= parameters_PDE['X_min'], parameters_PDE['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
       
    f0 = np.exp(-np.power(np.absolute(X-parameters_PDE['x_mean0']),2)/parameters_PDE['sigma0'])# initial density 
    rho0=(parameters_PDE['b_r'])/(parameters_PDE['C'])
    f0=f0/(np.sum(f0)*parameters_PDE['dX'])*rho0

    f = np.empty([nT, nX])        #densities for all times
    f[0]=f0
    
    # Computing constant death and mutation kernels
    Death = parameters_PDE['d_r']*np.power(np.absolute(X),parameters_PDE['d_e'])
    Mutation_kernel=np.exp(-(X/(parameters_PDE['sigma']))**2*1/2)*parameters_PDE['dX']
    Mutation_kernel=parameters_PDE['b_r']*Mutation_kernel/np.sum(Mutation_kernel)
    
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

def Next_Generation_PDE(f,parameters_PDE, pre_init_values):
    dX, T_max, dT = parameters_PDE['dX'], parameters_PDE['T_max'], parameters_PDE['dT']
    X_min, X_max= parameters_PDE['X_min'], parameters_PDE['X_max']
    Death, Mutation_kernel = pre_init_values['Death'], pre_init_values['Mutation_kernel']

    X_min_new=int(-X_min/dX)          # Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-2*X_min+X_max)/dX)
    I=range(pre_init_values['nX'])    # all the indexes
    
    rho = np.sum(f)*parameters_PDE['dX'] # I don't like the idea of comparing float with an integer, but okay
    
    death_term = Death + rho*parameters_PDE['C']
    birth_part = np.convolve(Mutation_kernel,f)[X_min_new:X_max_new]
    if rho>np.power(10.,-7):

        X_large=np.arange(2*X_min,2*X_max,dX)
        transfer_kernel=  parameters_PDE['tau']*(np.tanh(np.divide(X_large,parameters_PDE['delta']))-np.tanh(-np.divide(X_large,parameters_PDE['delta'])))*dX

        X_min_new=int(-2*X_min/dX)#Bounds for the new indexes that must be kept after the convolution
        X_max_new=int((-3*X_min+X_max)/dX)
        T = np.convolve(transfer_kernel,f/rho)[X_min_new:X_max_new]# that's the transfer term    
    else :
        T = 0
    new_f = f + dT*(np.multiply(f,(-death_term + T))+birth_part)
    return new_f



##########################################################################################
##############################    FUNCTIONS  HJ        ######################################
##########################################################################################


def Pre_Initialization_HJ(parameters_HJ):
    dX, T_max, dT = parameters_HJ['dX'], parameters_HJ['T_max'], parameters_HJ['dT']
    X_min, X_max= parameters_HJ['X_min'], parameters_HJ['X_max']
    
    nT = int(T_max/dT)            #number of times
    T = np.arange(0,T_max,dT)     #space of time

    nX = int((X_max-X_min)/dX)    #number of traits
    X = np.arange(X_min,X_max,dX) #space of traits
    var_0=parameters_HJ['sigma0']**2 #*parameters_HJ['eps']   
    # u0 = np.exp(-np.power(np.absolute(X-parameters_HJ['x_mean0']),2)/((1+np.power(X,2))/2*(parameters_HJ['sigma0']*parameters_HJ['eps'])**2) # initial density 
    u0 = -((np.abs(X-parameters_HJ['x_mean0'])<=1)*(np.power(X-parameters_HJ['x_mean0'],2)/(2*var_0))
           +(np.abs(X-parameters_HJ['x_mean0'])>1)*(np.abs(X-parameters_HJ['x_mean0'])/var_0-1/(2*var_0)))
    #u0= -np.divide(np.power(X-parameters_HJ['x_mean0'],2),1+2*np.power(X-parameters_HJ['x_mean0'],2))/(2*parameters_HJ['sigma0']**2)
    # Helene's trick:
    # u=(abs(x)<=1).*x.^2/2+(abs(x)>1).*(abs(x)-1/2) # first we initialize u
    rho0=np.sum(np.exp(u0/parameters_HJ['eps']))*dX  # then we initialize rho
    u0=u0+parameters_HJ['eps']*np.log(parameters_HJ['rho0']/rho0)
    f0 = np.exp(u0/parameters_HJ['eps']) # then we initialize f

    f = np.empty([nT, nX])        #densities for all times
    U = np.empty([nT, nX])        #
    U[0]=u0 # matrix for u
    f[0]=f0 # matrix for f
    Rho=np.empty([nT])
    Rho[0] = parameters_HJ['rho0']
    
    # Computing constant death and mutation kernels
    Death = parameters_HJ['d_r']*np.power(np.absolute(X),parameters_HJ['d_e'])
    
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


def Next_Generation_AP(u, rho, parameters_HJ, pre_init_values):
    ### That's our draft for the AP scheme
    dX, T_max, dT = parameters_HJ['dX'], parameters_HJ['T_max'], parameters_HJ['dT']
    X_min, X_max, eps = parameters_HJ['X_min'], parameters_HJ['X_max'],  parameters_HJ['eps']
    Death, X = pre_init_values['Death'], pre_init_values['X']

    # now we define a grid to evaluate an integral, cutting all the unnecessary values
    ymax = 8*parameters_HJ['sigma']
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
    transfer_kernel=  parameters_HJ['tau']*(np.tanh(np.divide(X_large,parameters_HJ['delta']))-np.tanh(-np.divide(X_large,parameters_HJ['delta'])))*dX
    # transfer_kernel = parameters_HJ['tau']*(np.heaviside(-np.divide(X_large,parameters_HJ['delta']),1)-np.heaviside(np.divide(X_large,parameters_HJ['delta']),1))*dX
    
    X_min_new=int(-2*X_min/dX)#Bounds for the new indexes that must be kept after the convolution
    X_max_new=int((-3*X_min+X_max)/dX)
    T = np.convolve(transfer_kernel,fsurrhou)[X_min_new:X_max_new]# that's the transfer term


    mut_kern = np.exp(-np.power(Yx,2)/(2*parameters_HJ['sigma']**2))
    cste=np.sum(mut_kern[:,0])*dy
    m_sq = mut_kern/cste*np.ones([len(X)])

    H = parameters_HJ['b_r']*np.exp(flux - np.power(np.ones([len(X)])*Yx,2)/2)/cste #*m_sq #*m(x) add a gaussian mutation kernel!
    H = np.sum(H,axis=0)*dy # that's the birth term

    A = u + dT*(-pre_init_values['Death'] + H + T) # that's the birth-death-transfer term
    max_A = np.max(A)
    E = eps*np.log(dX) + max_A + eps*np.log(np.sum(np.exp((A-max_A)/eps)))
    E2=dX*np.sum(np.exp((A-max_A)/eps))
    
    func = lambda x: E - eps*np.log(x)-parameters_HJ['C']*dT*x#the unknown rho is the root of this function
    invd_func = lambda x: -x/(eps+dT*x)#inverse of the derivative of the above function
    #fprime= lambda x: -eps/x-dT
    #fprime2= lambda x:eps/(x**2)

    func2= lambda x: x*np.exp((dT*parameters_HJ['C']*x-max_A)/parameters_HJ['eps'])-E2
    invd_func2= lambda x: parameters_HJ['eps']/(1+parameters_HJ['C']*dT*x)*np.exp(-(dT*parameters_HJ['C']*x-max_A)/parameters_HJ['eps'])
    # Compute rho:
    tol=np.power(10.,-7)
    tol0=np.power(10.,-10)
    threshold=np.power(10.,-10)
    err = 10
    i = 0
    if rho>threshold:      
        try:
            #r0 = np.exp(max_u/eps)*np.sum(np.exp((u-max_u)/eps))*dX # "dummy" rho
            r0 = np.sum(np.exp(u/eps))*dX
        except:
            r0=10
            print('erreur dummy rho')
    else:
        r0=1
        #print('BELOW THRESHOLD. rho = '+str(rho))

    # rho=scipy.optimize.newton(func, r0, fprime=fprime, args=(), tol=tol, maxiter=1000, fprime2=fprime2)
    r = r0
    while err>tol and r>=tol0 and i<1000:
        i+=1
        if r>threshold:
            newr=r-invd_func(r)*func(r)
        else:
            newr= r-invd_func2(r)*func2(r)
            #print('methode 2')
        err=np.abs(newr-r)
        #print('for i= '+str(i)+' we have r= '+str(r)+' and newr = '+str(newr))
        r=newr

    rho = max(r,0)
    u = -parameters_HJ['C']*dT*rho + A
    #u= A
    #u=u-np.maximum(np.max(u),0)
    return u, rho
















##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################  COMPUTATIONS STOCH  ######################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



#Initialization

X0 = np.random.normal(parameters_S['x_mean0'], parameters_S['sigma0'], int(parameters_S['rho0']*parameters_S['K'])) # Initial population
#X = [None]*int(parameters_S['T_max']/parameters_S['dT'])  # history of all populations up to time T_max
X = np.sort(X0)
Size_S=[X0.size/parameters_S['K']]
Mean_S=[np.mean(X0)]

#Initialization of histogram
X_max=parameters_S['X_max']
X_min=parameters_S['X_min']
bins=int((X_max-X_min)/parameters_HJ['dX'])
rg=(X_min,X_max)

Abs_Stoch=np.histogram(X0,bins,rg,True)[1]
U_S_0=np.histogram(X0,bins,rg,True)[0]*parameters_S['rho0']

U_S=np.zeros((int(parameters_S['T_max']/parameters_S['dT']),bins))
U_S[0]=U_S_0


nbr_iteration= int(parameters_S['T_max']/parameters_S['dT'])
for i in range(nbr_iteration-1):
    if i%500==0:
        print('STOCH: '+str(i)+'th iteration out of '+str(nbr_iteration)+'.')
    X=Next_Generation(X, parameters_S)
    U_S[i+1]=np.histogram(X,bins,rg,False)[0]/parameters_S['K']
    Size_S=Size_S+[X.size/parameters_S['K']]
    Mean_S=Mean_S+[np.mean(X)]




gc.collect()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################  COMPUTATIONS PDE  ######################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################




pre_init_values = Pre_Initialization_PDE(parameters_PDE) 
U_P, nT = pre_init_values['f'], pre_init_values['nT']

for i in range(nT-1):
    U_P[i+1] = Next_Generation_PDE(U_P[i],parameters_PDE, pre_init_values)
    if i%500==0:
        print('PDE: T= '+str(i*parameters_PDE['dT']), flush = True)




gc.collect()



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################  COMPUTATIONS HJ  ######################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################




pre_init_values = Pre_Initialization_HJ(parameters_HJ) 
U, Rho, nT = pre_init_values['U'], pre_init_values['Rho'], pre_init_values['nT']
U_H=U
for i in range(nT-1):
    U[i+1], Rho[i+1] = Next_Generation_AP(U[i], Rho[i], parameters_HJ, pre_init_values)
    if i%500==0: 
        print('HJ : '+str(i)+"-th iteration, over "+str(nT))
    A=np.max(U[i+1])
    U_H[i+1]=np.exp(A/parameters_HJ['eps'])*np.exp((U[i+1]-A)/parameters_HJ['eps'])

#U_H=np.exp(U_H/parameters_HJ['eps'])




gc.collect()



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##############################  ANIMATION ################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################





# we generate some synthetic data.
#we will make a list of X, Y lists.
frames1 = []
frames2 = []
frames3 = []


for i in range(0,frameNumber-1):
    xData = np.arange(parameters['X_min'],parameters['X_max'],parameters['dX'])
    yData1 = U_S[c*i]
    yData2= U_P[c*i]
    yData3= U_H[c*i]
    frames1.append( [xData, yData1] )
    frames2.append( [xData, yData2] )    
    frames3.append( [xData, yData3] )    
    
    
# now we put them together to make a movie! let's set up the movie writer
framerate = 60 # 24 frames per second
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='HT transfer Simulations', artist='Samuel Nordmann', comment='')
writer = FFMpegWriter(fps=framerate, metadata=metadata)


# figure setup

Abs_m=parameters_HJ['X_min']
Abs_M=parameters_HJ['X_max']
Ord_m=-0.1
Ord_M=1.5

fig = plt.figure()

fig1=plt.subplot(3,1,1)
fig1.set_xlim(Abs_m, Abs_M)
fig1.set_ylim(Ord_m, Ord_M)


fig2=plt.subplot(312, sharex=fig1)
fig2.set_xlim(Abs_m, Abs_M)
fig2.set_ylim(Ord_m, 4)

fig3=plt.subplot(313, sharex=fig1)
fig3.set_xlim(Abs_m, Abs_M)
fig3.set_ylim(Ord_m, 10)




firstFrameX1, firstFrameY1 = frames1[0]
firstFrameX2, firstFrameY2 = frames2[0]
firstFrameX3, firstFrameY3 = frames3[0]
l1, = fig1.plot(firstFrameX1, firstFrameY1, '-',label='$Stoch$')
l3, = fig3.plot(firstFrameX3, firstFrameY3, '-',color='green',label='$HJ$')
l2, = fig2.plot(firstFrameX2, firstFrameY2, '-',color='red',label='$PDE$')


plt.setp(fig1.get_xticklabels(), visible=False)
plt.setp(fig2.get_xticklabels(), visible=False)


l1.set_xdata(xData)
fig.legend(loc=0,fontsize='large')

plt.xlabel( "$x$" )
plt.xlim(parameters_HJ['X_min'],parameters_HJ['X_max'])
#plt.title()


# let's write to the file!
with writer.saving(fig, 'test2'+'.mp4', dpi=100):
    for i in range(frameNumber-1):
        if i%25==0:
            print('plot_i='+str(i)+' over '+str(frameNumber))
        x1, y1 = frames1[i]
        l1.set_data( x1, y1)
        
        x2, y2 = frames2[i]
        l2.set_data( x2, y2)
        
        x3, y3 = frames3[i]
        l3.set_data( x3, y3)
        
        writer.grab_frame()
