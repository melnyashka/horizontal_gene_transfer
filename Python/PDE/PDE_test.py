parameters = dict(T_max = 500, # maximal time 
                  dT = 0.01, # Discretization time 
                  K = 1000, # Maximal capacity of the system
                  N0 = 1000,    # Initial number of population
                  sigma0 = 0.1,  #Initial standard variation of the population
                  x_mean0 = 1.,
                  C = 0.5,    # competition
#                 p = 0.03,      # Probability of mutation
                  b_r = 1,     # birth rate
                  d_r = 1,      # death rate
                  d_e = 2,   #exponetial power
                  beta = 0, 
                  mu = 1,
                  sigma = 0.1,
                  tau = 0.8,  # transfer rate
                  L = 3, #length of the numerical interval of traits (for PDE!)
                  dX = 0.01, #discretization of the space of traits
                  eps = 1
                  )

L, dX, T_max, dT = parameters['L'], parameters['dX'], parameters['T_max'], parameters['dT']
X_min, X_max= -L, L
nT = int(T_max/dT)
nX = int((X_max-X_min)/dX) # number of traits
XT = np.empty([nT, nX]) # initial density 

pre_values = compute_things(parameters)                
XT[0] = pre_values['init_density']

for i in range(nT-1):
    XT[i+1] = Next_Generation_PDE(XT[i],parameters, pre_values)
    if i%10==0:
        print('T= '+str(i*parameters['dT']))





#PLOT !!!!!!!!!!!!!!!!!

import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import cm

figure = plt.figure()
im = imshow(XT,cmap=cm.coolwarm,aspect='auto')
colorbar(im)

