# import stochastic_continuous # Check if it works on your machine! 
parameters = dict(T_max = 1000, # maximal time 
                dT = 0.1, # Discretization time 
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
                tau = 0.28 # transfer rate
                )

# Change some parameters if needed!
# parameters['tau'] = 0.17
# idea: create a grid of parameters we want to check, and then run the experiments inside the loop! 

X0 = np.random.normal(parameters['x_mean0'], parameters['sigma0'], parameters['N0']) # Initial population
#X = [None]*int(parameters['T_max']/parameters['dT'])  # history of all populations up to time T_max
X = np.sort(X0)
l= np.empty(X0.size, dtype=np.object)
for i in range(X0.size):
    l[i]=[]
    
Abs=[]#Save all the individuals
Ord=[]

for i in range(int(parameters['T_max']/parameters['dT']-1)):
    if i%5000==0:
        print('T= '+str(i*parameters['dT']))
    for x in X:
        Abs.append(i*parameters['dT'])
        Ord.append(x)
    X,l=Next_Generation_lineage(X,l,i*parameters['dT'], parameters)
    
    
    
    
    
    
    
    
    
#Plot
figure = plt.figure()

plt.hist2d(Abs,Ord,bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)



n_plot=15
l_plot=np.random.choice(l,n_plot,replace=False)



for e in l_plot:
    Abs_plot=[]
    Ord_plot=[]
    for j in e:
        Abs_plot.append(j[0])
        Ord_plot.append(j[1])
    plt.plot(Abs_plot,Ord_plot)




    
    

par_str = '' # create a string of parameters to pass into plots
for k, v in parameters.items():
    if k == 'N0' or k == 'b_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth


plt.xlabel('Time. Sample of size '+str(n_plot))
plt.ylabel('trait');
plt.title(par_str)
#plt.show()
current_time = datetime.now().time()
figure.savefig(str("Figures_Lineage/plot_" + str(current_time)[0:8]+".pdf")) # Possibly different delimeter on Linux and Windows!


gc.collect()
