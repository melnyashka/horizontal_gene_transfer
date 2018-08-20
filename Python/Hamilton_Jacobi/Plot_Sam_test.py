#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:02:38 2018

@author: samuelnordmann
"""

for t in np.arange(13,25,0.5):
    plt.clf()
    figure = plt.figure()
    #t=19.5
    X_petit=X[int((-X_min-0.2)/dX):int((-X_min+1.75)/dX)]
    U_petit= U[int(t/parameters['dT']),int((-X_min-0.2)/parameters['dX']):int((-X_min+1.75)/parameters['dX'])]
    
    string= 't='+str(t)+', eps='+str(parameters['eps'])+', tau='+str(parameters['tau'])
    
    plt.plot(X_petit,U_petit)
    plt.title(string)
    plt.xlabel('trait')
    plt.ylabel('u')
    figure.savefig("Figures/APscheme/"+"histogram_"+"t="+str(t)+".pdf", bbox_inches ='tight')

