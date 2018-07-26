#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 16:10:45 2018

@author: samuelnordmann
"""
import matplotlib.pyplot as plt

Abs=[]   
Ord=[]

for i in range(int(T_max/dT-1)):
    X[i+1]=Next_Generation(X[i])
    for x in X[i]:
        Abs.extend([i*dT])
        Ord.extend([x])




#plt.scatter(Abs,Ord,s=0.01,alpha=0.2)
plt.hist2d(Abs,Ord,bins=K,cmap=plt.cm.Blues)
plt.colorbar()

plt.xlabel('time')
plt.ylabel('trait');

titre='T_max = '+str(T_max)+' ; K = '+str(K)+' ; N0 = '+str(N0)+' ; sigma0 = '+str(sigma0)+'\nC = '+str(C)+' ; p = '+str(p)+' ; b_r = '+str(b_r)+'\nd_r = '+str(d_r)+' ; beta = '+str(beta)+' ; mu = '+str(mu)+' ; sigma = '+str(sigma) +' ; tau = '+str(tau)          
plt.title(titre)



plt.show()





