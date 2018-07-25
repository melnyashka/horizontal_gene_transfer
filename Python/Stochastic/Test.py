#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 16:10:45 2018

@author: samuelnordmann
"""
import matplotlib.pyplot as plt

for i in range(int(T_max/dT-1)):
    X[i+1]=Next_Generation(X[i])
    
#X2=np.asarray(X)
#imgplot = plt.imshow(X2)