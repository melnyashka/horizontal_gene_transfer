#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:37:24 2018

@author: samuelnordmann
"""
import matplotlib.pyplot as plt
import numpy as np

#import plotly.plotly as py

for t in np.arange(200,500,1):
    plt.clf()
    figure = plt.figure()
    Ord_new=Ord[Abs==t]
    #plt.plot(Ord_new)
    
    plt.hist(Ord_new,bins=20)
    plt.title("t="+str(t)+', tau='+str(parameters['tau']))
    plt.xlabel("Trait")
    plt.ylabel("Number of individuals")
    figure.savefig("Figures/"+"histogram_"+"t="+str(t)+".pdf", bbox_inches ='tight')

#fig = plt.gcf()