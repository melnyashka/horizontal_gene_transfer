#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 16:26:49 2018

@author: samuelnordmann
"""

import matplotlib.pyplot as plt
from numpy import arange

for t in range(nT-1):
    N[t+1]=next_time(N[t])



from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt


# the function that I'm going to plot
levellines=np.concatenate((arange(10,100,50),arange(100,1000,250)))
im = imshow(N,cmap=cm.coolwarm) # drawing the function
# adding the Contour lines with labels
#cset = contour(N,levellines,linewidths=2,cmap=cm.Set2)
#clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
colorbar(im) # adding the colobar on the right
# latex fashion title
title('Population')
show()

