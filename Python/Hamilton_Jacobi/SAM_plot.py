#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 11:14:57 2018

@author: samuelnordmann
"""

from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


figure = plt.figure()
#im = imshow(u.transpose(),cmap=cm.coolwarm)
im = imshow(u.transpose()[::-1],cmap=cm.coolwarm,aspect='auto',extent=(0,parameters['T_max'],parameters['X_min'],parameters['X_max']),vmin=-30)
colorbar(im)
par_str = '' # create a string of parameters to pass into plots
for k, v in parameters.items():
    if k == 'N0' or k == 'd_r': 
        smth = ",\n" 
    else: 
        smth = ", "
    par_str += k + "=" + str(v) + smth

plt.xlabel('time')
plt.ylabel('trait');


plt.title(par_str)
plt.tight_layout()



current_time = datetime.now().time()
figure.savefig(str("Figures/APscheme/plot_" + str(current_time)[0:8].replace(':','_')+".pdf")) # Possibly different delimeter on Linux and Windows!

