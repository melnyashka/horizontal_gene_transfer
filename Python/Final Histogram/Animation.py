#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:34:42 2018

@author: samuelnordmann
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


# we generate some synthetic data.
#we will make a list of X, Y lists.
frames1 = []
frames2 = []
frames3 = []
frameNumber = 400
c=25 

for i in range(0,frameNumber):
    xData = np.arange(parameters_HJ['X_min'],parameters_HJ['X_max'],parameters_HJ['dX'])
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
Ord_M=2.5

fig = plt.figure()

fig1=plt.subplot(3,1,1)
fig1.set_xlim(Abs_m, Abs_M)
fig1.set_ylim(Ord_m, Ord_M)


fig2=plt.subplot(3,1,2)
fig2.set_xlim(Abs_m, Abs_M)
fig2.set_ylim(Ord_m, 5)

fig3=plt.subplot(3,1,3)
fig3.set_xlim(Abs_m, Abs_M)
fig3.set_ylim(Ord_m, Ord_M)




firstFrameX1, firstFrameY1 = frames1[0]
firstFrameX2, firstFrameY2 = frames2[0]
firstFrameX3, firstFrameY3 = frames3[0]
l1, = fig1.plot(firstFrameX1, firstFrameY1, '-',label='$Stoch$')
l2, = fig2.plot(firstFrameX2, firstFrameY2, '-',color='red',label='$PDE$')
l3, = fig3.plot(firstFrameX3, firstFrameY3, '-',color='green',label='$HJ$')




l1.set_xdata(xData)
fig.legend(loc=0,fontsize='large')

plt.xlabel( "$x$" )
plt.xlim(parameters_HJ['X_min'],parameters_HJ['X_max'])
#plt.title()


# let's write to the file!
with writer.saving(fig, 'test'+'.mp4', dpi=100):
    for i in range(frameNumber):
        if i%100==0:
            print('plot_i='+str(i)+' over '+str(frameNumber))
        x1, y1 = frames1[i]
        l1.set_data( x1, y1)
        
        x2, y2 = frames2[i]
        l2.set_data( x2, y2)
        
        x3, y3 = frames3[i]
        l3.set_data( x1, y1)
        
        writer.grab_frame()
