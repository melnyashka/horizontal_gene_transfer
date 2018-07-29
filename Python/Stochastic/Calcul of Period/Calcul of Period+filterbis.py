#!/usr/bin/env python3
# -*- coding: utf-8 -*-






##WHOLE GRAPH
t_plot=int(1000/parameters['dT'])
T_plot=int(4000/parameters['dT'])

t_calculus=int(1000/parameters['dT'])
T_calculus=int(parameters['T_max']/parameters['dT'])


#mask=np.logical_and(Abs>t_plot,Abs<T_plot)
#
#plt.subplot(3, 2, 1)
#plt.hist2d(Abs[mask],Ord[mask],bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)
#par_str = '' # create a string of parameters to pass into plots
#for k, v in parameters.items():
#    if k == 'N0' or k == 'b_r': 
#        smth = ",\n" 
#    else: 
#        smth = ", "
#    par_str += k + "=" + str(v) + smth
##figure = plt.figure()
##plt.hist2d(Abs,Ord,bins=3/2*parameters['K'],cmap=plt.cm.bone_r,alpha=1,cmax=2*parameters['K'],cmin=parameters['K']/100)
##plt.colorbar()
##plt.xlabel('time')
##plt.ylabel('trait');
#plt.title(par_str)







#FILTERING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from scipy import signal

freq_ech=parameters['T_max']/parameters['dT']
freq_expected=parameters['T_max']/100
delta=0.5

low = freq_expected*(1-delta) / ( freq_ech)
high = freq_expected*(1+delta) / ( freq_ech)

typefilter='bandpass'

if typefilter=='bandpass':
    title_filter='Bandpass ('+str(low)+','+str(high)+').'
    band=[low,high]
else:
    title_filter='Lowpass '+str(low)+'.'
    band=high

#critical_f=0.0001
b, a = signal.butter(3, [low,high], btype=typefilter, analog=False)














#TOTAL SIZE OF THE POPULATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Pop_tot=Ord_tot[t_calculus:T_calculus]
Pop_tot=Pop_tot-Pop_tot.mean()
Pop_tot=signal.filtfilt(b, a, Pop_tot)



plt.subplot(2, 2, 1)
plt.plot(Abs_tot[t_plot+1:T_plot],Pop_tot[t_plot:T_plot])
plt.title('total Population. Filter = '+title_filter)


Fourier=np.fft.rfft(Pop_tot)
freq=np.fft.rfftfreq(Pop_tot.shape[-1])

plt.subplot(2, 2, 2)
plt.plot(freq,Fourier)
plt.title('Fourier transform of the total population')








#MEAN TRAIT OF THE POPULATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Mean_trait=Ord_mean[t_calculus:T_calculus]
Mean_trait=signal.filtfilt(b, a, Mean_trait)
Mean_trait=Mean_trait-Mean_trait.mean()

plt.subplot(2, 2, 3)
plt.plot(Abs_mean[t_plot+1:T_plot],Mean_trait[t_plot:T_plot])
plt.title('mean trait')


Fourier_trait=np.fft.rfft(Mean_trait)
freq_trait=np.fft.rfftfreq(Mean_trait.shape[-1])

plt.subplot(2, 2, 4)
plt.plot(freq_trait,Fourier_trait)
plt.title('Fourier transform of the mean trait')







#Count the fundamental





plt.show()





