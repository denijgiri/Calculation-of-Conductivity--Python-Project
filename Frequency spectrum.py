# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:53:26 2024

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt

doped_td=np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19073/2023-07-25T15-36-11.677721-100avg-sample-X_25.000 mm-Y_5.000 mm.txt")

fourie_doped = np.fft.fft(doped_td[:,1])
n= doped_td[:,0].size
timestep=5*10**-14

freq = np.fft.fftfreq(n,d=timestep)
pos_freq = freq>0
freq = freq[pos_freq]
fourie_doped = fourie_doped[pos_freq]

db_doped = 20*np.log10(np.abs(fourie_doped))

plt.figure()
plt.plot(freq/10**12,db_doped, color = 'Blue', label ='Doped')
plt.xlim(0,5)
plt.title('Spectrum of the sample')
plt.xlabel('Frequency (THz)')
plt.ylabel('Amplitude[a.u]')
plt.legend()
plt.show()