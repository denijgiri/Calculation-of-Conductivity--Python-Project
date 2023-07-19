# -*- coding: utf-8 -*-
"""
Created on Wed May 31 07:38:31 2023

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt

doped_td=np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-44-05.321852-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
sam_td=np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-20-30.305680-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
ref_td=np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-20-10.178856-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt")




fourie_doped = np.fft.fft(doped_td[:,1])
n= doped_td[:,0].size
timestep=5*10**-14

freq = np.fft.fftfreq(n,d=timestep)
pos_freq = freq>0
freq = freq[pos_freq]
fourie_doped = fourie_doped[pos_freq]

fourie_sam = np.fft.fft(sam_td[:,1])
n = sam_td[:,0].size
timestep=5*10**-14

freq = np.fft.fftfreq(n,d=timestep)
pos_freq = freq>0
freq = freq[pos_freq]
fourie_sam = fourie_sam[pos_freq]

fourie_ref = np.fft.fft(ref_td[:,1])
n = ref_td[:,0].size
timestep=5*10**-14

freq = np.fft.fftfreq(n,d=timestep)
pos_freq = freq>0
freq = freq[pos_freq]
fourie_ref = fourie_ref[pos_freq]



db_doped = 20*np.log10(np.abs(fourie_doped))
db_sam = 20*np.log10(np.abs(fourie_sam))
db_ref = 20*np.log10(np.abs(fourie_ref))

plt.figure()
plt.plot(freq/10**12,db_doped, color = 'Blue', label ='Doped')
plt.plot(freq/10**12,db_sam, color = 'Red', label = 'Substrate')
plt.plot(freq/10**12,db_ref, color = 'Orange', label = 'Reference')
plt.xlim(0,5)
plt.title('Spectrum of the sample')
plt.xlabel('Frequency (THz)')
plt.ylabel('Amplitude[a.u]')
plt.legend()
plt.show()