# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 16:36:06 2023

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import signal
from Measured import average

def drude_model(freqs,mobility,N,mass_eff,e):
    sigma_dc = N * mobility *e
    tau = mobility * mass_eff /(e*10**4)
    omega = 2 * np.pi *freqs
    sigma_cc = (sigma_dc /(1 - (1j*omega*tau)))
    
    plt.figure()
    plt.plot(freqs,sigma_cc.imag,color='red', label = 'Imaginary part(Drude Model)')
    plt.plot(freqs, average.imag,color='blue', label='Imaginary part')
    plt.xlabel('Frequency(THz)')
    plt.ylabel('Conductivity[Ω-1cm-1]')
    plt.title('Conductivity of film')
    plt.figure()
    plt.plot(freqs,sigma_cc.real,color='green',label = 'Real Part(Drude Model)')
    plt.plot(freqs, average.real,color='orange', label='Real part')
    plt.xlabel('Frequency[Hz]')
    plt.ylabel('Conductivity[Ω-1cm-1]')
    plt.title('Conductivity of film')
    plt.ylim((-50, 680))
    plt.legend()
    plt.show()
    
freqs = np.linspace(0.5*10**12,3*10**12,1000)
drude_model(freqs,2500,1.3*10**18,0.067*9.1*10**-31,1.6*10**-19)





