# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:27:49 2023

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["axes.grid"] = True

def phase_correction(freq,phase):
    freq_range=np.array([1,2])*10**12
    freq_slice_idx = (freq>=freq_range[0]) * (freq <= freq_range[1])
    
    z=np.polyfit(freq[freq_slice_idx],phase[freq_slice_idx],1)
    
    # phase_corrected = phase - z[1]
    phase_corrected = freq*z[0]
    
    freq_cp = freq.copy()/10**12
    
    plt.figure()
    #plt.plot(freq_cp,phase,label="phase b4 correction")
    #plt.plot(freq_cp,phase_corrected,label="phase after correction")
    #plt.plot(freq_cp[freq_slice_idx],phase[freq_slice_idx],label="fit range")
    plt.plot(freq_cp, freq*z[0],label="fit")
    plt.ylabel("Unwrapped phase(rad)")
    plt.xlabel("Frequency(Thz")
    plt.legend()

    return phase_corrected

def fourie_transform(data_td, ret_phase=True):
    fourie=np.fft.fft(data_td[:,1])
    
    n=data_substrate_sam[:,1].size
    timestep=5*10**-14
    freq=np.fft.fftfreq(n,d=timestep)
    
    #freq_range = (freq>=0.0*10**12)*(freq<=2.2*10**12)
    freq_range = (freq>=0*10**12)*(freq<=2.5*10**12)
    
    freq=freq[freq_range]
    fourie = fourie[freq_range]
    
    if not ret_phase:
        return fourie, freq
    
    phase = np.angle(fourie)
    delta = np.unwrap(phase)
    
    plt.figure()
    plt.plot(freq/10**12,phase,label="wrapped phase")
    plt.ylabel("phase(rad)")
    plt.xlabel("Frequency(Thz")
    plt.legend()
    
    return delta, freq



data_substrate_sam = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-18-31.067459-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
data_substrate_ref = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-18-11.052862-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt")

data_sam = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-42-06.856860-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
data_ref = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-41-46.977074-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt")



n=data_substrate_sam[:,1].size
timestep=5*10**-14
freq_full_axis=np.fft.fftfreq(n,d=timestep)
freq_range = (freq_full_axis>=0*10**12)*(freq_full_axis<=2.5*10**12)


fft_substrate_sam = np.fft.fft(data_substrate_sam[:,1])
fft_substrate_ref = np.fft.fft(data_substrate_ref[:,1])


fft_sam = np.fft.fft(data_sam[:,1])
fft_ref = np.fft.fft(data_ref[:,1])

phi_ref_sub = np.unwrap(np.angle(fft_substrate_ref))
phi_ref_sam = np.unwrap(np.angle(fft_substrate_sam))

phi_sam = np.unwrap(np.angle(fft_sam))
phi_ref = np.unwrap(np.angle(fft_ref))



phi_ref_sub_corrected = phase_correction(freq,phi_ref_sub) 
phi_ref_sam_corrected = phase_correction(freq,phi_ref_sam)
phi_sub_corrected = phase_correction(freq,phi_ref)
phi_sam_corrected = phase_correction(freq,phi_sam)


#delta_sam_substrate_offset_corrected = phase_correction(freq, delta_substrate_sam)
#delta_ref_offset_corrected = phase_correction(freq,delta_ref)
#delta_sam_offset_corrected = phase_correction(freq,delta_sam)

r_sam = np.abs(fft_sam)
r_sub = np.abs(fft_ref)
r_ref_sub = np.abs(fft_substrate_ref)
r_ref_sam = np.abs(fft_substrate_sam)

ref_sub_fd = r_ref_sub * np.exp(1j*phi_ref_sub_corrected)
sam_sub_fd = r_ref_sam * np.exp(1j*phi_sub_corrected)

ref_sam_fd = r_sam * np.exp(1j*phi_ref_sam_corrected)
sub_sam_fd = r_sub * np.exp(1j*phi_sam_corrected)


T_substrate_sam = sam_sub_fd / ref_sub_fd
T_sam = sub_sam_fd / ref_sam_fd


#T_substrate_sam = fft_substrate_sam / fft_substrate_ref
#T_sam = fft_sam / fft_ref

T_substrate_sam = T_substrate_sam[freq_range]
T_sam = T_sam[freq_range]


delta_substrate_ref, freq = fourie_transform(data_substrate_ref)
delta_substrate_sam, freq = fourie_transform(data_substrate_sam)
delta_sam,freq = fourie_transform(data_sam)
delta_ref,freq = fourie_transform(data_ref)

delta_ref_substrate_offset_corrected = phase_correction(freq,delta_substrate_ref)
delta_sam_substrate_offset_corrected = phase_correction(freq,delta_substrate_sam)



phase_diff_corrected = delta_ref_substrate_offset_corrected - delta_sam_substrate_offset_corrected
#phase_diff_corrected2 = delta_ref_offset_corrected - delta_sam_offset_corrected
#phase_diff_corrected = (phase_diff_corrected1 + phase_diff_corrected2)/2



c= 3*10**8
d_sub = 700*10**-9
omega=2*np.pi*freq
epsilon_0 = 8.854157 * 10**-12
n=1+c*((phase_diff_corrected))/(omega*d_sub)
print(n)

plt.figure()
plt.plot(freq/10**12,n,color='blue')
plt.plot(freq/10**12,n.imag,color='Brown',label='Imaginary part')
plt.plot(freq/10**12,n.real,color='pink',label='Real part')
plt.xlabel('frequency(THz)')
plt.ylabel('Refractive index')
plt.title('Refractive index of substrate')



sigma = (T_substrate_sam - T_sam)*epsilon_0*c*(1+n)/(T_sam * d_sub)
print(sigma)
plt.figure()
plt.plot(freq/10**12,sigma,color='blue')
plt.plot(freq/10**12,sigma.imag,color='black',label='Imaginary part')
plt.plot(freq/10**12,sigma.real,color='orange',label='Real part')
plt.xlabel('frequency(THz)')
plt.ylabel('Cnductivity')
plt.title('conductivity of film')


#plt.ylim(0.95,3.05)
plt.show()










