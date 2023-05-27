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
from scipy import signal

mpl.rcParams["axes.grid"] = True

def phase_correction(data_fd):
    freq_range = [0.5, 1.0]
    freqs = data_fd[:, 0].real
    freq_slice_idx = (freqs >= freq_range[0]) * (freqs <= freq_range[1])
    phase = np.unwrap(np.angle(data_fd[:, 1]))
    
    z = np.polyfit(freqs[freq_slice_idx], phase[freq_slice_idx], 1)
    
    # phase_corrected = phase - z[1]
    phase_corrected = freqs*z[0]
    
    plt.figure()
    #plt.plot(freq_cp,phase,label="phase b4 correction")
    #plt.plot(freq_cp,phase_corrected,label="phase after correction")
    #plt.plot(freq_cp[freq_slice_idx],phase[freq_slice_idx],label="fit range")
    plt.plot(freqs, freqs*z[0], label="fit")
    plt.ylabel("Unwrapped phase(rad)")
    plt.xlabel("Frequency(Thz")
    plt.legend()

    return phase_corrected

def do_fft(data_td):
    n = data_td[:, 1].size
    timestep = 0.05
    freqs = np.fft.fftfreq(n, d=timestep)
    pos_slice = freqs > 0
    
    Y = np.fft.fft(data_td[:, 1])
    
    data_fd = np.array([freqs[pos_slice], Y[pos_slice]]).T
    
    return data_fd

def time_window(data_td, win_start=None, win_width=15, type_=""):
    t = data_td[:, 0].real - data_td[0, 0].real
    dt = np.mean(np.diff(t))
    
    c = data_td[0,1]
    data_td[:,1] = data_td[:,1] - c  
    
    plt.figure()
    plt.plot(t, data_td[:,1], label=f"Unwindowed {type_}")
    plt.xlabel("Time (ps)")
    plt.ylabel("Amplitude (arb. u.)")
    
    win_width = int(win_width/dt)
    
    if win_start is not None:
        win_start = int(win_start / dt)
    else:
        win_start = np.argmax(np.abs(data_td[:,1])) - win_width // 2
        if win_start < 0:
            win_start = 0
    
    window_axis = signal.windows.tukey(win_width, alpha=0.70)
    zero_pad0 = np.zeros(win_start)
    window_axis = np.concatenate((zero_pad0, window_axis))
    zero_pad1 = np.zeros(len(t) - win_width - win_start)
    window_axis = np.concatenate((window_axis, zero_pad1))
    
    data_td[:,1] = data_td[:,1] * window_axis
    
    plt.plot(t, window_axis*np.max(data_td[:,1]), label="Window")
    plt.plot(t, data_td[:,1], label=f"Windowed {type_}")
    plt.legend()
    plt.show()
    
    return data_td
    
    
    

sub_sam_td = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-24-26.874922-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
sub_ref_td = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-24-07.080916-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt")
doped_sam_td = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-48-05.490852-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt")
doped_ref_td = np.loadtxt(r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-47-45.781852-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt")


sub_sam_td = time_window(sub_sam_td, type_="sub_sam_td")
sub_ref_td = time_window(sub_ref_td, type_="sub_ref_td")

doped_sam_td = time_window(doped_sam_td, type_="doped_sam_td")
doped_ref_td = time_window(doped_ref_td, type_="doped_ref_td")

sub_sam_fd = do_fft(sub_sam_td)
sub_ref_fd = do_fft(sub_ref_td)

doped_sam_fd = do_fft(doped_sam_td)
doped_ref_fd = do_fft(doped_ref_td)

freqs = sub_sam_fd[:, 0].real

phi_sub_sam = phase_correction(sub_sam_fd) 
phi_sub_ref = phase_correction(sub_ref_fd)
phi_doped_sam = phase_correction(doped_sam_fd)
phi_doped_ref = phase_correction(doped_ref_fd)

r_sub_sam = np.abs(sub_sam_fd[:, 1])
r_sub_ref = np.abs(sub_ref_fd[:, 1])
r_doped_sam = np.abs(doped_sam_fd[:, 1])
r_doped_ref = np.abs(doped_ref_fd[:, 1])

sub_sam_fd = np.array([freqs, r_sub_sam * np.exp(-1j*phi_sub_sam)]).T
sub_ref_fd = np.array([freqs, r_sub_ref * np.exp(-1j*phi_sub_ref)]).T
doped_sam_fd = np.array([freqs, r_doped_sam * np.exp(-1j*phi_doped_sam)]).T
doped_ref_fd = np.array([freqs, r_doped_ref * np.exp(-1j*phi_doped_ref)]).T

T_sub = np.array([freqs, sub_sam_fd[:, 1] / sub_ref_fd[:, 1]]).T
T_doped = np.array([freqs, doped_sam_fd[:, 1] / doped_ref_fd[:, 1]]).T

#T_substrate_sam = T_substrate_sam[freq_range]
#T_sam = T_sam[freq_range]

delta_phi_sub = phase_correction(T_sub)

c = 3*10**8
d_film = 700*10**-9
d_sub = 500 * 10 ** -6
omega = 2*np.pi*freqs.real*10**12
epsilon_0 = 8.854157 * 10**-12
n = 1 + c*delta_phi_sub/(omega*d_sub)

print(n)
alpha = -(2/(100*d_sub))*np.log((n+1)**2*np.abs(T_sub[:, 1])/(4*n))
kappa = c*alpha / (2*omega)

plt.figure()
plt.plot(freqs, n, color='Blue', label='Real part')
plt.plot(freqs, kappa, color='Red', label='Imaginary part')
plt.xlabel('Frequency (THz)')
plt.ylabel('Refractive index')
plt.title('Refractive index of substrate')
plt.legend()

plt.figure()
plt.title("Absorption coefficient")
plt.plot(freqs, alpha, color='Red')
plt.xlabel('Frequency (THz)')
plt.ylabel('Absorption coefficient ($cm^{-1}$)')
plt.legend()

sigma = (T_sub[:, 1] - T_doped[:, 1])*epsilon_0*c*(1+n)/(T_doped[:, 1] * d_film)
print(sigma)

plt.figure()
plt.plot(freqs, sigma.real,color='orange', label='Real part')
plt.plot(freqs, sigma.imag,color='black', label='Imaginary part')
plt.xlabel('Frequency (THz)')
plt.ylabel('Conductivity (S/m)')
plt.title('Conductivity of film')
plt.legend()

#plt.ylim(0.95,3.05)
plt.show()










