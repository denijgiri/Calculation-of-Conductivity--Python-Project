# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 12:54:24 2023

@author: Denij
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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
    
    #plt.figure()
    #plt.plot(freq_cp,phase,label="phase b4 correction")
    #plt.plot(freq_cp,phase_corrected,label="phase after correction")
    #plt.plot(freq_cp[freq_slice_idx],phase[freq_slice_idx],label="fit range")
    #plt.plot(freqs, freqs*z[0], label="fit")
    #plt.ylabel("Unwrapped phase(rad)")
    #plt.xlabel("Frequency(Thz)")
    #plt.legend()

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
    
    #plt.figure()
    #plt.plot(t, data_td[:,1], label=f"GaAsTe_wafer19075")
    #plt.xlabel("Time (ps)")
    #plt.ylabel("Amplitude (db)")
    
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
    
    #plt.plot(t, window_axis*np.max(data_td[:,1]), label="Window")
    #plt.plot(t, data_td[:,1], label=f"GaAsTe_wafer19075")
    #plt.legend()
    #plt.show()
    
    return data_td
    
    
def analysis(ref_sub_path, sub_path, ref_sam_path, sam_path):
    sub_sam_td = np.loadtxt(sub_path)
    sub_ref_td = np.loadtxt(ref_sub_path)
    doped_sam_td = np.loadtxt(sam_path)
    doped_ref_td = np.loadtxt(ref_sam_path)
    
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
    
    delta_phi_sub = -phi_sub_sam + phi_sub_ref
    
    c = 3*10**8
    d_film = 700*10**-9
    d_sub = 500 * 10 ** -6
    omega = 2*np.pi*freqs.real*10**12
    epsilon_0 = 8.854157 * 10**-12
    n = 1 + c*delta_phi_sub/(omega*d_sub)
    
    print(n)
    alpha = -(2/(100*d_sub))*np.log((n+1)**2*np.abs(T_sub[:, 1])/(4*n))

    sigma = (T_sub[:, 1] - T_doped[:, 1])*epsilon_0*c*(1+n)/(T_doped[:, 1] * d_film)
        
    
    return freqs, sigma
    

sub_paths =  [r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-18-31.067459-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-19-10.484467-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-19-50.143466-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-20-30.305680-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-21-09.316397-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-21-48.975981-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-22-28.871915-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-23-07.948921-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-23-47.450916-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-24-26.874922-GaAs_sub_wafer25_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt"]

sub_ref_paths =  [r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAs_wafer25_sub/2023-01-19T15-18-11.052862-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-18-50.683459-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-19-30.325466-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-20-10.178856-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-20-49.811071-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-21-28.591924-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-22-08.573348-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-22-48.293922-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-23-27.558922-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAs_wafer25_sub\2023-01-19T15-24-07.080916-GaAs_sub_wafer25_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt"]

doped_paths = [r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-42-06.856860-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-42-46.419859-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-43-25.836851-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-44-05.321852-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-44-44.927851-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-45-25.068852-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-46-06.032854-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-46-46.674860-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-47-26.305859-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-48-05.490852-GaAsTe_wafer19075_2RH_100_avg-sample-X_24.000 mm-Y_14.000 mm.txt"]
              
doped_ref_paths = [r"C:/Users/Denij/OneDrive/Desktop/Semiconductors/GaAsTe_wafer19075/2023-01-19T15-41-46.977074-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-42-26.551860-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-43-06.189851-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-43-45.550850-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-44-24.990852-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-45-04.965860-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-45-45.107853-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-46-26.197863-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-47-06.391853-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt",r"C:\Users\Denij\OneDrive\Desktop\Semiconductors\GaAsTe_wafer19075\2023-01-19T15-47-45.781852-GaAsTe_wafer19075_2RH_100_avg-reference-X_-8.000 mm-Y_14.000 mm.txt"]

results = []
for i in range(10):
     sub_path = sub_paths[i]
     sub_ref_path = sub_ref_paths[i]
     doped_path = doped_paths[i]
     doped_ref_path = doped_ref_paths[i]
     freqs, sigma = analysis(sub_path,sub_ref_path,doped_path, doped_ref_path)

     results.append(sigma)
    
    
#print(sigma)
average = np.mean(results,axis=0)
average = average / 100
std = np.std(results,axis=0)
print(results)

plt.figure()
plt.plot(freqs, average.real, color = 'red')
plt.fill_between(freqs, average.real + std.real, average.real - std.real, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.figure()
plt.plot(freqs, average.imag, color = 'blue')
plt.fill_between(freqs,average.real + std.real, average.imag - std.real, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.show()
plt.xlabel('frequency')
plt.ylabel('sigma')
plt.legend()



#plt.ylim(0.95,3.05)
#plt.show()
