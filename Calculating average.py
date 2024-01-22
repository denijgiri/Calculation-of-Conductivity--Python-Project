# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:22:27 2024

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt


def calculate_average(data):
   averages = []

   for i in range(10):
         current_data = np.loadtxt(data[i])
         current_average = np.mean(current_data[:,1])
         averages.append(current_average)
   total_average = np.mean(averages)
   return total_average
        
         
data_sam = [r"D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-40-07.125625-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-40-24.115630-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-40-40.875626-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-40-57.645625-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-41-14.425624-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-41-31.215624-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-41-48.025631-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-42-04.865627-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-42-21.655631-REF--X_-5.000 mm-Y_15.000 mm.txt","D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-42-38.375626-REF--X_-5.000 mm-Y_15.000 mm.txt"]
results = calculate_average(data_sam)
print(results)

