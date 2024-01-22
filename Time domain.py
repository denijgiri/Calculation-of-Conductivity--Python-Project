# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 07:41:28 2024

@author: Denij
"""

import numpy as np
import matplotlib.pyplot as plt

data_sam = np.loadtxt(r"D:/Docs and Py code/MarielenaData/2021_08_24/GaAs_Te 19075/Reference/2021-08-24T10-40-07.125625-REF--X_-5.000 mm-Y_15.000 mm.txt")
c = data_sam[0,1]
data_sam = data_sam[:,1] - data_sam[0,1]
plt.plot(data_sam)      
plt.show()


