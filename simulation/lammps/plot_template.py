import numpy as np
import matplotlib.pyplot as plt

dataY = np.loadtxt("output_cnt_length.dat", skiprows=2 , usecols=1, dtype='double')
dataX = np.loadtxt("output_cnt_length.dat", skiprows=2 , usecols=0, dtype='int')
dataX = dataX / 1000

plt.plot(dataX, dataY)
plt.xlabel('t (ps)')
plt.ylabel('L (Angstroms)')
plt.show()
