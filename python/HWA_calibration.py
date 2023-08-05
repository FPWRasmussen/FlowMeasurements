import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

def HWA_calibration():
	velocity = np.arange(0, 20+1, 2) 
	
	mean_voltage = np.zeros(len(velocity))
	
	for i, u in enumerate(velocity):
		filename = "HWA/Calibration_"+f"0{u}".zfill(3)
		t, V = np.loadtxt(filename, skiprows=23).T
		mean_voltage[i] = V.mean()
	
	# fitting the line
	c = Polynomial.fit(mean_voltage, velocity, 4)
	V_fit = np.linspace(min(mean_voltage), max(mean_voltage), 100)
	
	if __name__ == "__main__":
		plt.figure(1)
		plt.plot(c(V_fit), V_fit, '-k')
		plt.plot(velocity,mean_voltage, 'mo', markersize=6)
		plt.grid(True)
		plt.xlabel('u [m/s]')
		plt.ylabel('E [V]')
		plt.xlim([-2.5, 22.5])
		plt.show()
	
	
	return c


if __name__ == "__main__":
	c =HWA_calibration()