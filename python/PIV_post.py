import numpy as np
import matplotlib.pyplot as plt
import os
from HWA import HWA


def PIV_post(deg):
	deg = deg
	dat = np.loadtxt(f"PIV/{deg}deg.dat", skiprows=3)
	dat_rms = np.loadtxt(f"PIV/{deg}deg2.dat", skiprows=3)
	
	x, y, Vx_rms, Vy_rms, isValid = dat_rms.T
	
	x, y, Vx, Vy, isValid = dat.T
	x, y = x[0:len(np.unique(x))], y[0::len(np.unique(x))]
	#%%
	X, Y = np.meshgrid(x,y)
	VX = Vx.reshape(X.shape)
	VY = Vy.reshape(X.shape)
	isValid = isValid.reshape(X.shape)
	isInvalid = np.logical_not(isValid)
	masked_VX = np.ma.masked_array(VX, mask=isInvalid)
	masked_VY = np.ma.masked_array(VY, mask=isInvalid)
	magnitude = np.sqrt(VX ** 2 + VY ** 2)
	masked_magnitude = np.ma.masked_array(magnitude, mask=isInvalid)
	
	VX_rms = Vx_rms.reshape(X.shape)
	VY_rms = Vy_rms.reshape(X.shape)
	masked_VX_rms = np.ma.masked_array(VX_rms, mask=isInvalid)
	masked_VY_rms = np.ma.masked_array(VY_rms, mask=isInvalid)
	magnitude_rms = np.sqrt(VX_rms ** 2 + VY_rms ** 2)
	masked_magnitude_rms = np.ma.masked_array(magnitude_rms, mask=isInvalid)
	
	
	
	# Subsample the data for cuiver
	subsample = 2
	X_sub = X[::subsample, ::subsample]
	Y_sub = Y[::subsample, ::subsample]
	VX_sub = VX[::subsample, ::subsample]
	VY_sub = VY[::subsample, ::subsample]
	magnitude_sub = magnitude[::subsample, ::subsample]
	
	# levels = np.linspace(0, np.max(magnitude), 10)  # Generate 10 levels from -1 to 1
	
	c, posy, voltage, velocity_rms = HWA(deg)
	
	HWA_posx = 150
	posx = np.ones(len(posy))*HWA_posx
	
	PIV_xindex = min(range(len(x)), key=lambda i: abs(x[i]-HWA_posx))
	
	# levels = np.linspace(0,15,1000)
	fig, ax = plt.subplots(1,3, figsize = [16,6],gridspec_kw={'width_ratios': [3, 1, 1]})
	
	shift = 0
	
	contour = ax[0].contourf(x, y, masked_magnitude, cmap = "jet", levels=1000)
	ax[0].quiver(X_sub, Y_sub, VX_sub, VY_sub, magnitude_sub, cmap = "gray", angles="xy", scale_units="xy", scale=1, width=0.001)
	ax[0].scatter(posx, posy+shift, c=c(voltage), cmap = "jet", edgecolor ="black", label ="HWA", vmin=masked_magnitude.min(), vmax=masked_magnitude.max())
	ax[0].set_xlabel("x [mm]")
	ax[0].set_ylabel("y [mm]")
	ax[0].legend(loc="upper left")
	fig.colorbar(mappable = contour, cmap = "jet", label = "Velocity [m/s]")
	
	ax[1].plot(c(voltage), posy+shift, label = "HWA", color='black', ls='dashed', marker='o')
	ax[1].plot(masked_magnitude[:, PIV_xindex], y, label = "PIV", color='blue', ls='dotted')
	ax[1].set_xlabel("Velocity [m/s]")
	ax[1].set_ylabel("y [mm]")
	ax[1].legend()
	ax[1].grid()
	
	ax[2].plot((velocity_rms), posy+shift, label = "HWA", color='black', ls='dashed', marker='o')
	ax[2].plot(masked_magnitude_rms[:, PIV_xindex], y, label = "PIV", color='blue', ls='dotted')
	ax[2].set_xlabel("RMS [m/s]")
	ax[2].set_ylabel("y [mm]")
	ax[2].legend()
	ax[2].grid()
	
	plt.tight_layout()
	
	
	
	
	# curl = np.gradient(v, x, axis=1) - np.gradient(u, y, axis=0)
	
	# # Plot the curl
	# plt.figure(figsize=(8, 6))
	# plt.contourf(X, Y, curl, cmap='viridis')
	# plt.colorbar()
	# plt.xlabel('X')
	# plt.ylabel('Y')
	# plt.title('Curl of the Flow Field')
	# plt.show()
	
	return x, y, masked_magnitude

if __name__ == "__main__":
	x, y, masked_magnitude = PIV_post(15)