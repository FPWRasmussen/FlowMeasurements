import numpy as np
import matplotlib.pyplot as plt
from HWA_calibration import HWA_calibration

def HWA(angle):
	angle = str(angle).zfill(2)
	pos = np.arange(-36, 44, 4)
	voltage = np.zeros(len(pos))
	velocity_rms = np.zeros(len(pos))
	
	c = HWA_calibration()

	for i, p in enumerate(pos):
		if p >= 0:
			pos_str = "Meausurement_"+"p"+f"{p}".zfill(2)+f"_{angle}".zfill(2)
		else:
			pos_str = "Meausurement_"+"m"+f"{-1*p}".zfill(2)+f"_{angle}".zfill(2)
		dt = np.loadtxt("HWA/"+pos_str, skiprows=23)
	
		voltage[i] = np.mean(dt[:,1])
		velocity_rms[i] = np.sqrt(np.mean((c(dt[:,1])-c(voltage[i]))**2))

	
	
	return c, pos, voltage, velocity_rms
	
	
if __name__ == "__main__":
	c, pos, voltage, voltage_rms = HWA(15)

	plt.figure(figsize=(7,5))
	plt.plot(c(voltage),pos, color = "black", linestyle = "dotted", linewidth=2)
	plt.scatter(c(voltage), pos, color = "blue", marker = "^")
	plt.ylabel("Position [mm]")
	plt.xlabel("Velocity [m/s]")
	plt.grid()
	plt.tight_layout()
	
	plt.figure(figsize=(7,5))
	plt.plot( (voltage_rms), pos,color = "black", linestyle = "dotted", linewidth=2)
	plt.scatter((voltage_rms), pos, color = "blue", marker = "^")
	plt.ylabel("Position [mm]")
	plt.xlabel("RMS [m/s]")
	plt.grid()
	plt.tight_layout()