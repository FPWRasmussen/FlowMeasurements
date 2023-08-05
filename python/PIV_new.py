import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imread
from scipy.signal import correlate
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import zoom
from PIV_post import PIV_post

# load the images
im = imread("PIV/0deg/B00001.tif")

x,y,vx,vy, mask = np.loadtxt("PIV/0deg/0deg.dat",skiprows=3).T
x, y = x[0:len(np.unique(x))], y[0::len(np.unique(x))]

mid_point = len(im) // 2 

a = im[:mid_point,:]
b = im[mid_point:,:]


#%%
if False:
	fig, ax = plt.subplots()
	im = ax.imshow(a, cmap=plt.cm.gray)
	current_image = a
	break_loop = [False]
	
	def on_close(event):
		break_loop[0] = True
		
	fig.canvas.mpl_connect('close_event', on_close)
	
	while not break_loop[0]:
		# Display the current image
		im.set_data(current_image)
		plt.draw()
		plt.pause(1)
		if current_image is a:
			current_image = b
		else:
			current_image = a
		

# Coordinate 
cord = {"A": [10,5]}

def vel_field(curr_frame, next_frame, win_size):
	x,y,vx,vy, mask = np.loadtxt("PIV/0deg/0deg.dat",skiprows=3).T
	ys = np.arange(0, curr_frame.shape[0], win_size)
	xs = np.arange(0, curr_frame.shape[1], win_size)
	dys = np.zeros((len(ys), len(xs)))
	dxs = np.zeros((len(ys), len(xs)))
	SNR = np.zeros((len(ys), len(xs)))
	T = 73.7
	
	x_scaling = np.abs((x[-1]-x[0])/(xs[-1]-xs[0]))/(T)*1000
	y_scaling = np.abs((y[-1]-y[0])/(ys[-1]-ys[0]))/(T)*1000
	
	for iy, y in enumerate(ys):
		for ix, x in enumerate(xs):
			int_win = curr_frame[y : y + win_size, x : x + win_size]
			search_win = next_frame[y : y + win_size, x : x + win_size]
			cross_corr = correlate(search_win - search_win.mean(), int_win - int_win.mean(), method="fft")
			dys[iy, ix], dxs[iy, ix] = (np.unravel_index(np.argmax(cross_corr), cross_corr.shape)- np.array([win_size, win_size]) + 1)
			
			dys[iy, ix] *= y_scaling
			dxs[iy, ix] *= x_scaling
			
# 			largest_value = np.partition(cross_corr.flatten(), -1)[-1]
# 			second_value = np.partition(cross_corr.flatten(), -2)[-2]   

			# Find the maximum value
			max_value = np.amax(cross_corr)
			
			# Create a mask to exclude values adjacent to the maximum
			mask = np.ones_like(cross_corr, dtype=bool)
			max_indices = np.argwhere(cross_corr == max_value)
			for i, j in max_indices:
			    mask[max(i - 1, 0):min(i + 2, cross_corr.shape[0]), max(j - 1, 0):min(j + 2, cross_corr.shape[1])] = False
			
			mask = np.logical_not(mask)
			masked_array = np.ma.array(cross_corr, mask=mask)
			second_max_value = np.ma.max(masked_array)

			SNR[iy, ix] = max_value/second_max_value				  
			
			for point in cord:
				if False:
					if cord[point] == [iy,ix]:
						fig = plt.figure()
# 						ax = fig.add_subplot(projection="3d")
						Y, X = np.meshgrid(np.arange(cross_corr.shape[0]), np.arange(cross_corr.shape[1]))
# 						ax.plot_surface(Y, X, cross_corr, cmap='rainbow', linewidth=0.2)  # type: ignore
						plt.contourf(Y,X, cross_corr, cmap="gray")
						plt.colorbar(label="Correlation")
						max_index = np.unravel_index(np.argmax(cross_corr), cross_corr.shape)
						second_max_index = np.unravel_index(np.argmax(masked_array), masked_array.shape)
						plt.scatter(max_index[1], max_index[0],marker="x", label = "Largest correlation")
						plt.scatter(second_max_index[1], second_max_index[0],marker="x", label = "Second largest correlation")
						plt.title(f"Correlation map example")
						plt.legend()
						plt.show()
				
				
		  
	# draw velocity vectors from the center of each window
	ys = ys + win_size / 2
	xs = xs + win_size / 2
	
	return xs, ys, dxs, dys, SNR

win_size = 32

ys = np.arange(0, a.shape[0], win_size)
xs = np.arange(0, a.shape[1], win_size)
Vx = np.zeros((len(ys), len(xs)))
Vy = np.zeros((len(ys), len(xs)))
SNR_avg = np.zeros((len(ys), len(xs)))

n = np.arange(1,21,1)

for i in n:
	im = imread(f"PIV/0deg/B00"+str(i).zfill(3)+".tif")
	mid_point = len(im) // 2 
	a = np.flip(im[:mid_point,:],axis=0)
	b = np.flip(im[mid_point:,:],axis=0)
# 	a = np.ma.masked_array(a, mask=scaled_mask).filled()
# 	b = np.ma.masked_array(b, mask=scaled_mask).filled()
	xs, ys, dxs, dys, SNR = vel_field(a, b, win_size)
	Vx += dxs
	Vy += dys
	SNR_avg += SNR

Vx /= len(n)
Vy /= len(n)
SNR_avg /= len(n)

mask = np.flip(np.logical_not(np.reshape(mask, [len(y),len(x)])),axis=0)
scaled_mask = zoom(mask, (Vx.shape[0] / mask.shape[0], Vx.shape[1] / mask.shape[1]),order=0)
SNR_mask = SNR < 2
total_mask = np.logical_or(scaled_mask, SNR_mask)

norm_drs = np.sqrt(Vx ** 2 + Vy ** 2)
norm_drs = np.ma.masked_array(norm_drs, mask=total_mask)

#%%
xs = np.linspace(x[0],x[-1],len(xs))
ys = np.linspace(y[-1],y[0],len(ys))

# fig, ax = plt.subplots(1,2, figsize = [16,6])
# contour1 = ax[0].contourf(xs,ys, norm_drs, alpha = 1, cmap = "jet", levels = 1000)
# # ax[0].set_title(f"Flow Field")
# ax[0].quiver(xs, ys, Vx,Vy, norm_drs, cmap="gray", angles="xy",  scale_units="xy",scale=2)
# ax[0].set_aspect("equal")
# ax[0].set_xlabel("x position [mm]")
# ax[0].set_ylabel("y position [mm]")
# fig.colorbar(mappable = contour1, cmap = "jet", label = "Velocity [m/s]")
# contour2 = ax[1].contourf(xs,ys, SNR, alpha = 1, cmap = "jet", levels = 1000)
# ax[1].set_aspect("equal")
# ax[1].set_xlabel("x position [mm]")
# ax[1].set_ylabel("y position [mm]")
# fig.colorbar(mappable = contour2, cmap = "jet", label = "Signal-to-noise ratio")
# plt.tight_layout()


x, y, masked_magnitude = PIV_post(0)
masked_magnitude = zoom(masked_magnitude, (norm_drs.shape[0] / masked_magnitude.shape[0], norm_drs.shape[1] / masked_magnitude.shape[1]),order=0)
mask = np.logical_not(masked_magnitude > 0)
masked_magnitude = np.ma.masked_array(masked_magnitude, mask=mask)
masked_magnitude = np.flip(masked_magnitude,axis=0)

diff = masked_magnitude-norm_drs 

plt.figure()
plt.contourf(xs, ys,diff , cmap = "jet", levels = 10)
plt.colorbar(cmap = "jet", label = "Difference [m/s]")
plt.xlabel("x position [mm]")
plt.ylabel("y position [mm]")
plt.tight_layout()

