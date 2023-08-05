import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate

data = np.loadtxt("HWA/CorrelationTest", skiprows=23)

# Split the data into t and u
t, u = data.T

meanvoltage = u.mean()
dev = u - meanvoltage

corr = correlate(dev, dev, mode='full')
corr2 = corr[corr.size//2:]
corr3 = corr2 / corr2[0]
plt.plot(corr3[:1000])
plt.grid()
plt.xlabel("Time [ms]")
plt.ylabel("Autocorrelation coefficient")
plt.tight_layout()
plt.show()
