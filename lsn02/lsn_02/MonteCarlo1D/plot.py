import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#plot uniform

Nthrows, mean, error = np.loadtxt("data1.dat", usecols=(0,1,2), unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.xlabel('#throws')
plt.ylabel('Integral')
plt.grid(True)
plt.title('Montecalro Integration')
plt.show()

#plot importance sampling 1

Nthrows, mean, error = np.loadtxt("data2.dat", usecols=(0,1,2),unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.xlabel('#throws')
plt.ylabel('Integral')
plt.grid(True)
plt.title('Montecalro Integration')
plt.show()

#plot importance sampling 2

Nthrows, mean, error = np.loadtxt("data3.dat", usecols=(0,1,2),unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.xlabel('#throws')
plt.ylabel('Integral')
plt.grid(True)
plt.title('Montecalro Integration')
plt.show()
