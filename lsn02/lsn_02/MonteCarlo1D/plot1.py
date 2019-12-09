import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#plot uniform

Nthrows, mean, error = np.loadtxt("data1.dat", usecols=(0,1,2), unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.plot(Nthrows, np.zeros(len(Nthrows)), label='Expected value')
plt.legend()
plt.xlabel('#throws')
plt.ylabel('$I_{MC} - I$')
plt.grid(True)
plt.title('Monte Carlo Integration')
plt.show()
a=np.array([1,2,34,3,25,4,13])

print(error[-1])
#plot importance sampling 1

Nthrows, mean, error = np.loadtxt("data2.dat", usecols=(0,1,2),unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.plot(Nthrows, np.zeros(len(Nthrows)), label='Expected value')
plt.legend()
plt.xlabel('#throws')
plt.ylabel('$I_{MC} - I$')
plt.grid(True)
plt.title('Monte Carlo Integration')
plt.show()

print(error[-1])

#plot importance sampling 2

Nthrows, mean, error = np.loadtxt("data3.dat", usecols=(0,1,2),unpack=True)

plt.errorbar(Nthrows ,mean-1, yerr=error)
plt.plot(Nthrows, np.zeros(len(Nthrows)), label='Expected value')
plt.legend()
plt.xlabel('#throws')
plt.ylabel('$I_{MC} - I$')
plt.grid(True)
plt.title('Monte Carlo Integration')
plt.show()

print(error[-1])
