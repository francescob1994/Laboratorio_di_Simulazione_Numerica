#leggo il file finalH.dat e trovo il minimo di H e i corrispondenti parametri mu e sigma.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mu, sigma, H ,errH= np.loadtxt("finalH.dat", usecols=(0,1,2,3), delimiter='   ', unpack='true')

minH = H[0]
ind = 0
pos = 1
while pos < len(H):
   if H[pos] < minH:
      minH = H[pos]
      ind = pos
   pos = pos + 1

MU = mu[ind]
SIGMA = sigma[ind]
err = errH[ind]

print("Minimum Energy = "+str(minH)+"+/-"+str(err))
print("Mu = "+str(MU)+"   Sigma = "+str(SIGMA))


file = open("opt_parameters.txt", "w")
file.write(str(MU)+"   "+str(SIGMA))
file.close()

norm = 2*SIGMA*np.sqrt(np.pi)*( 1 + np.exp(-(MU**2)/(SIGMA**2)) );

x = np.arange(-2,2, 0.02)
psi = ( np.exp(-((x+MU)**2)/(2*SIGMA**2)) + np.exp(-((x-MU)**2)/(2*SIGMA**2)) )
psi2 = (psi**2)/norm
V = x**4 - (5./2)*x**2

plt.plot(x,psi,label = "psi")
plt.plot(x,psi2,label = "psi2")
plt.plot(x,V, label = "V" )
plt.xlabel('x')
plt.grid(True)
plt.legend()
plt.show()



