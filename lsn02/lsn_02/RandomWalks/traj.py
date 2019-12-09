import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

x,y,z = np.loadtxt("lattice1000.dat", usecols=(0,1,2), unpack=True)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(x,y,z)
ax.scatter(x[0],y[0],z[0], label='starting point',color='r')
plt.legend()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Random Walk on a 3D cubic lattice')
plt.show()


x,y,z = np.loadtxt("continuum1000.dat", usecols=(0,1,2), unpack=True)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(x,y,z)
ax.scatter(x[0],y[0],z[0],color='r', label='starting point')
plt.legend()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.title('Random Walk in the continuum')
plt.show()




def f(x,k):
   return k*np.sqrt(x)



t, r2, sigma = np.loadtxt("diffusion.dat", usecols=(0,1,2), unpack=True)

r2 = np.sqrt(r2)
sigma = np.sqrt(sigma)

plt.errorbar(t,r2,yerr=sigma, label='simulation data')
p_opt, p_cov = curve_fit(f, t, r2)
y_fit = f(t,p_opt[0])
plt.plot(t,y_fit, label = 'fitting function', linewidth=3)
plt.title('Random Walk on a 3D cubic lattice')
plt.xlabel('t')
plt.ylabel('$<|r_t|^2>$')
plt.legend()
plt.grid(True)
plt.show()

k=round(p_opt[0],6)
sigmak=round(np.sqrt(p_cov[0][0]),6)
print("Lattice: "+str(k)+" +/- "+str(sigmak))




t, r2, sigma = np.loadtxt("diffusion.dat", usecols=(0,3,4), unpack=True)

r2 = np.sqrt(r2)
sigma = np.sqrt(sigma)

plt.errorbar(t,r2,yerr=sigma, label='simulation data')
p_opt, p_cov = curve_fit(f, t, r2)
y_fit = f(t,p_opt[0])
plt.plot(t,y_fit, label = 'fitting function', linewidth=3)
plt.title('Random Walk in the continuum')
plt.xlabel('t')
plt.ylabel('$<|r_t|^2>$')
plt.grid(True)
plt.legend()
plt.show()

k=round(p_opt[0],6)
sigmak=round(np.sqrt(p_cov[0][0]),6)
print("Continuum: "+str(k)+" +/- "+str(sigmak))








