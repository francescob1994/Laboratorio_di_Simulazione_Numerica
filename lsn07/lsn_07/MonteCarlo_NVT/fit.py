import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(t,a,b,c):  # defining the fitting function
    return a * np.exp(-t/b) + c


mu, mp = np.loadtxt("up.mean", usecols=(0,1), unpack=True)
u0ui = np.zeros(len(mu)-1)
p0pi = np.zeros(len(mp)-1)
for i in range(len(mu)-1):  #len(mu)-1
   u0ui[i] = mu[0]*mu[i+1]
   p0pi[i] = mp[0]*mp[i+1]

covu, covp = np.loadtxt("cov.mean", usecols=(0,1), unpack=True)
coru = np.zeros(len(covu))
corp = np.zeros(len(covp))
for i in range(len(covu)):  #len(covu)
   coru[i] = covu[i]-u0ui[i]
   corp[i] = covp[i]-p0pi[i]

t = np.arange(len(covu))#len(covu)


p_opt, p_cov = curve_fit(f, t, coru)  
coru_fit = f(t,p_opt[0],p_opt[1],p_opt[2])
print("optimized parameters [a,b,c] =")
print(p_opt)
print("parameters uncertainty =")
print(np.sqrt(np.diagonal(p_cov)))



plt.subplot(121)
plt.plot(t,coru_fit,label='fit') # plotting fitted function
plt.plot(t,coru)
plt.xlabel('t')
plt.ylabel('$<u_0\cdot u_t> - <u_0><u_t>$')
plt.legend()
plt.title('Energy correlations')
plt.grid(True)


p_opt, p_cov = curve_fit(f, t, corp)  
corp_fit = f(t,p_opt[0],p_opt[1],p_opt[2])
print("optimized parameters [a,b,c] =")
print(p_opt)
print("parameters uncertainty =")
print(np.sqrt(np.diagonal(p_cov)))

plt.subplot(122)
plt.plot(t,corp_fit,label='fit') # plotting fitted function
plt.plot(t,corp)
plt.xlabel('t')
plt.ylabel('$<p_0\cdot p_t> - <p_0><p_t>$')
plt.title('Pressure correlations')
plt.legend()
plt.grid(True)


plt.show()





