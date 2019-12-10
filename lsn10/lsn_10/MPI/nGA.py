import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# show best path

data = np.loadtxt("FinalSeq.dat", usecols=(0))

b = data[0]
ll = data[1]
l= data[2]
seq = np.zeros(len(data)-3)
for i in range(len(data)-3):
   seq[i] = data[i+3]


if ll==1:
   L="L2"
else: 
   L="L1"


if b==1:
   file = "square.dat"
   cit="square"
else:
   file = "circumference.dat"
   cit="circumference"
C = np.loadtxt(file, usecols=(0,1))
C1 = np.zeros((C.shape[0]+1, C.shape[1]))

for i in range(len(seq)):
   ind = int(seq[i])
   C1[i,:] = C[ind,:]

C1[len(seq),:]=C1[0,:]

plt.scatter(C1[:,0],C1[:,1], label='cities')
plt.plot(C1[:,0],C1[:,1])
plt.legend()
plt.grid(True)
plt.title('Best path: '+L+' = '+str(l))
plt.show()
   

# find timescale optimization

t, lmin, lave = np.loadtxt("length.dat", usecols=(0,1), unpack = True )
def f(x,a,l,b):
   return a*np.exp(-x/l)+b

p_opt, p_cov = curve_fit(f, t, lmin)
y_fit = f(t,p_opt[0],p_opt[1],p_opt[2])
plt.subplot(121)
plt.plot(t,y_fit,label='fitting function')
plt.plot(t,lmin, label='l_min')
plt.legend()
plt.title("tau = "+str(int(p_opt[1])))
plt.xlabel('generation')
plt.ylabel('l_min')
plt.grid(True)
plt.show()

