import matplotlib.pyplot as plt
import numpy as np

u, p = np.loadtxt("up.dat", usecols=(0,1), unpack=True)
L = np.array([10,20,25,50,100,125,200,250,500,625,1000,1250, 2000,2500,5000])
Ntot = 10000 #per sim

nsim = int(len(u)/Ntot)
u = np.reshape(u, (nsim,Ntot)) #reshape ditribuisce i dati sulle righe
p = np.reshape(p, (nsim,Ntot))
print(u.shape)

def Err(x,l,Ntot):
   Nblk = Ntot/l
   glob_ave = 0
   glob_ave2 = 0
   for i in range(int(Nblk)):
      sum = 0
      for j in range(l):
         k=j+i*l
         sum += x[k]

      glob_ave += sum/l
      glob_ave2 += (sum/l)**2
   
   glob_ave /= Nblk 
   glob_ave2 /= Nblk
   return np.sqrt( (glob_ave2 - glob_ave**2)/(Nblk-1) )


erru = np.zeros(len(L))
errp = np.zeros(len(L))

for k in range(nsim):
   for l in range(len(L)):
      erru[l] += Err(u[k,:],L[l],Ntot)
      errp[l] += Err(p[k,:],L[l],Ntot)

erru = erru/nsim
errp = errp/nsim




plt.plot(L,erru, label='u')
plt.plot(L,errp, label='p')
plt.legend()
plt.grid(True)
plt.title("Mean error vs Block dimension")
plt.xlabel('L')
plt.ylabel('$\sigma$')
plt.show()
