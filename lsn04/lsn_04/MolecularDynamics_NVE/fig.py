import numpy as np
import matplotlib.pyplot as plt

kb = 1.380649*(10**(-23))
amu = 1.66054*(10**(-27))

sigAr = 0.34*(10**(-9))
epsAr = 120*kb
mAr = 39.948*amu

sigKr = 0.364*(10**(-9))
epsKr = 164*kb
mKr = 83.798*amu

#Kinetic energy

#Nblocks, K, errK = np.loadtxt("gasK.out", usecols=(0,2,3), unpack=True)
Nblocks, K, errK = np.loadtxt("ave_ekin.out", usecols=(0,2,3), unpack=True)

K_Ar = K*epsAr
errK_Ar = errK*epsAr

K_Kr = K*epsKr
errK_Kr = errK*epsKr

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.errorbar(Nblocks ,K_Ar, yerr=errK_Ar)
plt.xlabel('#blocks')
plt.ylabel('$K/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,K_Kr, yerr=errK_Kr)
plt.xlabel('#blocks')
plt.ylabel('$K/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Kinetic energy per particle')
plt.show()

#Potential energy

#Nblocks, U, errU = np.loadtxt("gasU.out", usecols=(0,2,3), unpack=True)
Nblocks, U, errU = np.loadtxt("ave_epot.out", usecols=(0,2,3), unpack=True)

U_Ar = U*epsAr
errU_Ar = errU*epsAr

U_Kr = U*epsKr
errU_Kr = errU*epsKr

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.errorbar(Nblocks ,U_Ar, yerr=errU_Ar)
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,U_Kr, yerr=errU_Kr)
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Potential energy per particle')
plt.show()

#Total energy

#Nblocks, E, errE = np.loadtxt("gasE.out", usecols=(0,2,3), unpack=True)
Nblocks, E, errE = np.loadtxt("ave_etot.out", usecols=(0,2,3), unpack=True)

E_Ar = E*epsAr
errE_Ar = errE*epsAr

E_Kr = E*epsKr
errE_Kr = errE*epsKr

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.errorbar(Nblocks ,E_Ar, yerr=errE_Ar)
plt.xlabel('#blocks')
plt.ylabel('$E/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,E_Kr, yerr=errE_Kr)
plt.xlabel('#blocks')
plt.ylabel('$E/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Total energy per particle')
plt.show()

#Temperature

#Nblocks, T, errT = np.loadtxt("gasT.out", usecols=(0,2,3), unpack=True)
Nblocks, T, errT = np.loadtxt("ave_temp.out", usecols=(0,2,3), unpack=True)

T_Ar = T*epsAr/kb
errT_Ar = errT*epsAr/kb

T_Kr = T*epsKr/kb
errT_Kr = errT*epsKr/kb

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.errorbar(Nblocks ,T_Ar, yerr=errT_Ar)
plt.plot(Nblocks, np.array([1.2*epsAr/kb for i in range(len(Nblocks))]), label='fixed temperature')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$T$ (K)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,T_Kr, yerr=errT_Kr)
plt.plot(Nblocks, np.array([1.2*epsKr/kb for i in range(len(Nblocks))]), label='fixed temperature')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$T$ (K)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Temperature')
plt.show()


#Pressure

#Nblocks, P, errP = np.loadtxt("gasP.out", usecols=(0,2,3), unpack=True)
Nblocks, P, errP = np.loadtxt("ave_pres.out", usecols=(0,2,3), unpack=True)

P_Ar = P*epsAr/sigAr**3
errP_Ar = errP*epsAr/sigAr**3

P_Kr = P*epsKr/sigKr**3
errP_Kr = errP*epsKr/sigKr**3

plt.figure(figsize=(16,8))

plt.subplot(121)
plt.errorbar(Nblocks ,P_Ar, yerr=errP_Ar)
plt.xlabel('#blocks')
plt.ylabel('$P$ (Pa)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,P_Kr, yerr=errP_Kr)
plt.xlabel('#blocks')
plt.ylabel('$P$ (Pa)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Pressure')
plt.show()




