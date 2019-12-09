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




#solid state
r, g, errg = np.loadtxt("gr_solid.dat", usecols=(0,1,2), unpack = True)
r_md, g_md, errg_md = np.loadtxt("gdir_solid.dat", usecols=(0,1,2), unpack = True)

rAr = r*sigAr
rAr_md = r_md*sigAr 
rKr = r_md*sigKr
rKr_md = r_md*sigKr

plt.subplot(121)
plt.errorbar(rAr,g,yerr=errg,label="MC")
plt.errorbar(rAr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Argon")

plt.subplot(122)
plt.errorbar(rKr,g,yerr=errg,label="MC")
plt.errorbar(rKr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Kripton")

plt.suptitle("Solid State")
plt.show()




#liquid state
r, g, errg = np.loadtxt("gr_liquid.dat", usecols=(0,1,2), unpack = True)
r_md, g_md, errg_md = np.loadtxt("gdir_liquid.dat", usecols=(0,1,2), unpack = True)

rAr = r*sigAr
rAr_md = r_md*sigAr 
rKr = r_md*sigKr
rKr_md = r_md*sigKr

plt.subplot(121)
plt.errorbar(rAr,g,yerr=errg,label="MC")
plt.errorbar(rAr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Argon")

plt.subplot(122)
plt.errorbar(rKr,g,yerr=errg,label="MC")
plt.errorbar(rKr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Kripton")

plt.suptitle("Liquid State")
plt.show()



#gaseous state
r, g, errg = np.loadtxt("gr_gas.dat", usecols=(0,1,2), unpack = True)
r_md, g_md, errg_md = np.loadtxt("gdir_gas.dat", usecols=(0,1,2), unpack = True)

rAr = r*sigAr
rAr_md = r_md*sigAr 
rKr = r_md*sigKr
rKr_md = r_md*sigKr

plt.subplot(121)
plt.errorbar(rAr,g,yerr=errg,label="MC")
plt.errorbar(rAr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Argon")

plt.subplot(122)
plt.errorbar(rKr,g,yerr=errg,label="MC")
plt.errorbar(rKr_md,g_md,yerr=errg_md,label="MD")
plt.legend()
plt.grid(True)
plt.xlabel('$r$ (m)') 
plt.ylabel('$g(r)$')
plt.title("Kripton")

plt.suptitle("Gaseous State")
plt.show()


# U
#solid

Nblocks, U, errU = np.loadtxt("pot_solid_MC.dat", usecols=(0,2,3), unpack=True)
Nblocks_md, U_md, errU_md = np.loadtxt("solidU.out", usecols=(0,2,3), unpack=True)

U_Ar = U*epsAr
errU_Ar = errU*epsAr
U_Ar_md = U_md*epsAr
errU_Ar_md = errU_md*epsAr
U_Kr = U*epsKr
errU_Kr = errU*epsKr
U_Kr_md = U_md*epsKr
errU_Kr_md = errU_md*epsKr

plt.figure(figsize=(12,6))

plt.subplot(121)
plt.errorbar(Nblocks ,U_Ar, yerr=errU_Ar, label='MC')
plt.errorbar(Nblocks_md ,U_Ar_md, yerr=errU_Ar_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,U_Kr, yerr=errU_Kr)
plt.errorbar(Nblocks_md ,U_Kr_md, yerr=errU_Kr_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Potential energy per particle')
plt.show()

#liquid

Nblocks, U, errU = np.loadtxt("pot_liquid_MC.dat", usecols=(0,2,3), unpack=True)
Nblocks_md, U_md, errU_md = np.loadtxt("liquidU.out", usecols=(0,2,3), unpack=True)

U_Ar = U*epsAr
errU_Ar = errU*epsAr
U_Ar_md = U_md*epsAr
errU_Ar_md = errU_md*epsAr
U_Kr = U*epsKr
errU_Kr = errU*epsKr
U_Kr_md = U_md*epsKr
errU_Kr_md = errU_md*epsKr

plt.figure(figsize=(12,6))

plt.subplot(121)
plt.errorbar(Nblocks ,U_Ar, yerr=errU_Ar, label='MC')
plt.errorbar(Nblocks_md ,U_Ar_md, yerr=errU_Ar_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,U_Kr, yerr=errU_Kr)
plt.errorbar(Nblocks_md ,U_Kr_md, yerr=errU_Kr_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Potential energy per particle')
plt.show()


#gas

Nblocks, U, errU = np.loadtxt("pot_gas_MC.dat", usecols=(0,2,3), unpack=True)
Nblocks_md, U_md, errU_md = np.loadtxt("gasU.out", usecols=(0,2,3), unpack=True)

U_Ar = U*epsAr
errU_Ar = errU*epsAr
U_Ar_md = U_md*epsAr
errU_Ar_md = errU_md*epsAr
U_Kr = U*epsKr
errU_Kr = errU*epsKr
U_Kr_md = U_md*epsKr
errU_Kr_md = errU_md*epsKr

plt.figure(figsize=(12,6))

plt.subplot(121)
plt.errorbar(Nblocks ,U_Ar, yerr=errU_Ar, label='MC')
plt.errorbar(Nblocks_md ,U_Ar_md, yerr=errU_Ar_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,U_Kr, yerr=errU_Kr)
plt.errorbar(Nblocks_md ,U_Kr_md, yerr=errU_Kr_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$U/N$ (J)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Potential energy per particle')
plt.show()



# P
#solid

Nblocks, P, errP = np.loadtxt("pre_solid_MC.dat", usecols=(0,2,3), unpack=True)
Nblocks_md, P_md, errP_md = np.loadtxt("solidP.out", usecols=(0,2,3), unpack=True)

P_Ar = P*epsAr/sigAr**3
errP_Ar = errP*epsAr/sigAr**3
P_Ar_md = P_md*epsAr/sigAr**3
errP_Ar_md = errP_md*epsAr/sigAr**3

P_Kr = P*epsKr/sigKr**3
errP_Kr = errP*epsKr/sigKr**3
P_Kr_md = P_md*epsKr/sigKr**3
errP_Kr_md = errP_md*epsKr/sigKr**3

plt.figure(figsize=(12,6))

plt.subplot(121)
plt.errorbar(Nblocks ,P_Ar, yerr=errP_Ar, label='MC')
plt.errorbar(Nblocks_md ,P_Ar_md, yerr=errP_Ar_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$P$ (Pa)')
plt.grid(True)
plt.title('Argon')

plt.subplot(122)
plt.errorbar(Nblocks ,P_Kr, yerr=errP_Kr)
plt.errorbar(Nblocks_md ,P_Kr_md, yerr=errP_Kr_md, label='MD')
plt.legend()
plt.xlabel('#blocks')
plt.ylabel('$P$ (Pa)')
plt.grid(True)
plt.title('Kripton')


plt.suptitle('Pressure')
plt.show()

















