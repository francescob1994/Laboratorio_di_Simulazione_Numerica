import numpy as np
import matplotlib
import matplotlib.pyplot as plt


MU, SIGMA = np.loadtxt("opt_parameters.txt", usecols=(0,1), unpack = 'true' )


def psif(x):
    return ( np.exp(-((x+MU)**2)/(2*SIGMA**2)) + np.exp(-((x-MU)**2)/(2*SIGMA**2)) ) 

def Vpot(x):
    return (x**2 - 2.5)*x**2

norm = 2*SIGMA*np.sqrt(np.pi)*( 1 + np.exp(-(MU**2)/(SIGMA**2)) );
x = np.arange(-2,2, 0.02)
Psi = psif(x)
psi2 = (Psi**2)/norm
V = Vpot(x)


# 1) per i MU e SIGMA trovati stampo il datablock , max_rows = 1000

blk, progH, errH = np.loadtxt("output.H.0", usecols=(0,2,3), unpack='true') # skiprows = ind*1000, max_rows = 1000, inutile perché risimulo per avere x

plt.errorbar(blk, progH, yerr = errH)
plt.xlabel('block')
plt.ylabel('<H>')
plt.grid(True)
plt.title('Variational Montecarlo: MU = '+str(MU)+', SIGMA = '+str(SIGMA))
plt.show()


# 2) stampo l'istogramma di x e lo confronto la psi2 e con quella calcolata risolvendo l'eqz shrodinger numericamente

sample_x = np.loadtxt("x.dat", usecols=(0), unpack = 'true' )
n_bins = 50
n, bins, patches = plt.hist(sample_x, n_bins, density = True, range=(-5,5) )

plt.plot(x,psi2,label = "psi2")

plt.xlabel('x')
plt.ylabel('$|\psi (x)|^2$')
plt.title('Histogram of x')

plt.legend()
plt.grid(True)
plt.show()


# 3) soluzione numerica equazione di Schrodinger
  
a = 10
N = 100 # number of iterations
hbar = 1
m = 1

# Step sizes
x = np.linspace(-a/2, a/2, N)
dx = x[1] - x[0] # the step size
V = Vpot(x)

# The central differences method: f" = (f_1 - 2*f_0 + f_-1)/dx^2


CDiff = np.diag(np.ones(N-1),-1)-2*np.diag(np.ones(N),0)+np.diag(np.ones(N-1),1)
# np.diag(np.array,k) construct a "diagonal" matrix using the np.array
# The default is k=0. Use k>0 for diagonals above the main diagonal, 
# and k<0 for diagonals below the main diagonal

# Hamiltonian matrix
H = (-(hbar**2)*CDiff)/(2*m*dx**2) + np.diag(V)


# Compute eigenvectors and their eigenvalues
E,psi = np.linalg.eigh(H) #chiarire come è fatta psi, se agli autovettori stanno su righe o su colonne!!!!!!!????????????????????

# Take the transpose & normalize
psi = np.transpose(psi) 
psi = psi/np.sqrt(dx)       #perché?????????????????????? Forse perchè la prob è psi2*dx. Chiedi a Seb se è giusto!!!

print("Ground state energy: ", E[0])
print("1st excited state energy: ", E[1])
print("2nd excited state energy: ", E[2])

# Plot a few things
plt.figure(figsize=(8,5))
scale = 0.3
plt.plot(x, scale*V, color="Black", label="Potential") # plot the potential
plt.plot(x,(psi[0])**2, label="GS")
plt.plot(x,(psi[1])**2, label="1excited")
plt.plot(x,(psi[2])**2, label="2excited")
plt.title("Potential & Probabilities")
plt.xlabel("x")
plt.legend()
plt.grid(True)
plt.xlim((-3,3))
plt.ylim((-0.6,0.6))
plt.show()


# 4) confronto soluzione di psi2 con VMC, differenze finite e istogramma 


plt.hist(sample_x, n_bins, density = True, range=(-5,5), facecolor='y', label = 'VMC x sample' )

PSI = psif(x)
psi2 = (PSI**2)/norm

plt.plot(x,psi2,label = "psi2 VMC")
plt.plot(x,(psi[0])**2, label = "psi2 FD" )
plt.title("Variational Montecarlo & Finite Differencies")
plt.xlabel('x')
plt.legend()
plt.grid(True)
plt.show()


# confronto con QMC_1D

xq, psi2q, errpsiq = np.loadtxt("probability_psi1.dat", usecols=(0,1,2), unpack = 'true' )
xq_ms, psi2q_ms, errpsiq_ms = np.loadtxt("probability_psiMS.dat", usecols=(0,1,2), unpack = 'true' )

plt.errorbar(xq,psi2q, yerr=errpsiq, label='PIGS PSI_T=1')
plt.errorbar(xq_ms,psi2q_ms, yerr=errpsiq_ms, label='PIGS PSI_T=PSI(MU,SIGMA)')
vmc = psif(xq)
vmc2 = (vmc**2)/norm
plt.plot(xq, vmc2, label='VMC' )
plt.plot(x,(psi[0])**2, label='Finite Differences')
plt.legend()
plt.grid(True)
plt.show()


#diverse T

x, psi2_1, err_1 = np.loadtxt("prob_T0.5.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_2, err_2 = np.loadtxt("prob_T1.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_3, err_3 = np.loadtxt("prob_T1.5.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_4, err_4 = np.loadtxt("prob_T2.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_5, err_5 = np.loadtxt("prob_T2.5.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_6, err_6 = np.loadtxt("prob_T3.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_7, err_7 = np.loadtxt("prob_T3.5.dat", usecols=(0,1,2), unpack = 'true' )
x, psi2_8, err_8 = np.loadtxt("prob_T4.dat", usecols=(0,1,2), unpack = 'true' )

#plt.errorbar(xq_ms,psi2q_ms, yerr=errpsiq_ms, label='T=0')
#plt.errorbar(x,psi2_1, yerr=err_1, label='T=0.5')
plt.errorbar(x,psi2_2, yerr=err_2, label='T=1')
#plt.errorbar(x,psi2_3, yerr=err_3, label='T=1.5')
plt.errorbar(x,psi2_4, yerr=err_4, label='T=2')
#plt.errorbar(x,psi2_5, yerr=err_5, label='T=2.5')
plt.errorbar(x,psi2_6, yerr=err_6, label='T=3')
#plt.errorbar(x,psi2_7, yerr=err_7, label='T=3.5')
plt.errorbar(x,psi2_8, yerr=err_8, label='T=4')
plt.xlabel('x')
plt.ylabel('$|\psi (x)|^2$')
plt.legend()
plt.title('PIMC')
plt.grid(True)
plt.show()



























