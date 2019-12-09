import numpy as np
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py number")
    sys.exit(2)

nome_script, n = sys.argv #n contiene il numero del file
u, p = np.loadtxt("u_p."+str(n), usecols=(0,1), unpack = 'true')

cov_u = np.zeros(len(u)-1)  
cov_p = np.zeros(len(u)-1)  

U = u[0]
P = p[0]

for i in range(len(u)-1):  
   cov_u[i] = U*u[i+1]
   cov_p[i] = P*p[i+1]

file = open("cov."+str(n), "w")
for i in range(len(u)-1):
   file.write(str(cov_u[i])+"      "+str(cov_p[i])+"\n")

file.close()








