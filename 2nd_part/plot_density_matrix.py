import matplotlib.pylab as plt
from numpy import*
from qutip import *
from matplotlib import cm

N=800
r = loadtxt("reconst_state_N800_maxlik200.dat")
rho = r[:,0] + 1j*r[:,1]
rho = reshape(rho, (N,N))
N=10
rho = rho[0:N,0:N]
print(shape(rho))
rho_q = Qobj(rho)


rho_neg=rho.diagonal()[rho.diagonal()<0]
print(rho_neg)
print(rho[0,3],rho[3,0])

print("trace rho",trace(rho))

#matrix_histogram(rho_q)
hinton(rho_q)

plt.show()
