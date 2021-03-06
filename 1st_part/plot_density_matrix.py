import matplotlib.pylab as plt
from numpy import*
from qutip import *
from matplotlib import cm

N=50
r = loadtxt("reconst_state.dat",skiprows=1)
rho = r[:,0] + 1j*r[:,1]
rho = reshape(rho, (N,N))

rho = rho[0:20,0:20]

rho_q = Qobj(rho)

matrix_histogram(rho_q)
#hinton(rho_q)

plt.show()
