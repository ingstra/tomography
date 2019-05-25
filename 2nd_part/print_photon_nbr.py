import matplotlib.pylab as plt
from numpy import*
from qutip import *
from matplotlib import cm

N=500
r = loadtxt("reconst_state_N500_maxlik300.dat",skiprows=1)
rho = r[:,0] + 1j*r[:,1]
rho = reshape(rho, (N,N))

rho_q = Qobj(rho)

n_op = num(N)


n_avg = expect(n_op,rho_q)
print("nbr of photons",n_avg)
