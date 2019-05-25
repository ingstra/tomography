import matplotlib.pylab as plt
from numpy import*
from qutip import *
from matplotlib import cm
import matplotlib as mpl
mpl.rc('text', usetex=False)
N=500
r = loadtxt("reconst_state.dat")
rho = r[:,0] + 1j*r[:,1]

rho = reshape(rho, (N,N))
print(shape(rho))
rho_q = Qobj(rho)

n_op = num(N)
n_avg = expect(n_op,rho_q)
print("nbr of photons",n_avg)


x=linspace(-5,5,500)
p=x

W = wigner(rho_q, x, p)

fig = plt.figure()
ax2 = fig.add_subplot(111)
plt2 = ax2.contourf(x, p, W, 100, cmap=cm.coolwarm)
cb1 = fig.colorbar(plt2, ax=ax2)

#plt.savefig("wigner300_G.png")
#print("saved fig")


plt.show()
