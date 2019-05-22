import matplotlib.pylab as plt
from numpy import*
from qutip import *
from matplotlib import cm

N=800
r = loadtxt("reconst_state.dat",skiprows=1)
rho = r[:,0] + 1j*r[:,1]

rho = reshape(rho, (N,N))
print(shape(rho))
rho_q = Qobj(rho)


x=linspace(-10,10,500)
p=x

W = wigner(rho_q, x, p)

fig = plt.figure()
ax2 = fig.add_subplot(111)
plt2 = ax2.contourf(x, p, W, 100, cmap=cm.coolwarm)
cb1 = fig.colorbar(plt2, ax=ax2)

plt.savefig("wigner.png")
print("saved fig")

#plt.show()
