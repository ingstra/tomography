import matplotlib.pylab as plt
from numpy import*


d = loadtxt("reconst_diag.dat")

t = loadtxt("thermal.dat")

plt.plot(d,'-',linewidth=4,label='reconstructed')

plt.plot(t,'--',linewidth=4,label='ideal thermal state with 35 photons')

print(sum(d))

plt.ylim([0,0.05])
plt.legend()
plt.xlabel(r"diagonal elements of $\rho$")


plt.show()
