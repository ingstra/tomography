import matplotlib.pylab as plt
from numpy import*

m, d = loadtxt("frobenius.dat").T

plt.semilogy(m,d,'.-')

plt.title("Frobenius norm")


plt.show()
