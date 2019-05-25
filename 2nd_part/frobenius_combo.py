import matplotlib.pylab as plt
from numpy import*

m, d = loadtxt("frobenius.dat").T
n,b = loadtxt("frobenius_N500_maxlik1600_init1.dat").T

m=append(n,m)
d=append(b,d)

print(m,n)
plt.semilogy(m,d,'.-')



plt.title("Frobenius norm")


plt.show()
