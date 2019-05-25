import matplotlib.pylab as plt
from numpy import*


d = loadtxt("reconst_diag.dat")

t = loadtxt("thermal.dat")

plt.plot(d,'.-')

plt.plot(t,'-',linewidth=3)
print(sum(d))

plt.ylim([0,0.05])



plt.show()
