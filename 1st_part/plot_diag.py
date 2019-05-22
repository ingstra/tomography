import matplotlib.pylab as plt
from numpy import*

d = loadtxt("diagonals.dat")

plt.plot(d,'.')

#plt.ylim([0.5,1.5])

plt.show()
