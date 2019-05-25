import matplotlib as mpl
import matplotlib.pylab as plt
from numpy import*
mpl.rc('text', usetex=False)
d = loadtxt("diagonals.dat")

plt.plot(d,'.')

#plt.ylim([0.5,1.5])

plt.show()
