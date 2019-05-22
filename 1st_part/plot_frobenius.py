import matplotlib.pylab as plt
from numpy import*

m, d = loadtxt("frobenius_fakeOFF48.dat").T

plt.plot(m,d,'.-')

plt.title("Frobenius norm")
plt.xlim([-1,200])

plt.show()
