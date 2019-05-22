from qutip import *
from numpy import *
from scipy import sparse
import matplotlib.pyplot as plt

N=10
nbr_photons=0.1

rho = thermal_dm(N,nbr_photons)

nbr_grid_pts = 200

xlim = 2

x = linspace(-xlim, xlim, nbr_grid_pts)
p=x

X, P = meshgrid(x,p)

grid = X + 1j*P

disp = vectorize(displace)

D = disp(N, grid)

D_dag = [d.dag() for idx,d in ndenumerate(D)]
D_dag = reshape(D_dag,(nbr_grid_pts,nbr_grid_pts))

Pi = D*rho*D_dag/pi

sum_Pi = Pi.sum()

hinton(sum_Pi)
print(sum_Pi.data)
#plt.show()
