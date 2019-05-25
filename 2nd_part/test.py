from numpy import *
from math import factorial
from qutip import*


N=3
n=1

#rho=thermal_dm(N,n,method='analytic')
rho = array([[1,0,0],[0,1,0],[0,0,1]])/3
rho = Qobj(rho)
#print(rho)

grid_spacing =3.081013e-01

def proj(xi,pj,N):
    alpha = xi + 1j*pj
    D = displace(N,alpha)
    #print(D)
    #print(rho*D.dag())
    #print(D*rho*D.dag())
    Pi = D*rho*D.dag()*grid_spacing**2/pi
    return Pi

x =12.185949
p = 12.185949
#print(displace(3,0.25))


Pi=proj(x,p,N)

print(Pi)
