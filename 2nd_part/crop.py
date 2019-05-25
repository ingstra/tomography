from numpy import *

H = loadtxt("Bin143_start_-22.156_delta_0.312056_ON.txt")

size = 143
skip = 32
delta = 0.312056

x_start = 22.156
x_new = x_start - delta*skip
print(x_new)

h = H[skip:size-skip,skip:size-skip]

savetxt("Bin79_testON.dat",h, fmt="%f")

print(shape(h))
