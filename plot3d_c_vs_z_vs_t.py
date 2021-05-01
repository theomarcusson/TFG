from mpl_toolkits import mplot3d
from numpy import *
import matplotlib.pyplot as plt
import scipy.special as ss

#cB*(1-exp(-(v*z)/D)) - (cB/2) * ((erfc((z+v*t)/(sqrt(4*D*t)))) - exp(-((v*z)/D)) * (2 - (erfc((z-v*t)/(sqrt(4*D*t))))))

cB = float(5)
D = float(4.3e-11)
v = float(1.6e-7)
T = float((4*D)/(pi*v**2))

def c(t, z):
    return cB*(1-exp(-(v*z)/D)) - (cB/2) * ((ss.erfc((z+v*t)/(sqrt(4*D*t)))) - exp(-((v*z)/D)) * (2 - (ss.erfc((z-v*t)/(sqrt(4*D*t))))))
                                            
t = linspace(1e-15, 4000, 300)
z = linspace(0, 1e-3, 300)

X, Y = meshgrid(t, z)
Z = c(X, Y)
fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2, cmap='Blues', linewidth=1, antialiased=False)
ax.contour3D(X, Y, Z,300, cmap='binary')
ax.set_xlabel('t')
ax.set_ylabel('z')
ax.set_zlabel('c')

plt.show()
