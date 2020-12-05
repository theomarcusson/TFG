#!/usr/bin/env python
# coding: utf-8

from numpy import *
import scipy.special as ss

#Define c(z,t) = cB*(1-exp(-vs*z/D)) - cB/2 * (erfc((z+vs*t)/(sqrt(4*D*t)))) - exp(-(vs*z/D)) * (2 - (erfc((z-vs*t)/(sqrt(4*D*t)))))

#Define c0(z) = cB*(1-exp(-vs*z/D)) 

"""
Define: cB [mg/L] (from "Xerrada_Porto"), D [m**2/s] (r=6nm) (from "Diffusion"), vm [m/s] (from "Diffusion")

cb = 5
di = (1.38e-23*300)/(6*3.1416*0.8513e-3*6e-9)
vm = 0.16e-6

print ('cB=%g [mg/L], D=%g m²/s, vm=%g m/s' % (cb, di, vm))
"""


cB = float(5)
D = float(4.3e-11)
v = float(1.6e-7)
T = float((4*D)/(pi*v**2))

def c(z, t):
 return cB*(1-exp(-(v*z)/D)) - (cB/2) * (ss.erfc((z+v*t)/(sqrt(4*D*t)))) - exp(-((v*z)/D)) * (2 - (ss.erfc((z-v*t)/(sqrt(4*D*t)))))

dt = 10 
t = 1e-15
while t <= T:
 t *= dt
 dz = 10
 z = 1e-9
 while z <= 0.1:
        z *= dz
        print ('time %g     z %g      concentration %g' % (t, z, c(z, t)))


