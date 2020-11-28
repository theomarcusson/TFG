from math import sqrt, erfc, exp, pi

#Define: cB [mg/L] (from "Xerrada_Porto"), D [m**2/s] (r=6nm) (from "Diffusion"), vm [m/s] (from "Diffusion")

"""
Calculating diffussion and velocity

cb = 5
di = (1.38e-23*300)/(6*3.1416*0.8513e-3*6e-9)
vm = 0.16e-6

print ('cB=%g [mg/L], D=%g m²/s, vm=%g m/s' % (cb, di, vm))
"""

cB = float(5)
D = float(4.3e-11)
v = float(1.6e-7)
T = float((4*D)/(pi*v**2))

#T = 2138.6445477973434

#Define J(t) = cB * v * (1 - 0.5 * erfc( sqrt( t / (pi*T) ) ) + 0.5 * sqrt(T/t) * exp( - t / (pi*T) ) )

def J(t):
 J = cB*v*(1 - 0.5 * erfc( sqrt( t / (pi*T) ) ) + 0.5 * sqrt(T/t) * exp( - t / (pi*T) ) )
 return J

#Define time

dt = 10 
t = 0.000000000000001
while t <= T:
 print ('time %g flow %g' % (t, J(t)))
 t *= dt
 
print (J(t))
