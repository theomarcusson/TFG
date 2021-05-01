from math import sqrt, erfc, exp, pi
from numpy import *
from matplotlib.pyplot import *


#CAMBIANDO VELOCIDAD Y DIAMETRO

di = (1.38e-23*300)/(6*3.1416*0.8513e-3*6e-9) 
cB = float(5)
D = float(di)
v = float(1.6e-7)
T = float((4*D)/(pi*v**2))


#Variar t entre 1e-15 y 9e4

t0 = float (1e-15)
t1 = float (60)
t2 = float (15 * 60)
t3 = float (3600)
t4 = float (86400)
t5 = float (30 * 86400)

def q(z):
 return cB*(1-exp(-(v*z)/D)) - (cB/2) * ((erfc((z+v*t)/(sqrt(4*D*t)))) - exp(-((v*z)/D)) * (2 - (erfc((z-v*t)/(sqrt(4*D*t))))))

z = linspace(0, 1e-3, 400)
y0 = zeros(len(z))
y1 = zeros(len(z))
y2 = zeros(len(z))
y3 = zeros(len(z))
y4 = zeros(len(z))


for i in range (len(z)):
 t = t0
 y0[i] = q(z[i]) 
 t = t1
 y1[i] = q(z[i]) 
 t = t2
 y2[i] = q(z[i])    
 t = t3
 y3[i] = q(z[i]) 
 t = t4
 y4[i] = q(z[i])   


plot(z, y4, label = 't = 1 day')
plot(z, y3, label = 't = 1 h')
plot(z, y2, label = 't = 15 min')
plot(z, y1, label = 't = 1 min')
plot(z, y0, label = 't = 0 s')

xlabel ('z [m]')
ylabel('c(z,t) [mg/L]')
axis([0, 1e-3, -0.1, 5.1])
legend()
show()
    
