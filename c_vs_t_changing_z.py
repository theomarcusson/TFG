from math import sqrt, erfc, exp, pi
from numpy import linspace, zeros
from matplotlib.pyplot import *
 
    
    

#CAMBIANDO VELOCIDAD Y DIAMETRO
    
di = (1.38e-23*300)/(6*3.1416*0.8513e-3*6e-9)   
cB = float(5)
D = float(di)
v = float(1.6e-7)
T = float((4*D)/(pi*v**2))

#Variar z  entre 0 y 1e-3

z0 = float (0)
z1 = float (5e-6)
z2 = float (5e-5)
z3 = float (1e-4)
z4 = float (3e-4)
z5 = float (1e-3)

def q(t):
 return cB*(1-exp(-(v*z)/D)) - (cB/2) * ((erfc((z+v*t)/(sqrt(4*D*t)))) - exp(-((v*z)/D)) * (2 - (erfc((z-v*t)/(sqrt(4*D*t))))))

t = linspace(1e-15, 3e3, 200)
y0 = zeros(len(t))
y1 = zeros(len(t))
y2 = zeros(len(t))
y3 = zeros(len(t))
y4 = zeros(len(t))
y5 = zeros(len(t))

for i in range (len(t)):
 z = z0
 y0[i] = q(t[i])
 z = z1
 y1[i] = q(t[i]) 
 z = z2
 y2[i] = q(t[i])
 z = z3
 y3[i] = q(t[i])   
 z = z4
 y4[i] = q(t[i])   
 z = z5
 y5[i] = q(t[i])   
    
plot(t, y0, label = 'z = 0 m')
plot(t, y1, label = 'z = 5e-6 m')
plot(t, y2, label = 'z = 5e-5 m')
plot(t, y3, label = 'z = 1e-4 m')
plot(t, y4, label = 'z = 3e-4 m')
plot(t, y5, label = 'z = 1e-3 m')
xlabel ('t [s]')
ylabel('c(z,t) [mg/L]')
legend(['c(z,t) [mg/L] vs t [s]'])
axis([0, 3000, -0.1, 5.1])
legend()
show()

