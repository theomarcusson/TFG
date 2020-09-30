#solve: dy/dt = -k*y(t)
#init cond: y(0) = 5

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#the function dy/dt is defined as func

def func(y,t):
    k = 1.0
    dydt = -k*y
    return dydt

#init cond
y0 = 5

#time
t = np.linspace(0,10)

#solve EDO
y = odeint(func,y0,t)

#plot

plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()
