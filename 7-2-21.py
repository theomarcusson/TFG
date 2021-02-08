#7 febrero 2021


from fenics import *
import numpy as np
import matplotlib.pyplot as plt

#Define time
T = 2.0		# final time
num_steps = 50		# number of time steps
dt = T / num_steps 	# time step size

#Insert mesh

nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)

V = FunctionSpace(mesh, 'P', 1)


#Boundary Condition c(z=0, t) = 0  -> Homogenous Dirichlet boundary condition

u_D = Constant(0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)


#Define constants
c_B = Constant(5)
D = Constant(4.3*pow(10, -11))
w = Constant(6*pow(10, -7))
#Define initial value: c(z, t=0) = c_B

u_0 = Constant(c_B)
#u_n = project(u_0, V)
#or
u_n = interpolate(u_0, V)

#Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx - D*dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx - w*dt*grad(u)*v*dx
#Not sure about the notation here                             (^                 ^)
#Error: Can only integrate scalar expressions. The integrand is a tensor expression with value shape (2,) and free indices with labels ().	(ufl.log.UFLException:)
#This is a similar error to the attempts before January 2021


a, L = lhs(F), rhs (F)

# Create VTK file for saving solution
vtkfile = File('heat/solution.pvd')

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):	
	
	# Update current time
	t += dt		
	
	# Compute solution
	solve(a == L, u, bc)	
	
	# Save to file and plot solution
	vtkfile << (u, t)
	plot(u)

	# Update previous solution
	u_n.assign(u)

# Hold plot
plt.show()

