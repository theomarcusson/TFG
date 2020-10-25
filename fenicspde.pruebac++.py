from fenics import *
import matplotlib.pyplot as plt

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
# If not counterindicated, units are in SI

# Define diffusion D=kB*T/6*pi*η*r: T=300, kB=1,38065*10^-23, pi=3,1416, η=0,8513*10^-3, r=10 nm
D = 2.5812 * pow(10, -11)

# Define vs
v = pow(10, -8)

u_D = Expression("formula", degree=2)

"""
Aquí no me queda claro si debería escribir la ecuación diferencial o la solución analítica... A continuación te copio un fragmento del tutorial oficial de FEniCS. 

"formula":

The next step is to specify the boundary condition: u = u D on ∂Ω. This is done by: 

"bc = DirichletBC(V, u_D, boundary)"

where u_D is an expression defining the solution values on the boundary, and boundary is a function (or object) defining which points belong to the boundary. Boundary conditions of the type u = u D are known as Dirichlet conditions. For the present finite element method for the Poisson problem, they are also called essential boundary conditions, as they need to be imposed explicitly as part of the trial space (in contrast to being defined implicitly as part of the variational formulation). Naturally, the FEniCS class used to define Dirichlet boundary conditions is named DirichletBC. The variable u_D refers to an Expression object, which is used to represent a mathematical function. The typical construction is:

"u_D = Expression(formula, degree=1)"

where formula is a string containing a mathematical expression. The formula must be written with C++ syntax and is automatically turned into an efficient, compiled C++ function.

degree=2: "[...] To obtain optimal
(order of) accuracy in computations, it is usually a good choice to use the same degree as for the space V that is used for the trial and test functions. However, if an Expression is used to represent an exact solution which is used to evaluate the accuracy of a computed solution, a higher degree must be used for the expression (one or two degrees higher). [...] u_D = Expression(’1 + x[0]*x[0] + 2*x[1]*x[1]’, degree=2). We set the degree to 2 so that u_D may represent the exact quadratic
solution to our test problem."

Pongo degree=2 pero realmente no entiendo bien qué conlleva este parámetro. Si no computa bien lo podría cambiar.

"""

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
plot(u)
plot(mesh)

plt.show()

# Save solution to file in VTK format
vtkfile = File("poisson/solution.pvd")
vtkfile << u

# Compute error in L2 norm
error_L2 = errornorm(u_D, u, "L2")

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
print("error_L2 =", error_L2)
print("error_max =", error_max)

# Hold plot

