#! /usr/bin/env python3
#
from dolfin import *
from mshr import *

def diffusion ( my_grid, my_range ):

#*****************************************************************************80
#
## diffusion simulates a simple 2D diffusion problem.
#
#  Discussion:
#
#      uxx = f in Omega = the unit interval.
#      u(0) = 0
#      u(1) = 1
#
#    where
#      x' = x/L and t'= t*v_M/L
#      (D= diffusion coef, L=size of the cuvette)
#
#      u_exact = x
#
#      f(x) = 0
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    12 March 2021
#
#  Author:
#
#    John Burkardt
#    Modified by Jordi Faraudo & Theo Marcusson
#
#  Reference:
#
#    Anders Logg, Kent-Andre Mardal,
#    Lectures on the Finite Element Method.
#
#  Parameters:
#
#    Input, integer my_grid, the resolution on the unit interval; float my_range [0,1], the interval of resolution
#

  import matplotlib.pyplot as plt

#
#  Set the mesh.
#
# 
#  Define the region for the equation to be solved
#
  west = Constant(0.0)
  south = Constant(0.0)
  east = Constant(1.0)
  north = Constant(1.0)

  domain = Rectangle(Point(west, south), Point(east, north)) 
  
# Define subdomain generation for the boundary layer
  
  n = 5
  m = 6
  series = [n**-i for i in range(m)]
  rectang = [Rectangle(Point(west, south), Point(east, v)) for v in series]
  
  for (i, bl) in enumerate(rectang):
   domain.set_subdomain(1 + i, bl)
  
# Generate the mesh
  
  my_mesh = generate_mesh(domain, my_grid)
    
  print ( '' )
  print ( '  Rectangular domain with mesh n = %d' % ( my_grid ) )

  
#
#  Set the function space.
# 
  V = FunctionSpace ( my_mesh, 'CG', 1 )
#
#  Set the trial and test functions.
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
#
#  Set velocity vector and constants (mu) 
#
  c_B = Constant(5)	#[mg/L]
  L = Constant(0.05)
  D = Constant(4.3*pow(10, -11))
  v_M = Constant(6*pow(10, -7))
  w = as_vector((0.0, (1)))
  my_mu = Constant( D / (L * v_M) )
# my_mu = 3.58333*pow(10, -6) 
#
#  Set the source term.
#
  f = Constant ( 0.0 )
#
#  Define the bilinear form.
#
  F = my_mu*dot(grad(u),grad(v))*dx - dot(w,grad(u))*v*dx - f*v*dx
  
  a, L = lhs(F), rhs(F)
#
#  Define the exact solution.
#
  u_expr = Expression ( "x[1]", degree = 10 )
#
#  Define the boundary condition (comment Jordi)
# Take care not to use == as explained in Fenics manual but to allow for tolerance
#
  def boundary(x, on_boundary):
#   check whether the point is at the boundary of the domain
    if on_boundary:
#       apply boundary conditions at y=0 or y=1
        if near(x[1], south, DOLFIN_EPS) or near(x[1], north, DOLFIN_EPS):
            return True
        else:
            return False
    else:
        return False

# The boundary condition uses the exact expression that provides automatically u(0)=0 and u(1)=1
# but any other function that has this property will work (comment Jordi)
  bc = DirichletBC ( V, u_expr, boundary )
#
#  Solve
#
  uh = Function ( V )
#
#  Solve the system.
#
  solve ( a == L, uh, bc )
#
# Plot the mesh
#
  plot(my_mesh)   
  plt.xlim([0, my_range])
  plt.ylim([0, my_range])
  plt.title ( 'trial_mesh_2d, grid %d' % ( my_grid ) )
  filename = ( 'trial_mesh_2d_grid_%d_range_%g.png' % ( my_grid, my_range ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  print ( '  Max resolution in Boundary Layer: %g m' % (n**(-m+1)))

#
#  Plot the solution (concentration).
#
  fig = plt.figure()
  #2D plot with a color bar
  ax = plt.subplot ( 111 )
  cs = plot (uh)
  cbar = fig.colorbar(cs)
  plt.xlim([0, my_range])
  plt.ylim([0, my_range])
  filename = ( 'diffusion_solutions_grid%d_range_%g.png' % ( my_grid, my_range ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def diffusion_test ( ):

#*****************************************************************************80
#
## convection_diffusion_test tests convection_diffusion.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 November 2018
#
#  Author:
#
#    John Burkardt
#    Modified Jordi Faraudo & Theo Marcusson
#
  import dolfin
  import platform
  import time

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'Stationary Solution of convection diffusion Equation:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  FENICS version %s'% ( dolfin.__version__ ) )
  print ( '  Diffusion problem on the unit interval.' )
  print ( '  - ux = f in Omega = the unit interval.' )
  print ( '  u(0) = 0, u(1) = 1' )
  
#  
#
# Set grid
#  
  my_grid = 20

# Set range from y'(y=0)=0 to y'(y=L)=1; 
# my_range = [0, 1]
  
  for my_range in ( 0.01 , 0.1 , 1.0 ):
    diffusion ( my_grid , my_range )

#
#  Terminate.
#
  print ( '' )
  print ( 'diffusion_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  diffusion_test ( )
