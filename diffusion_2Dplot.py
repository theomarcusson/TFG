#! /usr/bin/env python3
#
from dolfin import *
from mshr import *

def diffusion ( my_grid ):

#*****************************************************************************80
#
## diffusion simulates a simple 1D diffusion problem.
#
#  Discussion:
#
#      uxx = f in Omega = the unit interval.
#      u(0) = 0
#      u(1) = 1
#
#    where
#      x' = x/L and t'= Dt/L*L
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
#    Modified by Jordi Faraudo
#
#  Reference:
#
#    Anders Logg, Kent-Andre Mardal,
#    Lectures on the Finite Element Method.
#
#  Parameters:
#
#    Input, integer my_grid, the resolution on the unit interval.
#
  #import matplotlib as mpl
  import matplotlib.pyplot as plt
  #from mpl_toolkits.mplot3d import Axes3D
  #from matplotlib import ticker, cm

#
#  Set the mesh.
#
#
#  Define the region for the equation to be solved
#
  sw = Point ( 0.0, 0.0 )
  ne = Point ( 1.0, 1.0 )

#  domain = Rectangle ( sw, ne )
#
#  Mesh the region.
#
#  grid_x = my_grid
#  grid_y = my_grid
  my_mesh = RectangleMesh ( sw, ne, my_grid ,my_grid)

  print ( '' )
  print ( '  Rectangualr domain with mesh n = %d' % ( my_grid ) )
  #my_mesh = UnitIntervalMesh ( my_grid )
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
#  Set the right hand side.
#
  f = Constant ( 0.0 )
#
#  Define the bilinear form and right hand side
#
  #a = (  u.dx(0) * v.dx(0) ) * dx
  a = dot( grad(u) , grad(v) ) * dx
  L = f * v * dx
#
#  Define the exact solution.
#
  u_expr = Expression ( "x[0]", degree = 10 )
#
#  Define the boundary condition (comment Jordi)
# Take care not to use == as explained in Fenics manual but to allow for tolerance
#
#  def boundary ( x ):
#   Boundaries are x=0 and x=1.0
#    value = x[0] < DOLFIN_EPS or 1.0 - DOLFIN_EPS < x[0]
#    return value

  def boundary(x, on_boundary):
#   check whether the point is at the boundary of the domain
    if on_boundary:
#       apply boundary conditions at x=0 or x=1
        if near(x[0], 0, DOLFIN_EPS) or near(x[0], 1, DOLFIN_EPS):
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
#  Project the exact solution.
#
  u_exact = interpolate ( u_expr, V )
  flux_expr = Constant ( 1.0 )
  flux_exact = interpolate (flux_expr,V)

# Calculate the flux of particles which is even harder than the concentration
# by Jordi Faraudo
  #diffusive flux in dimensionless units
  flux = project(uh.dx(0),FunctionSpace(my_mesh, 'CG', 1))

#
#  Plot the solution (concentration).
#
  fig = plt.figure()
  #2D plot with a color bar
  ax = plt.subplot ( 111 )
  cs = plot (uh)
  #plot (u_exact, label = 'Exact' )
  #ax.legend ( )
  ax.grid ( True )
  cbar = fig.colorbar(cs)
  plt.title ( 'diffusion solutions, grid %d' % ( my_grid ) )
  filename = ( 'diffusion_solutions_grid%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )

#
#  Plot also the flux of particles 
#
  fig = plt.figure ( )
  ax = plt.subplot ( 111 )
  #ax.set(xlim=(0.0, 1.0), ylim=(-1.5, 1.5))
  cs=plot (flux, label = 'Computed' )
  #plot (flux_exact, label = 'Exact (1.0)' )
  ax.legend ( )
  ax.grid ( True )
  cbar = fig.colorbar(cs)
  plt.title ( 'Flux, grid %d' % ( my_grid ) )
  filename = ( 'Flux_grid%d.png' % ( my_grid ) )
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
#    Modified Jordi Faraudo 
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

  for my_grid in ( 10,100):
    diffusion ( my_grid )
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
