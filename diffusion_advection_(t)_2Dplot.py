
from dolfin import *
from mshr import *

def diffusion_advection ( my_grid, my_range ):

#*****************************************************************************80
#
## diffusion_advection simulates a 2D diffusion-advection time-dependant problem.
#
#  Discussion:
#
#      mu*(d²u/dx²) - du/dx = du/dt + f(x) in Omega = the unit interval.
#      u(0) = 0
#      u(1) = 1
#
#    where
#      x = z/L and t = rt*v/L
#      (x=dimensionless magnetoph. direction, z=magnetoph. direction, L=size of the cuvette(m))
#      (t=dimensionless time, rt=real time (s), v=magnetoph. velocity (m/s))
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
#    Input: integer my_grid, the resolution on the mesh; 
#    Input: float my_range [0,1], the range of observation (= (1/zoom))
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
  
  n = 1.3
  m = 30
  series = [n**-i for i in range(m)]
  rectang = [Rectangle(Point(west, south), Point(east, v)) for v in series]
    
  for (i, bl) in enumerate(rectang):
   domain.set_subdomain(1 + i, bl)
  
# Generate the mesh
  
  my_mesh = generate_mesh(domain, my_grid)
    
  print ( '' )
  print ( '  Rectangular domain with mesh n = %d' % ( my_grid ) )
  #my_mesh = UnitIntervalMesh ( my_grid )
  
#
#  Set the function space.
# 
  V = FunctionSpace ( my_mesh, 'CG', 1 )
#
#  Set the trial and test functions.
#
  u = Function ( V )
  u1 = Function ( V )
  v = TestFunction ( V )
#
#  Set velocity vector and constants (mu) 
#
  c_B = Constant(5)				#[mg/L]
  L = Constant(0.05)				#[m]
  D = Constant(4.3*pow(10, -11))		#[m²/s]
  v_M = Constant(1.6*pow(10, -7))		#[m/s]
  w = as_vector((0.0, (1)))
  my_mu = Constant( D / (L * v_M) )
  
#     my_mu = 5.375*pow(10, -3) 
#     1/my_mu ~ 186 
#     t -> dimensionless time; rt -> real time, in seconds
#     t = rt * ( v_M / L )
#     t = rt * 3.2*e-6
#     t * 312500 = rt

#
#  Set the source term.
#
  f = Constant ( 0.0 )
#
#  Set the time-stepping
#
  T = (1)                            #Final dimensionless time
  RT = Constant( T * ( L / v_M ))    #Real time, in seconds
  num_steps = 312500                 #number of time steps
  dt = (T / num_steps)               #
  
#
#  Define the initial condition u(t=0) = 1
#  
  u_0 = Constant(1) 
  u1 = project(u_0, V)
#
#  Set the right hand and left hand side.
# 
  F = my_mu*dot(grad(u),grad(v))*dx - dot(w,grad(u))*v*dx + ((u - u1) / dt)*v*dx
#
  a, L = lhs(F), rhs(F)
#
#  Define the exact solution.
#
  u_expr = Expression ( "x[1]", degree = 10 )
#
#  Define the boundary condition
#
  def boundary(x, on_boundary):
#   check whether the point is at the boundary of the domain
    if on_boundary:
#       apply boundary conditions at x=0 or x=1
        if near(x[1], south, DOLFIN_EPS) or near(x[1], north, DOLFIN_EPS):
            return True
        else:
            return False
    else:
        return False
#  
#  
  print ( '  Max resolution in Boundary Layer: %g m' % (n**(-m+1)))
  
# The boundary condition uses the exact expression that provides automatically u(0)=0 and u(1)=1
# but any other function that has this property will work (comment Jordi)
  bc = DirichletBC ( V, u_expr, boundary )
#
#  Solve
#
  t = 0			#dimensionless time

  for j in range(num_steps):
    t+=dt
    solve(F == 0, u, bc)
    plt.xlim([0, my_range])
    plt.ylim([0, my_range])
    zoom = 1/my_range
    u1.assign(u)
    rt = t * 312500		#real time, in seconds.
    filename = ('difu_adv_zoom_%g_rt_%.1f_t_%.6f.png' % ( zoom, rt, t ))
    ro = round((t-dt)*num_steps)

# This is made to save results as a png file just once for a given number of iterations
# given by the number following the % sign in the if statement. ro is made to round to
# integers and make the algorithm work with no problems of tolerance.   
    
    if rt <= 60:
      if (ro) % 1 == 0:
        pu = plot(u)
        cb = plt.colorbar(pu)
        plt.savefig(filename)
        plt.close()
    elif 60 < rt <= 1000:
      if (ro) % 10 == 0:
        pu = plot(u)
        cb = plt.colorbar(pu)
        plt.savefig(filename)
        plt.close()
    elif 1000 < rt <= 20000:
      if (ro) % 100 == 0:
        pu = plot(u)
        cb = plt.colorbar(pu)
        plt.savefig(filename)
        plt.close()
    else:
      if (ro) % 10000 == 0:
        pu = plot(u)
        cb = plt.colorbar(pu)
        plt.savefig(filename)
        plt.close()

#
#  Terminate.
#
  return

def diffusion_advection_test ( ):

#*****************************************************************************80
#
## diffusion_advection_test tests diffusion_advection.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 May 2021
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
  print ( '  Diffusion-advection time dependant problem on 2D mesh.' )
  print ( '  mu*(d²u/dx²) - du/dx = du/dt + f(x)' )
  print ( '  u(0) = 0, u(1) = 1' )
  
#  
#
# Set grid
#  
  my_grid = 30

# Set range
# from x(i.e. z=0) = 0 to x(i.e. z=L) = 1; my_range = [0, 1]
  
  for my_range in ( 0.02 , 1.0 ):
    diffusion_advection ( my_grid , my_range )

#
#  Terminate.
#
  print ( '' )
  print ( 'diffusion_advection_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  diffusion_advection_test ( )
