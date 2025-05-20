# Code description

Code computes the potential (current-free) magnetic field on a rectilinear mesh in a Cartesian box, given
the normal component of the magnetic field on each boundary face. The solution is computed using geometric multrigrid applied
to a finite-difference scheme. The method computes the magnetic field via a vector potential formulation.  The code returns both the vector 
potential in Coulomb gauge and the corresponding magnetic field.

The code uses a second-order finite-difference scheme for the discretization. In principle, this 
means that the numerical truncation error should decrease as the square of the mesh spacing. 

The backbone of the code is a set of modules for solving Poisson's equation in N dimensions
using geometric multigrid. This multigrid solver was written first, and then the vector potential
code was added later.The name of the code NDSM is derived from the original set of modules, i.e.
N-Dimensional Solver Multigrid (NDSM). The vector-potential module, however, is specifically 
designed for 3D, but leverages the more general N-dimensional backend. 

Code is written in Fortran 2003 and tested using the gfortran 8.3.0 compiler. It has only been
tested on a Linux platform. 

An earlier version of this code was used and is decribed in the paper Yang K.E., Wheatland M.S., and Gilchrist S.A.: 2020, ApJ,984,151. 
Paper DOI: 10.3847/1538-4357/ab8810

# Mesh and dimensions

The mesh is rectilinear, i.e. it is described by three mesh vectors x,y,z. These are assumed
to have fixed spacing. The code makes no explicit assumptions about the units of either B or A, although the
length scales are non-dimensional. 

# Vector Potential Gauge 

The vector potential is computed in the Coulomb gauge. However, note that in a box,
the Coulomb gauge is not necessarily unique when the boundary conditions are on the normal component
of the magnetic field. The NDSM code makes a particular choice in resolving this ambiguity. See
the paper and the notes for more details. 

# Convergence 

The multigrid method arrives at a solution via iteration. 
When convergence is poor, the output of the code may not accurately represent the 
solution to the underlying boundary-value problem (see ndsm_notes.pdf for details of the BVP). 
This has several consequences. Firstly, the magnetic field may not be a potential field, and significant electric 
currents may exist within the volume. Secondly, the normal component of the magnetic 
field may not match the normal component specified as boundary conditions.

## Metric

Two metrics are available for measuring the convergence of the solution: 
the max or mean difference between iterations. By default, the max is used. 

The max is sensitive to failure of convergence at any point, and therefore may be inappropriate
for many practical problems, but is useful for testing. The mean is a 
more robust convergence metric and may be more appropriate for practical problems. 

Setting mean=True, will use the mean rather than the max. 

## Tolerances

The code has two tolerance parameters that determine when to stop
iterating. 

### vc\_tol

The V cycle iteration stops when the max/mean difference is less than
vc\_tol. If vc\_tol is large, the code may return quickly, but the solution
may not be an accurate solution of the BVP.

### ex\_tol

The multrigrid method requires solution of a BVP on the coarsest mesh. 
This is solved via relaxation. The relaxation stops when the change in 
solution is less than ex_tol. Setting this to a large value will
result in inefficient V cycles, because the BVP is not being accurately
solved at each V cycle iteration.

### Choice of vc\_tol and ex\_tol

When testing the code on analytic solutions, both vc\_tol and ex\_tol
can be set to very small values. The default values defined in ndsm.py 
reflect values used for testing.

For some practical problems, the change in solution between iterations
may never reach the desired value of vc_tol: the solution is not improving
with additional V cycles. In this case, the iteration
will run until ncycles\_max is reached. This may take a long time 
depending on how ncycles\_max is chosen. A warning will be printed 
if the code returns without achieving vc\_tol. 

Setting a large value for vc\_tol (and ex\_tol) may prevent the 
code from running to ncycles\_max, but a large value of vc\_tol
in particular will mean the solutions is poorly converged: the numerical
solution is not an accurate solution of the underlying boundary-value problem.


# Compile shared library

The core Fortran code builds a shared library. 

Running make will build the shared library, called ndsm.so by default.

# OpenMP

The code is parallelized using the OpenMP standard. However, it should compile and run without OpenMP,
it will just be very slow on a multicore machine. 

# REAL and INTEGER Types

The core Fortran modules are written with a real type defined in NDSM_ROOT as REAL(FP). By default
this is set to C_DOUBLE. This can be changed to any supported Fortran real type without
breaking anything in the Fortran modules, however the Python interface only works with a real type
that is intercompatible with C_DOUBLE. 

Similarly, the basic integer type used throughout the code is INTEGER(IT), with IT = C_INT64_T. 
This again can be changed without resulting compiler errors. However, making the int size too small
may lead to overflow if large meshes are used, since the total number of mesh points is stored as a signed
Fortran integer. In addition, changing IT will break the Python wrapper. 

# Python 

The shared library can be called via the ndsm.py module. The module calls the subroutines
in the shared library using the Python ctypes module. The shared library needs to be compiled
first and either exist in sys.path, or else the explicit path to the shared library needs to
be passed as an argument to the function (see the docstring). 

The basic Python module only requires numpy and ctypes. Some of the tests require more 
modules, e.g. matplotlib. 

# Tests 

The repository contains code for running a number of integration and unit tests. Some are 
written in Fortran, while others are written in Python. The main integration test is designed
to demonstrate that the truncation error has the correct scaling with mesh spacing. This is 
a basic test of correctness for the method. 

The truncation error is estimated by applying the code to a known analytic test case and computing metrics
for the difference between the numerical and analytic solutions. The error metrics used are the max. 
and mean magnitude of the difference between the numerical and analytic vector fields. For a correctly 
implemented second-order scheme, (generally) both these metrics should decrease with the square of the mesh spacing
(for a uniform mesh). The max. error in particular may not achieve second order scaling for certain problems. 

A more complete description of the testing and results is included in the notes. 




