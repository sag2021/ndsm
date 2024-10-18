# Code description

Code computes the potential (current-free) magnetic field in a box with a rectilinear mesh using 
geometric multrigrid given the normal component of the magnetic field on the boundary faces of the box. 
The method computes the potential field via a vector potential formulation.  The code returns both the vector 
potential in Coulomb gauge and the corresponding magnetic field.

The code uses a second-order finite-difference scheme for the discretization. In principle, this 
means that the numerical truncation error should decrease as the square of the mesh spacing. 

The backbone of the code is a set of module for solving Poisson's equation in N dimensions
using geometric multigrid. This multigrid solver was written first, and then the vector potential
code was added later.Thename of the code NDSM is derived from the original set of modules, i.e.
N-Dimension Solver Multrigrid (NDSM). The vector-potential module, however, is specifically 
designed for 3D, but leverages the more general N-dimensional backend. 

Code is written in Fortran 2003 and tested using the gfortran 8.3.0 compiler. It has only been
tested on a Linux platform. 

An earlier version of this code was used and is decribed in the paper Yang K.E., Wheatland M.S., and Gilchrist S.A.: 2020, ApJ,984,151. 
Paper DOI: 10.3847/1538-4357/ab8810

# Mesh and dimensions

The mesh is rectilinear, i.e. it is described by three mesh vectors x,y,z. These are assumed
to have fixed spacing. 

The code makes no explicit assumptions about the units of either B or A, although the
length scales are non-dimensional. 

# Vector Potential Gauge 

The vector potential is computed in the Coulomb gauge. However, note that in a box,
the Coulomb gauge is not necessarily unique, when the boundary conditions are on the normal component
of the magnetic field. The NDSM code makes a particular choice in resolving this ambiguity. See
the paper and the notes for more detail. 

# Compile shared library

The core Fortran code builds a shared library. 

Running make will build the shared library, called ndsm.so by default.

Without make, the code be be compiled as 

gfortran -O3 -fpic -shared -o ndsm.so ndsm_root.f90  ndsm_interp.f90 ndsm_multigrid_core.f90 ndsm_optimized.f90 ndsm_poisson.f90 ndsm_vector_potential.f90

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

The basic Python module only requires numpy and ctypes. The integration and unit tests require
more modules, e.g. matplotlib. 

# Tests 

The repository contains code for running a number of integration and unit tests. Some are 
written in Fortran, while others are written in Python. The main integration test is designed
to demonstrate that the truncation error has the correct scaling with mesh spacing. This is 
a basic test of correctness for the method. 

The truncation error is estimated by applying the code to a known analytic test case and computing metrics
for the difference between the numerical and analytic solutions. The error metrics used as the max. 
and mean magnitude of the difference between the numerical and analytic vector fields. For a correctly 
implemented second-order scheme, (generally) both these metrics should decrease with the square of the mesh spacing
(for a uniform mesh). The max. error in particular may not achieve second order scaling for certain problems. 

A more complete description of the testing and results is included in the notes. 




