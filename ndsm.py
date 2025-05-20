# BSD-3-Clause 
#
# Copyright 2024 S.A Gilchrist
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors 
# may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A !PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT !LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

#
# Python interface for NDSM library. Basically a wrapper for 
# the ndsm_vector_solve subroutine
# 

# Basic
import os 
import sys

# Numpy
import numpy as np
from numpy import ctypeslib as nct

# Ctypes
import ctypes 

def get_lib_path(libname):
  """
    Locate first instance of libname in system path

    Parameters:
    -----------
      libname: str
        Name of shared library

    Returns:
    --------
      libpaths: str
        Path to library libname if it is found, None otherwise. 	
  """
  libpaths = []
  for path in sys.path:
    for root,dirs,files in os.walk(path):
      if libname in files:
        libpaths.append(os.path.join(root,libname))
  # Cheap way to get unique values
  return list(set(libpaths))

# ---------------------------------------------------------------------

def vector_potential(x,y,z,b,niterex_max=10000,ncycles_max=1024,ex_tol=1e-13,vc_tol=1e-10,ms=5,mean=False,libname="ndsmf.so",libpath=None,debug=False):
  """

    For multigrid, it is necessary to compute an "exact" solution on the 
    coarsest grid/mesh. This solution is computed using the same 
    iterative relaxation method as for smoothing. This needs to be run
    to convergence for the method to work well. The ex_tol parameter 
    controls the convergence tolerance.
   

    Parameters:
    -----------
      x,y,z: (nx,),(ny,),(nz,)
        Mesh vectors
      b: (3,nz,ny,nx)
        Array to hold magnetic field. The function takes the normal
        component at the boundary as boundary conditions to compute 
        the vector potential and magnetic field. Interior values are not used.
      ex_tol: float
        Convergence tolerance used for solution on coarsest grid. Iteration
        on the coarsest grid stops when the max absolute difference between
        iterations is below this value. 
      vc_tol: float
        Convergence tolerance for V-cycles. V-cycling will halt 
        when this is reached and solution will be returned. 
      ms: int 
        Number of relaxation sweeps to perform before and after 
        interpolation/restriction. This should be a small number. 
      niterex_max: int
        Max. number of iterations used fot exact solution 
      ncycles_max: int
        Max. number of V cycles. Algorithm will stop at this many
        cycles and return a not converged flag
      libname: str
        Name of external Fortran shared library 
      libpath: str (optional)
        Path to folder containing shared library. If not set, then
        the sys.path is searched. An error will occur if either the library
        isn't found or if multiple versions are found.
      debug: bool
        Debug flag. If set, the Fortran code will print more 
        information 
      mean: bool
        Use mean difference as the metric to measure change in solution
        rather than the max

    Returns:
    --------
      ierr: int 
        Error flag. A value other than zero indicates a problem
        with convergence of the V cycle
      Apop:(3,nz,ny,nx):
        Vector potential 
      b: (3,nz,ny,nx)
        Magnetic field 
  """

  # Get path to .so
  if(libpath is None):
    libpath = get_lib_path(libname)
    if(len(libpath) is None):
      raise ValueError("Could not locate {:s}".format(libname))
    elif(len(libpath) > 1):
      print(libpath)
      raise ValueError("More than once instance of {:s} found\n {:s}".format(libname,str(libpath)))
    else:
      libpath = libpath[0]    

  # Try to load .so file
  try:
    libc = ctypes.cdll.LoadLibrary(libpath)
  except: 
    raise ValueError("Could not load library at "+libpath)

  # Arguments
  arg_nshape = nct.ndpointer(np.intc   ,ndim=1,flags=('C','A','W'))
  arg_ioptc  = nct.ndpointer(np.intc   ,ndim=1,flags=('C','A','W'))
  arg_ropt   = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))  
  arg_b      = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))
  arg_a      = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))
  arg_x      = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))
  arg_y      = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))
  arg_z      = nct.ndpointer(np.float64,ndim=1,flags=('C','A','W'))

  # Set arg and return types
  libc.ndsm_vector_solve.argtypes = [ctypes.c_size_t,arg_nshape,arg_ioptc,arg_ropt,
                                     arg_x,arg_y,arg_z,arg_a,arg_b]

  # Get length of options vector 
  IOPT_LEN = libc.get_iopt_len()

  # Get arguments
  #
  # Shape need to be in Fortran order: [nx,ny,nz,3]
  # 
  nshape = np.array(b.shape[::-1],dtype=np.intc)
  nsize  = ctypes.c_size_t(b.size)
  ioptc  = np.zeros(IOPT_LEN,dtype=np.intc)
  ropt   = np.zeros(IOPT_LEN,dtype=np.float64)
 
  # Get option flag positions. These correspond to positions in 
  # iopt and ropt. These indices are in the range [0,IOPT_LEN)
  iopt_ms      = libc.get_iopt_ms()
  iopt_ncycles = libc.get_iopt_ncycles()
  iopt_ctol    = libc.get_ropt_ctol()
  iopt_vtol    = libc.get_ropt_vtol()
  iopt_debug   = libc.get_iopt_debug()
  iopt_dumax   = libc.get_iopt_dumax()
  iopt_nmaxex  = libc.get_iopt_iopt_nmaxex()
  
  # Check bounds. This shouldn't occur in normal function. Only 
  # an bug in the Fortran lib should cause this
  indx = np.array([iopt_ms,iopt_ncycles,iopt_ctol,iopt_vtol],dtype=np.intc)
  if(np.any(indx<0) or np.any(indx>=IOPT_LEN)):
    raise Exception("Option vector (IOPT) index out of bounds. This shouldn't occur.")

  # Set options
  ioptc[iopt_ms]      = ms           # Smoothing steps
  ioptc[iopt_ncycles] = ncycles_max  # Max. V cycles
  ioptc[iopt_nmaxex]  = niterex_max  # Max. iterations for exact solution
  ropt[iopt_vtol]     = vc_tol       # V cycle tol
  ropt[iopt_ctol]     = ex_tol       # Tol. course

  # Debug flag
  if(debug):
    ioptc[iopt_debug] = libc.get_iopt_true()          
  else:
    ioptc[iopt_debug] = libc.get_iopt_false()     

  # If mean flag is set, then set dumax to false
  if(mean):
    ioptc[iopt_dumax] = libc.get_iopt_false()
  else:
    ioptc[iopt_dumax] = libc.get_iopt_true()
                

  # Data arrays
  Apot = np.zeros(b.size,dtype=np.float64)  
  b    = b.flatten()
  
  # Call external
  ierr = libc.ndsm_vector_solve(nsize,nshape,ioptc,ropt,x,y,z,Apot,b)

  # Return and reshape to [3,nz,ny,nx]
  return ierr,Apot.reshape(nshape[::-1]),b.reshape(nshape[::-1])

