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

#
# Integration test 1: 
#
# Computes the truncation error scaling for NDSM by applying 
# the code to a known (analytic) test case. NDSM is a second-order
# scheme, so ideally, estimates of the truncation error should 
# decrease with the square of the mesh spacing, h, for a uniform 
# mesh. In other words, the gamma output by this script should 
# be around 2. 
#
# Error printed to screen is in the form: Ea_max,Ea_avg,Eb_max,Eb_avg,time
# where a and b indicate error in vector potential and magnetic field
# respectively. Also reported is runtime in seconds
#
#

# Standard
import os
import sys
import numpy as np
import time

# Local
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../fortran/')))
from utests import power_law_fit
from ndsm import vector_potential

# ---------------------------------------------------------------------

def test_case1(X,Y,Z):
  """
    Potential field test case. Note that to serve as a test case 
    for the vector potential, the wave number, wn, must be chosen 
    as PI*N, where N is an integer. This is necessary so that Ax and Ay
    satisfy the boundary conditions assumed by NDSM. 
 
    Similarly, the mesh must be chosen correctly so that picking 
    wn = PI*N gives the correct boundary conditions. 

    Parameters:
    -----------
      X,Y,Z: (nz,ny,nx)
        Cartesian mesh. Should cover domain [0,1]x[0,1]x[0,1]
    Returns:
    --------
      A: (3,nz,ny,nx)
        Vector potential
      b: (3,nz,ny,nx)
        Vector magnetic field
  """

  # Get shape
  nz,ny,nz = X.shape

  # Set parameters
  wn = 1.*np.pi
  l  = np.sqrt(2*wn**2)

  # Output arrays
  b = np.zeros((3,)+X.shape)
  A = np.zeros((3,)+X.shape)

  # Magnetic field
  b[0,:,:,:] = +l*np.sin(wn*X)*np.cos(wn*Y)*np.exp(-l*Z)
  b[1,:,:,:] = +l*np.cos(wn*X)*np.sin(wn*Y)*np.exp(-l*Z)
  b[2,:,:,:] = +2*wn*np.cos(wn*X)*np.cos(wn*Y)*np.exp(-l*Z)

  # Vector potential
  A[0,:,:,:] = -np.cos(wn*X)*np.sin(wn*Y)*np.exp(-l*Z)
  A[1,:,:,:] = +np.sin(wn*X)*np.cos(wn*Y)*np.exp(-l*Z)

  return A,b

# ---------------------------------------------------------------------

# Base mesh
nshape_base = np.array([22,22,22])

# Scale up factors
scale_factors = [1,2,3,3.5,4,4.5,7.3,8,10]

# Error names
enames = ["Ea_max","Ea_avg","Eb_max","Eb_avg","Time"]

# array to hold errors and mesh spacing
dx    = np.zeros(len(scale_factors))
errors = np.zeros([len(enames),len(scale_factors)])

# Perform calculation at different resolutions
for i,scale in enumerate(scale_factors):

  # Scaleup number of points
  nshape = (scale*nshape_base).astype(int)

  # Mesh
  nz,ny,nx = nshape
  x        = np.linspace(0,1,nx)
  dx[i]    = x[1]-x[0]
  y        = np.arange(ny)*dx[i]
  z        = np.arange(nz)*dx[i]
  
  # Full mesh
  Z,Y,X = np.meshgrid(z,y,x,indexing='ij')

  # Build test case
  A1,b1 = test_case1(X,Y,Z)

  # Compute vector potential and magnetic field
  bc         = b1.copy()
  t1         = time.time()
  ierr,A2,b2 = vector_potential(x,y,z,bc)
  t2         = time.time()
  dt         = t2-t1

  # Error metrics
  Eb_avg1 = (np.linalg.norm(b1-b2,axis=0)).mean()
  Ea_avg1 = (np.linalg.norm(A1-A2,axis=0)).mean()
  Eb_max1 = (np.linalg.norm(b1-b2,axis=0)).max()
  Ea_max1 = (np.linalg.norm(A1-A2,axis=0)).max()

  # Save
  evec        = [Ea_max1,Ea_avg1,Eb_max1,Eb_avg1,dt]
  errors[:,i] = evec

  # Print results
  out = "{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}\t{:.5e}".format(dx[i],*evec)
  print(out)

# Fit power law and report indices
for i,name in enumerate(enames):
  gamma,A,ev = power_law_fit(dx,errors[i,:])
  print("Power-law index {:s}: {:g}".format(name,gamma))

