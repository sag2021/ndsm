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
# Help scripts for testing
#

import numpy as np

def power_law_fit(x,y,cov=False):
  """
    Fit a power law y = A*x^gamma
      
    Parameters:
    -----------
      x: array
        Array of x values
      y: array
        Array of y values

    Returns:
    --------
      gamma: float
        power-law index 
      A: float
        Coefficent of x^gamma
      ev: lambda
        Lambda that returns y = A*x^gamma for the given parameters

  """
  Lx    = np.log10(x)
  Ly    = np.log10(y)
  if(cov):
    p,C    = np.polyfit(Lx,Ly,1,cov=True,full=False)
  else:
    p    = np.polyfit(Lx,Ly,1)
  A     = 10.**p[1]
  ev    = lambda x : A*x**p[0]
 
  if(cov):
    return p[0],A,C,ev 
  else:
    return p[0],A,ev 

