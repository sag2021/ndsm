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
# Unit test for 2D solver. The corresponding Fortran code 
# needs to be run first. 
#
# Usage:
#   make 
#   ./unit_test_2D_solve 
#   python unit_test_2D_solve.py
#
#

import os
import sys
import numpy as np

# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

# Local
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from utests import power_law_fit

# Plot and text parameters
plt.rc('text' , usetex=True)
plt.rc('font' , family='serif')
plt.rc('xtick', labelsize= 16)  
plt.rc('ytick', labelsize= 14)  
mpl.rcParams['axes.linewidth'] = 1.5
LABEL_SIZE = 16

# Read data. 
res = np.loadtxt("res.txt")

# Extract
h    = res[:,0]
Emax = res[:,1]

# Fit power law
gamma,A,ev = power_law_fit(h,Emax)
hc         = np.logspace(-3,0,16)
E          = ev(hc) # Ideal Power-law error

# Report power-law index
print ("Power-law index: {:.12g}".format(gamma))

# LogLog plot
plt.loglog(h,Emax,'.',label="$E_{\\rm max}$",zorder=4,color='b')
plt.loglog(hc,E  ,'-',label="$E_{\\rm max} \propto h^\gamma : \gamma=%3.3f$"%gamma,zorder=4,color='r')

# Labels etc.
plt.title("NDSM 2D Poisson solver test")
plt.xlabel("Mesh spacing: $h$",fontsize=LABEL_SIZE)
plt.ylabel("Numerical Error ($E_{\\rm max}$)",fontsize=LABEL_SIZE)
plt.grid('on',which='both',color='.8',linestyle='-')
plt.legend(loc='upper left',fontsize=14)
plt.minorticks_on()
plt.xlim([1e-3,1])

  # Write plot
oname ="unit_test_2D_solve.pdf"
print ("Writing: "+oname)
plt.savefig(oname,bbox_inches="tight")

