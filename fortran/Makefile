#!/bin/mask
#
# Makefile for NDSM
#

# Fortran compiler and flags
FC      = gfortran
FFLAGS  = -O3  -fPIC  -std=f2003 -fopenmp 

# Source files
SRC     = ndsm_root.f90  ndsm_interp.f90 ndsm_multigrid_core.f90\
          ndsm_optimized.f90 ndsm_poisson.f90 \
          ndsm_vector_potential.f90 ndsm_python_wrapper.f90 

# Object files
OBJS = $(SRC:.f90=.o) 

# Shared lib. name
EXEC = ndsmf.so

$(EXEC)	: $(OBJS) 
	$(FC) $(OBJS) -o $(EXEC) $(FFLAGS) -shared

%.o : %.f90
	$(FC) $< -o $@ -c $(FFLAGS)      

clean :
	-\rm $(EXEC) *.o  *.mod 
