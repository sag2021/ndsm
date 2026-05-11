#!/bin/mask
#
# Makefile unit tests
#

VPATH = ./:../../fortran

# FORTRAN files and compiler
FC      = gfortran
FFLAGS  = -O3 -fopenmp  -std=f2003
SRC     = ndsm_root.f90  ndsm_interp.f90 ndsm_multigrid_core.f90\
          ndsm_optimized.f90 ndsm_poisson.f90 \
          ndsm_vector_potential.f90

# Object files
OBJS = $(SRC:.f90=.o) 

EXEC = test

test  : $(OBJS) dev_test2.o
	$(FC) $(OBJS) dev_test2.o -o $(EXEC) $(FFLAGS)

%.o : %.f90
	$(FC) $< -o $@ -c $(FFLAGS)      

clean :
	-\rm  *.o  *.mod $(EXEC) 

