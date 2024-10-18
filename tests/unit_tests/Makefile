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

all	: unit_test_galerkin.f90 unit_test_interp unit_test_2D_solve

unit_test_2D_solve  : $(OBJS) unit_test_2D_solve.o 
	            $(FC) $(OBJS) unit_test_2D_solve.o -o unit_test_2D_solve $(FFLAGS)

unit_test_galerkin  : $(OBJS) unit_test_galerkin.o 
	            $(FC) $(OBJS) unit_test_galerkin.o -o unit_test_galerkin $(FFLAGS)

unit_test_interp: $(OBJS) unit_test_interp.o
		$(FC) $(OBJS) unit_test_interp.o -o unit_test_interp $(FFLAGS)

%.o : %.f90
	$(FC) $< -o $@ -c $(FFLAGS)      

clean :
	-\rm  *.o  *.mod unit_test_galerkin unit_test_interp	unit_test_2D_solve

