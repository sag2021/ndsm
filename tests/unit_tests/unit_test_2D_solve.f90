! BSD-3-Clause 
!
! Copyright 2024 S.A Gilchrist
!
! Redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, 
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, 
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A !PARTICULAR PURPOSE ARE DISCLAIMED. 
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
! OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT !LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
! OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

!
!>@name unit_test_2D_solve
!>@brief Unit test for 2D solver
!>@details
!!
!! Tests 2D Poisson Neumann solver. Computes the solution to Poisson
!! equation with homogeneous Neumann boundary conditions. The source 
!! function is simple enough to write down an simple analytic solution.
!! The max. and mean absolute difference between the exact and numerical solution
!! is used to estimate truncation error. The calculation is performed
!! at different mesh sizes so that the error scaling can be estimated.
!!
!! The output is written to a file called res.txt in the format 
!! dx, Emax, Eavg
!!
!! Where dx is the mesh spacing, which is uniform. If the filedump flag
!! is set to true, then the solution will also be dumped to a crude
!! Fortran unformatted file.
!!
!
!
PROGRAM UNIT_TEST_2D_SOLVE

  USE NDSM_ROOT     
  USE NDSM_POISSON

  IMPLICIT NONE

  ! Mesh size
  INTEGER(IT),DIMENSION(2),PARAMETER :: nshape_base  = [27,36]        ! Mesh size [nx,ny]  
  INTEGER(IT)             ,PARAMETER :: BASE_GRID    = 2              ! Smallest grid size
  INTEGER(IT)             ,PARAMETER :: max_vcycles  = 256            ! Max. number of V-cycles

  ! SEED: Used for a1 and b1
  INTEGER,PARAMETER                :: seed = 2112
  INTEGER                          :: seed_len
  INTEGER,DIMENSION(:),ALLOCATABLE :: seed_vec

  ! Mesh extent
  REAL(FP),PARAMETER :: Lx = 1.D0

  ! Scale factors
  REAL(FP),DIMENSION(9),PARAMETER :: scalefac = [1.D0,1.5D0,2.D0,4.D0,5.5D0,10.D0,15.D0,20.D0,25.D0]
  
  ! Result of test
  REAL(FP)    :: res(3),a1,b1
  INTEGER(IT) :: ierr,nshape(2),i
  INTEGER(IT) :: u = 34 ! Output file
 
  ! Set random seed
  CALL RANDOM_SEED(SIZE=seed_len)
  ALLOCATE(seed_vec(seed_len))
  seed_vec = SEED
  CALL RANDOM_SEED(PUT=seed_vec)

  ! Random values for a1 and b1
  CALL RANDOM_NUMBER(a1)
  CALL RANDOM_NUMBER(b1)

  !
  ! Solve at different resolutions
  !
  PRINT *,"Output file: res.txt"
  PRINT *,"Solving..."
  OPEN(u,FILE='res.txt',ACTION='write')
  WRITE(u,*) "# Result dx,Emax,Eavg"

  ! Compute solution at different mesh spacings
  DO i=1,SIZE(scalefac)

    ! Scale mesh size and round up to INT
    nshape = CEILING(nshape_base*scalefac(i))

    ! Solve
    CALL solve_test_case(nshape,max_vcycles,1D-12,Lx,res,ierr,.False.)

    ! Check for errors
    IF(ierr .ne. IERR_POISSON_SUCCESS) PRINT *,"ERROR: FAILED TO CONVERGE"

    ! Write to file
    WRITE(u,*) res

  ENDDO
  CLOSE(u)

CONTAINS

!>@name solve_test_case
!>@brief Solves test case at given resolution 
!>@details
!!  Returns vector [dx,E_max,E_avg] 
!!  ierr returns the error flag for the Poisson solver
!!  
!!  If filedump is true, then the potential and the analytic 
!!  potential are dumped to a Fortran unformatted file called 
!!  dump.dat
!! 
!!
SUBROUTINE solve_test_case(nshape,max_vcycles,ex_tol,Lx,res,ierr,filedump) 

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT),DIMENSION(2),INTENT(IN) :: nshape
  INTEGER(IT)             ,INTENT(IN) :: max_vcycles
  REAL(FP)                ,INTENT(IN) :: Lx,ex_tol
  LOGICAL                 ,INTENT(IN) :: filedump

  ! OUTPUT
  REAL(FP),DIMENSION(3),INTENT(OUT) :: res
  INTEGER(IT)          ,INTENT(OUT) :: ierr

  ! Data arrays
  REAL(FP),DIMENSION(:,:),ALLOCATABLE :: u,ue,rhs
  INTEGER(IT)                         :: nsize
  REAL(FP)                            :: nmin
  INTEGER(IT)                         :: ndim = SIZE(nshape)

  ! Mesh 2D
  REAL(FP)                         :: dq
  TYPE(MG_PTR),DIMENSION(2),TARGET :: mesh
  REAL(FP)                         :: Ly

  INTEGER :: i,j

  REAL(FP),PARAMETER :: inv2 = REAL(1,FP)/REAL(2,FP)
  REAL(FP),PARAMETER :: inv3 = REAL(1,FP)/REAL(3,FP)

  ! BVP solver
  TYPE(MG_HANDLE) :: bvp
  REAL(FP)        :: du_last,mean_ue
  REAL(FP)        :: E_avg,E_max
  INTEGER(IT)     :: ngrids

  ! Alias
  REAL(FP),DIMENSION(:),POINTER :: x,y

  ! Mesh spacing: Uniform in all dimensions
  !
  dq = REAL(1,FP)/(nshape(1)-REAL(1,FP))
    
  ! Allocate memory
  DO i=1,SIZE(mesh)
    ALLOCATE(mesh(i)%val(nshape(i)))
  ENDDO
  
  ! Construct base mesh
  DO i=1,SIZE(mesh)
    mesh(i)%val = [(j,j=0,nshape(i)-1)]*dq
  ENDDO
  x => mesh(1)%val
  y => mesh(2)%val

  ! Get Ly
  Ly = MAXVAL(y)-MINVAL(y)

  ! Allocate data arrays
  ALLOCATE(u(nshape(1),nshape(2)))
  ALLOCATE(ue(nshape(1),nshape(2)))
  ALLOCATE(rhs(nshape(1),nshape(2)))
 
  ! Get total number of points
  nsize = SIZE(u)

  ! Build RHS
  DO j=1,SIZE(rhs,2)
    DO i=1,SIZE(rhs,1)
      rhs(i,j) = a1*(2*x(i)-Lx) + b1*(2*y(j)-Ly)
    ENDDO
  ENDDO

  ! Number of grids
  nmin   = MINVAL(nshape)
  ngrids = FLOOR(LOG(nmin/BASE_GRID)/LOG(REAL(2,FP)))

  ! New handle
  CALL new_mg_handle(bvp,ndim,nshape,ngrids,mesh)

  ! Set values
  bvp%ms        = 5
  bvp%ex_tol    = ex_tol
  bvp%copt(1:4) = "N"         

  ! Solve BVP 
  u = 0 
  CALL solve_poisson_bvp(bvp,nsize,ex_tol,max_vcycles,u,rhs,du_last,ierr) 

  ! Build RHS
  DO j=1,SIZE(rhs,2)
    DO i=1,SIZE(rhs,1)
      ue(i,j) = a1*x(i)**2*(inv3*x(i) - inv2*Lx) &
              + b1*y(j)**2*(inv3*y(j) - inv2*Ly)
    ENDDO
  ENDDO
 
  ! Remove mean
  mean_ue = SUM(ue)/REAL(SIZE(ue),FP)
  ue      = ue - mean_ue

  ! Compute mean and max error 
  res(1) = dq
  res(2) = MAXVAL(ABS(u-ue))
  res(3) = SUM(ABS(u-ue))/REAL(nsize,FP)

  ! Quick output dump
  IF(filedump) THEN
    PRINT *,"Dumping to file: dump.dat"
    OPEN(11,FILE='dump.dat',ACTION='write',FORM='unformatted')
      WRITE(11) nshape
      WRITE(11) u
      WRITE(11) ue
    CLOSE(11)
  ENDIF

END SUBROUTINE

END PROGRAM
