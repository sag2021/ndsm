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
!>@name ndsm_multigrid_core
!
!>@brief Core module for N-dimensional multigrid
!
!>@details
!!
!! This module defines the MG_HANDLE type which is used 
!! to hold data/meta-data for multigrid calculation.
!!
!! Also defined are the basic subroutines from moving between
!! the coarse and fine problem and V cycling.
!!
!! This module defines restriction/interpolation operators, but it 
!! does not defined relaxation/residual operations. The latter
!! tend to be more problem dependent, as so are passed as 
!! externally defined procedures, in order to facilitate 
!! customization.
!!
!! The MG_HANDLE type includes optional arguments iopt,ropt, and
!! copt. The internal algorithm does not depend on these variables,
!! and may be used by a user to pass data.
!!
!! The module also defines some basic subroutines for performing
!! copying and initialization in parallel. 
!!
!
!>@author S.A. Gilchrist

MODULE NDSM_MULTIGRID_CORE

 USE NDSM_ROOT
 
 IMPLICIT NONE
   
  ! Length of user defined option
  INTEGER(IT),PARAMETER :: MG_OPT_LEN = 256 !< Length of MG_HANDLE option vectors
    
  ! =================
  ! DEFINE: MG_HANDLE
  ! =================
  !
  ! 
  ! 
  !
  !>@name mg_handle
  !>@brief Multrigrid handle
  !>@details
  !!
  !! This class encapsulates the basic data structures needed for
  !! a multigrid calculation in NDIM dimensions
  !!
  !!
  !! The MG_HANDLE type includes optional arguments iopt,ropt, and
  !! copt. The internal algorithm does not depend on these variables,
  !! and may be used by a user to pass data.
  !
  !
  !
  TYPE :: MG_HANDLE
    INTEGER(IT)                             :: ms       = -1 !< Number of relaxation sweeps for smoothing
    INTEGER(IT)                             :: ngrids   = -1 !< Number of grids 
    REAL(FP)                                :: ex_tol   = -1 !< Tolerance for coarsest grid solution 
    INTEGER(IT)                             :: ndim     = -1 !< Number of dimensions
    INTEGER(IT),DIMENSION(:,:) ,ALLOCATABLE :: nshape        !< Sizes (in mesh points) of each grid
    INTEGER(IT),DIMENSION(:)   ,ALLOCATABLE :: nsize         !< Total number of points in each grid
    TYPE(MG_PTR),DIMENSION(:)  ,ALLOCATABLE :: u             !< Dependent variable
    TYPE(MG_PTR),DIMENSION(:)  ,ALLOCATABLE :: rhs           !< Right-hand side
    TYPE(MG_PTR),DIMENSION(:,:),ALLOCATABLE :: meshes        !< Meshes for each grid   
    INTEGER(IT)     ,DIMENSION(MG_OPT_LEN)  :: iopt          !< User-defined integer options
    REAL(FP)        ,DIMENSION(MG_OPT_LEN)  :: ropt          !< User-defined real options
    CHARACTER(LEN=1),DIMENSION(MG_OPT_LEN)  :: copt          !< User-defined character options
    LOGICAL                                 :: du_max        !< Use max as convergence metric
  END TYPE
  
  !
  ! Abstract interfaces for relaxation and residual operators
  !
  ABSTRACT INTERFACE
  
    !>@name mg_relax
    !>@brief Interface for relaxation subroutine
    !>@details
    !! Performs relaxation on a given grid. The grid id (g_id)
    !! is passed to select a particular grid stored in the 
    !! MG_HANDLE object 
    !!
    !! Array size/shape info is passed via the MG_HANDLE object
    !!
    SUBROUTINE MG_RELAX(this,g_id,u,rhs)
      IMPORT :: FP,IT,MG_HANDLE
      TYPE(MG_HANDLE)      ,TARGET,INTENT(INOUT)  :: this
      INTEGER(IT)                 ,INTENT(IN)     :: g_id
      REAL(FP),DIMENSION(*),TARGET,INTENT(IN)     :: rhs
      REAL(FP),DIMENSION(*),TARGET,INTENT(INOUT)  :: u      
    END SUBROUTINE

   !>@name mg_residual
   !>@brief Interface for residual subroutine
   !
   SUBROUTINE MG_RESIDUAL(this,g_id,u,rhs,r)
      IMPORT :: FP,IT,MG_HANDLE
      TYPE(MG_HANDLE)      ,TARGET,INTENT(INOUT)  :: this
      INTEGER(IT)                 ,INTENT(IN)     :: g_id 
      REAL(FP),DIMENSION(*),TARGET,INTENT(IN)     :: u,rhs
      REAL(FP),DIMENSION(*),TARGET,INTENT(OUT)    :: r      
    END SUBROUTINE
            
  END INTERFACE 

  ! Internal procedures
  PRIVATE :: fine_to_coarse
  PRIVATE :: coarse_to_fine
  PRIVATE :: copy,zero
  
 CONTAINS

! ----------------------------------------------------------------------
!
! CONSTRUCTOR/DESTRUCTOR
!
! ----------------------------------------------------------------------

!>@name new_handle_base
!
!>@brief Basic constructor for MG_HANDLE
!
!>@details
!!
!! Initializes MG_HANDLE type. 
!!
!! Allocates %u and %rhs, but not their fields, 
!! i.e. %u(1)%val is unallocated on return 
!! 
!! Does not initialize iopt,ropt,copt. These are user defined.
!!
!
SUBROUTINE new_mg_handle(this,ndim,nshape,ngrids,mesh,du_max)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                 ,INTENT(IN) :: ndim
  INTEGER(IT)                 ,INTENT(IN) :: ngrids
  INTEGER(IT),DIMENSION(ndim) ,INTENT(IN) :: nshape
  TYPE(MG_PTR),DIMENSION(ndim),INTENT(IN) :: mesh
  LOGICAL                     ,INTENT(IN) :: du_max
  
  ! OUTPUT
  TYPE(MG_HANDLE),INTENT(OUT) :: this
  
  ! LOCAL
  INTEGER(IT) :: i,j,n,nq
  REAL(FP)    :: Lq,qil,dq
  
  ! LOCAL PARAM.
  REAL(FP),PARAMETER :: inv2 = REAL(1,FP)/REAL(2,FP) ! 1/2
    
  ! Set conv. flag
  this%du_max = du_max

  ! Set total number of grids
  this%ngrids = ngrids
  
  ! Set number of dimensions
  this%ndim = ndim

  ! ================
  ! GRID SIZES 
  ! ================
  !
  ! Determine the number of grid points for 
  ! each grid level
  
  ! Allocate memory for sizes
  ALLOCATE(this%nshape(ndim,ngrids))
  ALLOCATE(this%nsize(ngrids))
  
  ! Finest grid 
  this%nshape(:,1) = nshape
    
  ! Compute grid sizes. Each coarser grid has half as many points
  ! as the previous one, but never goes below one.
  DO i=2,ngrids
    this%nshape(:,i) = MAX(FLOOR(this%nshape(:,i-1)*inv2),1)
  ENDDO  
  
  ! Construct the total size of each grid 
  DO i=1,ngrids
    this%nsize(i) = PRODUCT(this%nshape(:,i))
  ENDDO
  
  ! ==========
  ! MESH
  ! ==========
  !
  ALLOCATE(this%meshes(ndim,ngrids))

  ! Finest grid
  DO i=1,ndim
    ALLOCATE(this%meshes(i,1)%val(nshape(i)))
    
    DO j = 1,nshape(i)
      this%meshes(i,1)%val(j) = mesh(i)%val(j)  
    ENDDO
    
  ENDDO
  
  !
  ! Coarser meshes
  !
  DO n=2,ngrids
  
    DO i=1,ndim  
    
      nq = this%nshape(i,n)
    
      ! Allocate
      ALLOCATE(this%meshes(i,n)%val(nq))
    
      ! Extent
      qil = MINVAL(mesh(i)%val)
      Lq  = MAXVAL(mesh(i)%val) - qil  
         
      ! Generate
      DO j=1,nq
        this%meshes(i,n)%val(j) = (j-1)*Lq/REAL(nq-1,FP) + qil 
      ENDDO
   
    ENDDO
    
  ENDDO
    
  ! Data arrays. Not individual elements are not allocated, 
  ! i.e. this%u(1) is unallocated
  ALLOCATE(this%u(ngrids))
  ALLOCATE(this%rhs(ngrids))
  
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name delete_mg_handle
!
!>@brief Destructor for MG_HANDLE type
!
SUBROUTINE delete_mg_handle(this)

  IMPLICIT NONE
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(INOUT) :: this
  
  ! LOCAL
  INTEGER(IT) :: i,j
  
  ! Free: u 
  IF(ALLOCATED(this%u)) THEN
    DO i=1,SIZE(this%u)
      IF(ALLOCATED(this%u(i)%val)) DEALLOCATE(this%u(i)%val)
    ENDDO
  DEALLOCATE(this%u)
  ENDIF
  
  ! Free: RHS 
  IF(ALLOCATED(this%rhs)) THEN
    DO i=1,SIZE(this%rhs)
      IF(ALLOCATED(this%rhs(i)%val)) DEALLOCATE(this%rhs(i)%val)
    ENDDO
  DEALLOCATE(this%rhs)
  ENDIF
  
  IF(ALLOCATED(this%meshes)) THEN
    DO j=1,SIZE(this%meshes,2)
      DO i=1,SIZE(this%meshes,1)
        IF(ALLOCATED(this%meshes(i,j)%val)) DEALLOCATE(this%meshes(i,j)%val)
      ENDDO
    ENDDO
    DEALLOCATE(this%meshes)
  ENDIF
  
  ! Free size vectors
  IF(ALLOCATED(this%nshape)) DEALLOCATE(this%nshape)
  IF(ALLOCATED(this%nsize )) DEALLOCATE(this%nsize)

  ! Set field to invalid values. This is to try an force a 
  ! crash in the event that a stale (deleted) instance of the
  ! object is called
  !
  this%ms      = -HUGE(this%ms)
  this%ngrids  = -HUGE(this%ngrids)
  this%ex_tol  = -HUGE(this%ex_tol) 
  this%ndim    = -HUGE(this%ndim)
  this%ropt    = -HUGE(this%ropt)
  this%iopt    = -HUGE(this%iopt)
  this%copt    = "!"
   
END SUBROUTINE

! ----------------------------------------------------------------------
!
! MULTIGRID-CYCLES
!
! ----------------------------------------------------------------------

!>@name v_cycle
!
!>@brief Perform a single V cycle starting from grid_id
!
SUBROUTINE v_cycle(this,grid_id,relax,residual)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT),INTENT(IN) :: grid_id
  PROCEDURE(MG_RELAX)    :: relax
  PROCEDURE(MG_RESIDUAL) :: residual
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(INOUT) :: this
  
  ! LOCAL
  INTEGER(IT) :: j 
  
  ! ===================
  ! PERFORM A V CYCLE 
  ! ===================
 
  !CALL one_grid(this,INT(1,IT),relax)
  !PRINT *,"ONE GRID"
  !RETURN

  ! Go to coarsest grid
  DO j=grid_id,this%ngrids-1
    CALL fine_to_coarse(this,j,relax,residual)
  ENDDO
  
  ! Solve exactly on coarsest grid
  CALL solve_exact(this,this%ngrids,relax)
  
  ! Go back to fine grids and apply corrections
  DO j=this%ngrids,grid_id+1,-1   
    CALL coarse_to_fine(this,j,relax)
  ENDDO
 
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name two_grid
!
!>@brief Two-grid correction scheme. Useful for testing.
!
SUBROUTINE two_grid(this,grid_id,relax,residual)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT),INTENT(IN) :: grid_id
  PROCEDURE(MG_RELAX)    :: relax
  PROCEDURE(MG_RESIDUAL) :: residual 
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(INOUT) :: this

  ! LOCAL
  INTEGER(IT),PARAMETER :: i1 = 1
  INTEGER(IT),PARAMETER :: i2 = 2
  
  ! Go to coarse grid
  CALL fine_to_coarse(this,i1,relax,residual)
  
  ! Solve exactly on coarse grid
  CALL solve_exact(this,i2,relax)
   
  ! Go back to fine grid
  CALL coarse_to_fine(this,i2,relax)   
 
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name one_grid
!
!>@brief Computes solution using only one grid. Useful for testing relaxation/exact solution operator
!
!>@details
!!
!! Does not perform multigrid. Just calls the exact solution
!! subroutine on the finest mesh
!!

SUBROUTINE one_grid(this,grid_id,relax)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT),INTENT(IN) :: grid_id
  PROCEDURE(MG_RELAX)    :: relax
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(INOUT) :: this

  ! LOCAL
  INTEGER(IT),PARAMETER :: i1 = 1
   
  ! Solve exactly on coarse grid
  CALL solve_exact(this,i1,relax)
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
! GRID REFINE/COARSEN SUBROUTINES
!
! ---------------------------------------------------------------------

!>@name fine_to_coarse
!
!>@brief Move problem from fine to coarse grid
!
!>@details
!!
!! Moves problem from fine grid to coarse grid.
!!
!! Calculates the residual on the fine grid and then 
!! restricts it to a coarser grid. In the process
!! memory is allocated for the source and solution
!! on the coarser grid.
!!
!! Fine grid BVP: 
!!
!!     L[u_f] = rhs_f 
!!
!! Coarse grid BVP:
!!
!!     L[u_c] = rhs_c = R[r_f]
!!
!! Where R[.] is the restriction operator, and r_f = L[u_f] - rhs_f,
!! i.e. the residual on the fine grid.
!!
!! ACTIONS PERFORMED
!! =================
!!
!! STEP 1: Compute residual on fine grid:  r_f = L[u_f] - rhs_f 
!!
!! STEP 2: Restrict r_f to coarse grid: rhs_c = R[r_f]
!!
!
!
SUBROUTINE fine_to_coarse(this,grid_id,relax,residual)

  IMPLICIT NONE

  ! INPUT 
  INTEGER(IT),INTENT(IN) :: grid_id
  PROCEDURE(MG_RELAX)    :: relax
  PROCEDURE(MG_RESIDUAL) :: residual 

  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),TARGET,INTENT(INOUT)  :: this
    
  ! LOCAL
  INTEGER(IT) :: i 

  ! COARSE GRID VARIABLES
  INTEGER(IT)  :: id_c
  INTEGER(IT)  :: nsize_c             
  
  ! FINE GRID VARIABLES
  INTEGER(IT)                       :: id_f,nsize_f
  REAL(FP),DIMENSION(:),ALLOCATABLE :: r_f    ! Residual 
  REAL(FP),DIMENSION(:),POINTER     :: u_f    ! Solution
  REAL(FP),DIMENSION(:),POINTER     :: rhs_f  ! RHS
        
  ! Set "grid ID" of the coarse and fine grids. Finest grid 
  ! has ID 1. 
  id_f = grid_id
  id_c = grid_id + 1

  ! Determine the size of the fine grid
  nsize_f = this%nsize(id_f)
  nsize_c = this%nsize(id_c)
 
  ! Get points to the solution and the RHS on the fine grid
  u_f   => this%u(id_f)%val
  rhs_f => this%rhs(id_f)%val

  ! =================
  ! PRESMOOTH ON FINE
  ! =================
  DO i=1,this%ms
    CALL relax(this,id_f,u_f,rhs_f)
  ENDDO
 
  ! Allocate memory for residual on fine grid
  ALLOCATE(r_f(nsize_f))
  CALL zero(r_f,nsize_f)

  ! ========================
  ! COMPUTE RESIDUAL ON FINE
  ! ========================
  !
  ! Compute r = L(u) - rhs on fine grid
  !

  ! Calculate residual on fine grid
  CALL residual(this,id_f,u_f,rhs_f,r_f)
  
  ! Allocate memory for the 'f' on coarse grid initialize to zero
  ALLOCATE (this%rhs(id_c)%val(nsize_c))
  CALL zero(this%rhs(id_c)%val,nsize_c)

  ! =======================
  ! RESTRICT TO COARSE GRID
  ! =======================
  !
  ! Perform restriction: Compute RHS function on coarse grid from 
  ! fine-grid residual
  CALL mg_restrict(this,id_f,r_f,id_c,this%rhs(id_c)%val)

  ! Variable r is no longer needed
  DEALLOCATE(r_f)
  
  ! Allocate memory for 'u' on coarse grid.
  ALLOCATE(this%u(id_c)%val(nsize_c))
  CALL zero(this%u(id_c)%val,nsize_c)

END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name coarse_to_fine
!
!>@brief Moves problem from coarse to fine grid
!
!>@details
!!
!! Moves from coarse grid to fine grid. First solution u is pre-smoothed
!! on the course grid. Coarse grid u is interpolated onto the finer
!! grid to give the fine grid correction c. Finally the correction
!! c is added to the solution u on the fine grid and post-smoothed. 
!!
!! Fine grid BVP: 
!!
!!      L[u_f] = rhs_f 
!!
!! Coarse grid BVP:
!!
!!      L[u_c] = rhs_c = R[r_f]
!!
!! Where R[.] is the restriction operator, and r_f = L[u_f] - rhs_f,
!! i.e. the residual on the fine grid.
!!
!! ACTIONS PERFORMED
!! =================
!!
!! STEP 1: Compute correction cor_f = I[cor_c]
!!
!! STEP 2: Correct fine-grid problem: u_f -> u_f + cor_f
!
SUBROUTINE coarse_to_fine(this,grid_id,relax)

  IMPLICIT NONE 

  ! INPUT
  INTEGER(IT),INTENT(IN)  :: grid_id
  PROCEDURE(MG_RELAX)     :: relax
 
   ! INPUT/OUTPUT
  TYPE(MG_HANDLE),TARGET,INTENT(INOUT) :: this

  ! LOCAL 
  INTEGER(IT) :: i 

  ! COARSE GRID VARIABLES
  INTEGER(IT) :: id_c
  INTEGER(IT) :: nsize_c
  REAL(FP),DIMENSION(:),POINTER :: u_c         ! Solution
  REAL(FP),DIMENSION(:),POINTER :: rhs_c       ! RHS
  
  ! FINE GRID VARIABLES
  INTEGER(IT) :: id_f
  INTEGER(IT) :: nsize_f
  REAL(FP),DIMENSION(:),ALLOCATABLE  :: cor_f       ! Correction
  REAL(FP),DIMENSION(:),POINTER      :: u_f         ! Solution
  REAL(FP),DIMENSION(:),POINTER      :: rhs_f       ! RHS
  
  ! =============================
  ! SETUP POINTER ALIASES
  ! =============================
    
  ! Set "grid ID" of the coarse and fine grids. Finest grid 
  ! has ID 1. 
  id_f = grid_id - 1
  id_c = grid_id 

  ! Determine the size of the fine grid
  nsize_f = this%nsize(id_f)
  nsize_c = this%nsize(id_c)
  
  ! Get points to grids
  u_c   => this%u(id_c)%val
  rhs_c => this%rhs(id_c)%val
  u_f   => this%u(id_f)%val
  rhs_f => this%rhs(id_f)%val

  ! =============================
  ! PRE-SMOOTH ON COARSE GRID
  ! =============================
  DO i=1,this%ms
    CALL relax(this,id_c,u_c,rhs_c)
  ENDDO

  ! Deallocate memory for coarse grid RHS
  DEALLOCATE(this%rhs(id_c)%val)

  ! Allocate memory for the correction on the fine grid
  ALLOCATE(cor_f(nsize_f))
  CALL zero(cor_f,nsize_f)

  ! =============================
  ! INTERPOLATE TO FINE GRID
  ! =============================
  !
  ! Interpolate the correction to the fine grid
  !
  CALL mg_interp(this,id_f,cor_f,id_c,u_c)
     
  ! Deallocate memory for uc.
  DEALLOCATE(this%u(id_c)%val)

  ! =============================
  ! APPLY CORRECTION ON FINE GRID
  ! =============================
  !
  ! Correct the solution on the fine grid
  !
  ! u_f = u_f + cor_f
  !
  CALL add_correction(u_f,cor_f,nsize_f)
  
  ! Deallocate memory for correction
  DEALLOCATE(cor_f)

  ! ========================
  ! POST-SMOOTH ON FINE GRID
  ! ========================
  DO i=1,this%ms
    CALL relax(this,id_f,u_f,rhs_f)
  ENDDO

END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name add_correction
!
!>@brief Add coarse grid correction

SUBROUTINE add_correction(u_f,cor_f,nsize_f)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: nsize_f
  REAL(FP),DIMENSION(nsize_f),INTENT(IN) :: cor_f

  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(nsize_f),INTENT(INOUT) :: u_f
    
  ! LOCAL
  INTEGER(IT) :: i

  !$OMP PARALLEL DO PRIVATE(i)
  DO i=1,nsize_f
    u_f(i) = u_f(i) + cor_f(i)
  ENDDO
  !$OMP END PARALLEL DO
  
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name solve_exact
!
!>@brief Solve the BVP exactly, i.e. iterate until convergence
!
!>@details
!!
!! Solves the BVP by applying the relaxation operator until
!! convergence. Since the relaxation is very slow, it must be 
!! performed on the coarsest grid, in accordance with the 
!! multigrid philosophy 
!!

SUBROUTINE solve_exact(this,grid_id,relax)

  IMPLICIT NONE

  ! INPUT 
  INTEGER(IT),INTENT(IN) :: grid_id
  PROCEDURE(MG_RELAX)    :: relax

  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),TARGET,INTENT(INOUT) :: this
  
  ! LOCAL
  REAL(FP) :: du
  REAL(FP),DIMENSION(:),POINTER     :: ue,fe
  REAL(FP),DIMENSION(:),ALLOCATABLE :: u_sav
  INTEGER(IT) :: i,nsize
  INTEGER(IT) :: nsteps
  
  ! Get size 
  nsize = this%nsize(grid_id)
  
  ! Get points to the current grid.
  ue => this%u(grid_id)%val
  fe => this%rhs(grid_id)%val
  
  ! Allocate memory for temp. storage of u
  ALLOCATE(u_sav(nsize))
  CALL zero(u_sav,nsize)

  ! Initially, set the change to a large value  
  du = HUGE(REAL(1,FP))
    
  !
  ! Solve exactly using the relaxation operator
  !
  nsteps = 0
  SOLVE_LOOP : DO     
 
    ! Leave loop if the change is below ex_tol
    IF(du .le. this%ex_tol) EXIT SOLVE_LOOP

    ! Perform relaxation
    CALL relax(this,grid_id,ue,fe)
 
    ! Compute maximum change in the solution
    du = du_max(u_sav,ue)
          
    ! Save previous iteration 
    CALL copy(u_sav,ue)  
    
    nsteps = nsteps + 1
    
  ENDDO SOLVE_LOOP
    
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name du_max
!
!>@brief Return max difference between u1 and u2 
!
FUNCTION du_max(u1,u2)

  IMPLICIT NONE
  
  ! INPUT
  REAL(FP),DIMENSION(:),INTENT(IN) :: u1,u2
  
  ! RETURN VALUE 
  REAL(FP) :: du_max
  
  ! LOCAL
  INTEGER(IT) :: i,nsize
  REAL(FP)    :: du_thread,du_i
  
  ! Get size
  nsize = SIZE(u1)
  
  ! Initialise to zero
  du_thread = 0
  
  ! Determine max in parallel 
  !$OMP  PARALLEL DO PRIVATE(i,du_i) &
  !$OMP& REDUCTION(max:du_thread) 
  DO i=1,nsize
  
    ! Determine difference at point i
    du_i = ABS(u1(i) - u2(i))
        
    ! Update max 
    du_thread = MAX(du_thread,du_i)

  ENDDO
  !$OMP END PARALLEL DO 
  
  ! Return value
  du_max = du_thread
  
END FUNCTION

! ----------------------------------------------------------------------
!
! INTERPOLATION 
!
! ----------------------------------------------------------------------

!>@name mg_interp
!
!>@brief Default interpolation method

SUBROUTINE mg_interp(this,id_f,u_f,id_c,u_c)

  USE NDSM_INTERP

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)          ,INTENT(IN) :: id_f,id_c
  REAL(FP),DIMENSION(*),INTENT(IN) :: u_c

  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(INOUT),TARGET :: this
 
  ! OUTPUT
  REAL(FP),DIMENSION(*),INTENT(OUT) :: u_f
 
  ! LOCAL
  INTEGER(IT)                      :: i,n
  INTEGER(IT)                      :: nsize_c,nsize_f
  INTEGER(IT),DIMENSION(this%ndim) :: nshape_c,nshape_f,ivec
  REAL(FP)   ,DIMENSION(this%ndim) :: q0
     
  ! Total size
  nsize_c = this%nsize(id_c)
  nsize_f = this%nsize(id_f)
 
  ! Get size of coarse/fine grid shape
  nshape_c = this%nshape(:,id_c)  
  nshape_f = this%nshape(:,id_f)
    
  ! Zero f 
  CALL zero(u_f,nsize_f)
  
  ! Interpolate coarse grid (uc) onto finder grid (uf)
  !
  !$OMP PARALLEL DO PRIVATE(i,n,q0,ivec)
  DO n=1,nsize_f
    
    ! Convert n -> (i1,i2,i3,...)
    !
    ivec = 0
    ivec = lin2nd(this%ndim,nshape_f,n)
    
    ! Unpack coordinate
    DO i=1,this%ndim
      q0(i) = this%meshes(i,id_f)%val(ivec(i))
    ENDDO
       
    ! Interpolate from coarse to fine grid
    !
    u_f(n) = 0
    u_f(n) =  ninterp(this%ndim,nsize_c,nshape_c,this%meshes(:,id_c),q0,u_c)
    
  ENDDO
  !$OMP END PARALLEL DO
    
END SUBROUTINE

! ----------------------------------------------------------------------
!
! RESTRICTION 
!
! ----------------------------------------------------------------------

!>@name mg_restrict
!
!>@brief Default restriction operator 
!
!>@details
!!
!! Restriction operator u_c = R[u_f] for computing the values of a function
!! (u_c) on the coarse grid given the values of the function (u_f)
!! on the fine grid. This operator is the adjoin of the trilinear 
!! interpolation operator.
!!
!!
!!  o------o------o
!!  | .  . | .  . |
!!  |      |      |
!!  | .  . | .  . |
!!  o------x------o
!!  | .  . | .  . |
!!  |      |      |  
!!  | .  . | .  . |
!!  o------o------o
!!
!!  Figure 1: Fine grid points (.) that contribute to the value 
!!            at the coarse grid point x in 2D. Note that this 
!!            is a cartoon, and the actual number of fine grid points
!!            may vary.
!!            
!!
!!  x: Coarse grid point being computed
!!  o: Neighboring coarse grid points
!!  .: Fine grid points enclosed by neighboring coarse points
!!
!!
!!  To determine the value of the function at a coarse grid point 
!!  (xc,yc,zc) it is necessary to perform the restriction over all 
!!  the fine grid points that are enclosed by the adjacent coarse 
!!  grid cells. This is because, when interpolating, the point x 
!!  contributes to all the fine grid points in surrounding cells.
!!
!!
!!    x---------------------o
!!    |                     |
!!    |                     |
!!    |                     |
!!    |          *<-------->|
!!    |          ^          |
!!    |          |          |
!!    |          v          |  
!!    o---------------------o
!!
!!    Figure 2: The distances that appear in the weight that 
!!              determines the contribution of the point at * to
!!              the coarse grid function at x.
!! 
!!                                  
!!   The restriction operator has the  
!!   form 
!!                    _
!!   u_c(xc,yc,zc) =  \
!!                    /_ w(x',y',z') * u_f(x',y',z')
!!
!!   where the sum is over points that surround (xc,yc,zc). The 
!!   weights and the exact points used in the summation are 
!!   determined by requiring that R[] is the adjoin of the trilinear
!!   interpolation operator P[].
!!
!!   The contribution from a given fine grid point (xf,yf,zf) 
!!   (shown by * in Figure 2) depends on the distance from the 
!!   point to the cell boundary 
!!
!!   w(xf,yf,zf) = dx2*dy2*dz2
!!
!!   dx2 = dx_c - |xf-xc|
!!   dy2 = dy_c - |yf-yc|
!!   dz2 = dz_c - |zf-zc|
!!
!!   Here dx_c is the coarse cell spacing in x 
!!   (which is assumed uniform).
!
!

SUBROUTINE mg_restrict(this,id_f,u_f,id_c,u_c)

  USE NDSM_INTERP

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)          ,INTENT(IN) :: id_f,id_c
  REAL(FP),DIMENSION(*),INTENT(IN) :: u_f
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),INTENT(IN) :: this
 
  ! OUTPUT
  REAL(FP),DIMENSION(*),INTENT(OUT) :: u_c

  ! LOCAL
  INTEGER(IT)                      :: i,n
  INTEGER(IT)                      :: nsize_c,nsize_f
  INTEGER(IT),DIMENSION(this%ndim) :: nshape_c,nshape_f,ivec
  REAL(FP)   ,DIMENSION(this%ndim) :: q0

  ! Total size
  nsize_c = this%nsize(id_c)
  nsize_f = this%nsize(id_f)
 
  ! Get size of coarse/fine grid shape
  nshape_c = this%nshape(:,id_c)  
  nshape_f = this%nshape(:,id_f)

  ! Zero u_c 
  CALL zero(u_c,nsize_c)
  
  !$OMP PARALLEL DO PRIVATE(n,i,q0,ivec) 
  DO n=1,nsize_c

    ! Convert n -> (i1,i2,i3,...)
    !
    ivec = 0
    ivec = lin2nd(this%ndim,nshape_c,n)
    
    ! Get coordinate at point (ivec)
    DO i=1,this%ndim
      q0(i) = this%meshes(i,id_c)%val(ivec(i))
    ENDDO
        
    ! Perform restriction
    u_c(n) = nrestrict(this%ndim,nsize_c,nsize_f   , &
                                nshape_c,nshape_f  , &
                                this%meshes(:,id_c), &
                                this%meshes(:,id_f), &
                                q0,u_f)
  ENDDO
  !$OMP END PARALLEL DO
    
END SUBROUTINE

! ---------------------------------------------------------------------------
!
! UPDATE SUBROUTINES
!
! ---------------------------------------------------------------------------

!>@name update_u
!
!>@brief Updates u and computed convergence metrics

SUBROUTINE update_u(nsize,u_old,u_new,du_max,du_mean)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)              ,INTENT(IN) :: nsize
  REAL(FP),DIMENSION(nsize),INTENT(IN) :: u_old
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: u_new
  REAL(FP),OPTIONAL        ,INTENT(OUT) :: du_max,du_mean
  
  ! LOCAL
  INTEGER(IT) :: i
  REAL(FP)    :: du_max2,du_mean2,du
  
  ! Zero metrics
  du_max2  = 0
  du_mean2 = 0
    
  !$OMP  PARALLEL DO PRIVATE(i,du) & 
  !$OMP& REDUCTION(max:du_max2)    &  
  !$OMP& REDUCTION(+:du_mean2) 
  DO i=1,nsize
      
    ! Compute change
    du = ABS(u_new(i) - u_old(i))
       
    ! Compute metrics
    du_max2  = MAX(du,du_max2)
    du_mean2 = du_mean2 + du
              
    ! Update
    u_new(i) = u_old(i)
      
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Divide mean by the total 
  du_mean2 = du_mean2/REAL(SIZE(u_old),FP)
  
  ! Output if necessary
  IF(PRESENT(du_mean)) du_mean = du_mean2
  IF(PRESENT(du_max))  du_max  = du_max2
    
END SUBROUTINE

! ----------------------------------------------------------------------
!
! PARALLEL INITIALIZATION AND COPYING SUBROUTINES
!
! ----------------------------------------------------------------------

!>@name mg_zero
!
!>@brief Set array to zero using OpenMP
!
SUBROUTINE zero(array,nsize)

  IMPLICIT NONE 

  ! INPUT
  INTEGER(IT),INTENT(IN)  :: nsize
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: array

  ! LOCAL
  INTEGER(IT) :: i  

  !$OMP PARALLEL DO PRIVATE(i)
  DO i=1,nsize
    array(i) = 0
  ENDDO
  !$OMP END PARALLEL DO

END SUBROUTINE 

! ----------------------------------------------------------------------
!
!>@name copy
!
!>@brief Copy RHS to LHS
!
SUBROUTINE copy(lhs,rhs)

  IMPLICIT NONE
  
  ! INPUT
  REAL(FP),DIMENSION(:),INTENT(IN) :: rhs 
  
  ! OUTPUT
  REAL(FP),DIMENSION(:),INTENT(OUT) :: lhs
  
  ! LOCAL
  INTEGER(IT) :: i,nsize
  
  ! Determine input size
  nsize = SIZE(rhs)
  
  ! Copy in parallel
  !$OMP PARALLEL DO PRIVATE(i)
  DO i=1,nsize  
    lhs(i) = rhs(i)  
  ENDDO
  !$OMP END PARALLEL DO
  
END SUBROUTINE

! ---------------------------------------------------------------------  
!
!>@name mean
!
!>@brief Compute the mean in parallel 
!
!>@details
!!
!!  Computes the mean using the simple summation formula w = w + a(n).
!!  The OpenMP loop is ordered to force consistency with 
!!  serial execution
!!

FUNCTION mean(nsize,u)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)              ,INTENT(IN) :: nsize
  REAL(FP),DIMENSION(nsize),INTENT(IN) :: u

  ! RETURN
  REAL(FP) :: mean
  
  ! LOCAL
  INTEGER(IT) :: n
  
  mean = 0
  !$OMP PARALLEL DO PRIVATE(n) REDUCTION(+:mean) ORDERED
  DO n=1,nsize
    mean = mean + u(n)
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Divide by size
  mean = mean/REAL(nsize,FP)
  
END FUNCTION

END MODULE
