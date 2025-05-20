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


! ---------------------------------------------------------------------
!
!>@name NDSM_POISSON
!
!>@brief Defines the relaxation/residual operators for the Poisson equation

MODULE NDSM_POISSON

  USE NDSM_ROOT
  USE NDSM_MULTIGRID_CORE

  IMPLICIT NONE
  
  ! Boundary mask flags
  INTEGER,PARAMETER :: BND_NONE  = 0   !< Flag: Point is not at boundary
  INTEGER,PARAMETER :: BND_LOWER = 1   !< Flag: Point is at a lower boundary
  INTEGER,PARAMETER :: BND_UPPER = 2   !< Flag: Point is at an upper boundary

  ! ERROR FLAGS
  INTEGER,PARAMETER,PUBLIC :: IERR_POISSON_SUCCESS = 0
  INTEGER,PARAMETER,PUBLIC :: IERR_POISSON_COVFAIL = 1
      
 CONTAINS  

! ---------------------------------------------------------------------
!
!>@name solve_poission_bvp_ierr
!
!>@brief Solve the Poisson BVP defined by bvp
!>@details
!! 
!!
!! ierr = 0 OK
!! ierr = 1 Not OK
!!
 
SUBROUTINE solve_poisson_bvp(bvp,nsize,vc_tol,nmax,u,rhs,du_last,ierr) 

  IMPLICIT NONE
  
  ! SUBNAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB="solve_poisson_bvp"

  ! INPUT
  INTEGER(IT),INTENT(IN) :: nsize
  REAL(FP)   ,INTENT(IN) :: vc_tol
  INTEGER(IT),INTENT(IN) :: nmax
  
  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),TARGET   ,INTENT(INOUT) :: bvp
  REAL(FP),DIMENSION(nsize),INTENT(INOUT) :: u,rhs

  ! OUTPUT
  INTEGER(IT),INTENT(OUT) :: ierr
  REAL(FP)   ,INTENT(OUT) :: du_last
  
  ! LOCAL
  INTEGER(IT)                   :: i
  LOGICAL                       :: converged
  REAL(FP)                      :: du,du_max,du_mean
  REAL(FP),DIMENSION(:),POINTER :: u_p,rhs_p
  INTEGER(IT),PARAMETER         :: i1 = 1
  CHARACTER(LEN=32)             :: du_string

  ! Allocate memory for finest mesh
  IF(.NOT.ALLOCATED(bvp%u(1)%val  )) ALLOCATE(bvp%u(1)%val(nsize))
  IF(.NOT.ALLOCATED(bvp%rhs(1)%val)) ALLOCATE(bvp%rhs(1)%val(nsize))
  
  ! Construct pointers to finest mesh
  u_p   => bvp%u(1)%val
  rhs_p => bvp%rhs(1)%val

  ! Set finest mesh to the input values
  u_p   = u
  rhs_p = rhs
  
  ! Inititialize difference to largest possible value
  du = HUGE(du)
  
  ! Initialize convergence flag
  converged = .FALSE.

  ! Default
  ierr = IERR_POISSON_SUCCESS
  
  !
  ! Perform multigrid V cycles until convergence
  !
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Performing V cycles...")
  VCYCLE_LOOP : DO i=1,nmax

    ! Do V cycle
    CALL v_cycle(bvp,i1,ndsm_relax_wrapper,ndsm_residual_wrapper)

    ! Update and compute difference
    CALL update_u(nsize,u_p,u,du_max=du_max,du_mean=du_mean)
   
    ! Pick convergence metric
    IF(bvp%du_max) THEN
      du = du_max
    ELSE
      du = du_mean
    ENDIF

    ! Print change in solution
    WRITE(du_string,"(ES12.4)") du
    IF(DEBUG) CALL debug_msg(THIS_SUB,"Solution delta: "//TRIM(du_string))

    ! Exit loop if converged
    IF(du < vc_tol) THEN
      converged = .TRUE.
      EXIT VCYCLE_LOOP
    ENDIF 
  
  ENDDO VCYCLE_LOOP
  
  ! Return change in solution at last iteration
  du_last = du
    
  ! Warn the user if the V-cycle did not converge
  IF(.NOT.converged) THEN
    ierr = IERR_POISSON_COVFAIL
    PRINT *,"Warning: IOPT_NCYCLES exceeded. V-cycle iteration may not have converged"
  ENDIF
   
  ! Copy out
  u = u_p
 
END SUBROUTINE
 
! ---------------------------------------------------------------------
!
!>@name ndsm_residual_wrapper
!
!>@brief Class wrapper for residual

SUBROUTINE ndsm_residual_wrapper(this,g_id,u,rhs,r)

  USE NDSM_OPTIMIZED,ONLY: poisson_residual_3D

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                 ,INTENT(IN) :: g_id 
  REAL(FP),DIMENSION(*),TARGET,INTENT(IN) :: u,rhs

  ! INPUT/OUTPUT
  TYPE(MG_HANDLE),TARGET,INTENT(INOUT)  :: this

  ! OUTPUT
  REAL(FP),DIMENSION(*),TARGET,INTENT(OUT) :: r      
  
  ! LOCAL
  CHARACTER(LEN=1),DIMENSION(this%ndim,2) :: bcs
  INTEGER(IT)                             :: nc
  INTEGER(IT),PARAMETER                   :: i2 = 2
  
  nc  = 2*this%ndim
  bcs = RESHAPE(this%copt(1:nc),[this%ndim,i2])

  
  ! Call optimized routine for 3D
  IF(this%ndim == 3) THEN
    CALL poisson_residual_3D(bcs                   , &
                             this%nshape(1,g_id)    , &
                             this%nshape(2,g_id)    , &
                             this%nshape(3,g_id)    , &
                             this%meshes(1,g_id)%val, &
                             this%meshes(2,g_id)%val, & 
                             this%meshes(3,g_id)%val, &
                             this%rhs(g_id)%val     , & 
                             this%u(g_id)%val       , &
                             r                      )
    RETURN
  ENDIF
 
  
  !
  ! Call generic residual operator on grid with ID g_id
  !
  CALL residual(this%ndim          ,      &
                this%nsize(g_id)   ,      &
                this%nshape(:,g_id),      &
                this%meshes(:,g_id),      &
                bcs                ,      &
                this%u(g_id)%val   ,      &
                this%rhs(g_id)%val ,      & 
                r)
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name ndsm_relax_wrapper
!
!>@brief Wrapper for object class
 
SUBROUTINE ndsm_relax_wrapper(this,g_id,u,rhs) 
 
  USE NDSM_OPTIMIZED 
 
  IMPLICIT NONE 

  ! INPUT
  INTEGER(IT)                 ,INTENT(IN) :: g_id
  REAL(FP),DIMENSION(*),TARGET,INTENT(IN) :: rhs

  ! INPUT/OUTPUT
  TYPE(MG_HANDLE)        ,TARGET,INTENT(INOUT)  :: this
  REAL(FP),DIMENSION(*)  ,TARGET,INTENT(INOUT)  :: u      

  ! LOCAL
  CHARACTER(LEN=1),DIMENSION(this%ndim,2) :: bcs
  INTEGER(IT)                             :: nc
  INTEGER(IT),PARAMETER                   :: i2 = 2


  nc  = 2*this%ndim
  bcs = RESHAPE(this%copt(1:nc),[this%ndim,i2])

  ! Call optimized routine for 3D
  IF(this%ndim == 3) THEN
    CALL red_black_gauss_3D(bcs                    , &
                            this%nshape(1,g_id)    , &
                            this%nshape(2,g_id)    , &
                            this%nshape(3,g_id)    , &
                            this%meshes(1,g_id)%val, &
                            this%meshes(2,g_id)%val, & 
                            this%meshes(3,g_id)%val, &
                            this%rhs(g_id)%val     , & 
                            this%u(g_id)%val       )
    RETURN
  ENDIF
  
  !
  ! Call relaxation operator on grid with grid ID g_id
  !
  CALL relax(this%ndim          , &   
             this%nsize(g_id)   , &
             this%nshape(:,g_id), &
             this%meshes(:,g_id), &
             bcs                , &
             this%u(g_id)%val   , &
             this%rhs(g_id)%val)
  
END SUBROUTINE  
 
! ---------------------------------------------------------------------
!
!>@name residual
!
!>@brief Compute residual 
 
SUBROUTINE residual(ndim,nsize,nshape,qvecs,bcs,u,rhs,r)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                       ,INTENT(IN) :: ndim,nsize
  INTEGER(IT)     ,DIMENSION(ndim)  ,INTENT(IN) :: nshape
  TYPE(MG_PTR),DIMENSION(ndim)      ,INTENT(IN) :: qvecs 
  REAL(FP),DIMENSION(nsize)         ,INTENT(IN) :: u,rhs
  CHARACTER(LEN=1),DIMENSION(ndim,2),INTENT(IN) :: bcs
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: r
  
  ! LOCAL
  INTEGER(IT) ,DIMENSION(ndim) :: nstrides,ivec,bmask
  REAL(FP),DIMENSION(ndim)     :: dq,wc
  INTEGER(IT)                  :: i,n,dn(2)
  REAL(FP)                     :: laplace_u
  
  ! Compute strides
  nstrides = shape_to_strides(ndim,nshape)
      
  ! Compute mesh spacing and weights
  DO n=1,ndim
    dq(n) = qvecs(n)%val(2) - qvecs(n)%val(1)
    wc(n) = REAL(1,FP)/dq(n)**2
  ENDDO
    
  !
  ! Subtract Laplacian: r = rhs - Lap[u]
  !
  !$OMP PARALLEL DO PRIVATE(n,i,ivec,dn,bmask,laplace_u)
  ARRAY_LOOP : DO n=1,nsize
    
      ! Compute N-D vector
      ivec = 0
      ivec = lin2nd(ndim,nshape,n)
  
      ! Compute boundary mask
      bmask = 0
      CALL boundary_mask(ndim,ivec,nshape,bmask) 
  
      ! Check if point is at Dirichlet boundary. These points
      ! have zero residual
      IF(at_dirichlet_boundary(ndim,bmask,bcs)) THEN
        r(n) = 0
        CYCLE ARRAY_LOOP
      ENDIF

      ! ============
      ! STENCIL LOOP
      ! ============
      !
      laplace_u = 0
      
      STENCIL_LOOP : DO i=1,NDIM
    
        ! Strides
        dn = 0
        CALL stencil_stride(ndim,bmask,nstrides,i,dn)
        
        ! Apply stencil
        laplace_u = laplace_u + (u(n+dn(1)) - 2*u(n) + u(n+dn(2)))*wc(i) 
  
    ENDDO STENCIL_LOOP
    
      ! Compute residual
      r(n) =  rhs(n) - laplace_u
    
  ENDDO ARRAY_LOOP
  !$OMP END PARALLEL DO
 
END SUBROUTINE

! ---------------------------------------------------------------------  
!
!>@name at_dirichlet_boundary
!
!>@brief Returns true if bmask(i) is a Dirichlet boundary

LOGICAL FUNCTION at_dirichlet_boundary(ndim,bmask,bcs)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                       ,INTENT(IN) :: ndim
  INTEGER(IT),DIMENSION(ndim)       ,INTENT(IN) :: bmask
  CHARACTER(LEN=1),DIMENSION(ndim,2),INTENT(IN) :: bcs

  ! LOCAL
  LOGICAL :: lb,ub
  INTEGER(IT) :: i
  
  IF(ALL(bmask == BND_NONE)) THEN
    at_dirichlet_boundary = .FALSE.
    RETURN
  ENDIF
  
  DO i=1,NDIM
    lb = (bmask(i)==BND_LOWER) .AND. ( bcs(i,1) == "D" )
    ub = (bmask(i)==BND_UPPER) .AND. ( bcs(i,2) == "D" )
    IF(lb .OR. ub) THEN
      at_dirichlet_boundary = .TRUE.
      RETURN
    ENDIF
  ENDDO
  
  at_dirichlet_boundary = .FALSE.
  
END FUNCTION

! ---------------------------------------------------------------------  
!
!>@name boundary_mask
!
!>@brief Return flags indicating when ivec lies on a boundary
!
!>@details
!!
!! Returns a vector of integers (bmask) that flags when the 
!! point described by ivec lies on a boundar. Elements of he flag vector
!! (bmask) take three possible values:
!!
!! bmask(k) = BND_UPPER: Point is at the upper boundar of dimension k
!! bmask(k) = BND_LOWER: Point is at the lower boundary of dimension k
!! bmask(k) = BND_NONE:  Point is in the interior of dimension k
!!
!
PURE SUBROUTINE boundary_mask(ndim,ivec,nshape,bmask) 

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim 
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: ivec,nshape
  
  ! OUTPUT
  INTEGER(IT),DIMENSION(ndim),INTENT(OUT) :: bmask
  
  ! LOCAL
  INTEGER(IT) :: i
  
  DO i=1,ndim
    IF(ivec(i) == 1) THEN
      bmask(i) = BND_LOWER
    ELSEIF(ivec(i) == nshape(i)) THEN
      bmask(i) = BND_UPPER        
    ELSE
      bmask(i) = BND_NONE
    ENDIF 
  ENDDO
      
END SUBROUTINE

! ---------------------------------------------------------------------  
!
!>@name relax
!
!>@brief Perform a relaxation sweep in N dimensions
!
!>@details
!!
!! Performs a single Red-Black relaxation sweep relaxtion
!! for Poisson's equation
!!
!! Homogenous Neumman or Dirichlet boundary conditions are 
!! enforced at the boundaries. The input, bcs, determines
!! the specific 
!!
!
SUBROUTINE relax(ndim,nsize,nshape,qvecs,bcs,u,rhs)

  USE NDSM_INTERP
  
  IMPLICIT NONE

  ! INPUT 
  INTEGER(IT)                         ,INTENT(IN) :: ndim,nsize
  INTEGER(IT)         ,DIMENSION(ndim),INTENT(IN) :: nshape
  TYPE(MG_PTR)    ,DIMENSION(ndim)    ,INTENT(IN) :: qvecs
  REAL(FP)        ,DIMENSION(nsize)   ,INTENT(IN) :: rhs
  CHARACTER(LEN=1),DIMENSION(ndim,2)  ,INTENT(IN) :: bcs
  
  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(INOUT) :: u

  ! LOCAL
  INTEGER(IT)                 :: n,i
  INTEGER(IT),DIMENSION(ndim) :: nstrides,ivec,parity
  REAL(FP),DIMENSION(0:ndim)  :: wc
  REAL(FP)                    :: dq,u_mean
  LOGICAL                     :: is_red

  ! LOCAL LITS.
  INTEGER(IT),PARAMETER :: i2 = 2

  ! Compute strides
  !
  nstrides = shape_to_strides(ndim,nshape)

  ! Compute weights
  !
  wc(0) = 0
  DO i=1,NDIM
    dq    = qvecs(i)%val(2) - qvecs(i)%val(1)    
    wc(i) = REAL(1,FP)/dq**2  
    wc(0) = wc(0) + REAL(2,FP)*wc(i)
  ENDDO
  wc(0) = REAL(1,FP)/wc(0)
  
  !
  ! Relax: Red
  !
  !$OMP PARALLEL DO PRIVATE(n,ivec,parity,is_red)
  DO n=1,nsize
  
    ivec   = 0
    ivec   = lin2nd(ndim,nshape,n)
    parity = MOD(ivec,i2)

    is_red = ALL(parity == 0) .OR. ALL(parity == 1)
    
    IF(is_red) THEN
      u(n) = relax_stencil_update(ndim,nsize,nshape,nstrides,bcs,wc,n,u,rhs)    
    ENDIF
    
  ENDDO  
  !$OMP END PARALLEL DO

  !
  ! Relax: Black
  !
  !$OMP PARALLEL DO PRIVATE(n,ivec,parity,is_red)
  DO n=1,nsize  

    ivec   = 0
    ivec   = lin2nd(ndim,nshape,n)
    parity = MOD(ivec,i2)

    is_red = ALL(parity == 0) .OR. ALL(parity == 1)
    
    IF(.NOT.is_red) THEN
      u(n) = relax_stencil_update(ndim,nsize,nshape,nstrides,bcs,wc,n,u,rhs)    
    ENDIF

  ENDDO  
  !$OMP END PARALLEL DO

  ! 
  ! Neumann correction
  !
  ! For pure Neumann problems, the mean must be subtracted off
  !  
  IF(ALL(bcs == "N")) THEN

    ! Compute mean
    u_mean = mean(nsize,u)
    
    ! Subtract  mean
    !
    !$OMP PARALLEL DO PRIVATE(n)
    DO n=1,nsize
      u(n) = u(n) - u_mean
    ENDDO
    !$OMP END PARALLEL DO
    
  ENDIF
  
END SUBROUTINE  

! ---------------------------------------------------------------------  
!
!>@name relax_stencil_update
!
!>@brief

FUNCTION relax_stencil_update(ndim,nsize,nshape,nstrides,bcs,wc,n,u,rhs) RESULT(u_new)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                       ,INTENT(IN) :: ndim,n,nsize
  INTEGER(IT) ,DIMENSION(ndim)      ,INTENT(IN) :: nshape,nstrides
  REAL(FP),DIMENSION(0:ndim)        ,INTENT(IN) :: wc
  CHARACTER(LEN=1),DIMENSION(ndim,2),INTENT(IN) :: bcs
  REAL(FP),DIMENSION(nsize)         ,INTENT(IN) :: u,rhs
  
  ! RETURN
  REAL(FP) :: u_new
  
  ! LOCAL
  INTEGER(IT),DIMENSION(ndim) :: ivec,bmask
  INTEGER(IT),DIMENSION(2)    :: dn 
  INTEGER(IT)                 :: i
  
  ! Compute ND index of centre point 
  ivec = lin2nd(ndim,nshape,n)
  
  ! Check if point lies on a boundary
  CALL boundary_mask(ndim,ivec,nshape,bmask)
  
  ! ===============
  ! DIRICHLET CHECK
  ! ===============
  !
  ! Check if the point at n lies on a Dirichlet boundary.
  ! No relaxation is performed at Dirichlet boundaries, so 
  ! the function returns the input u(n) there.
  ! 
  !
  IF(at_dirichlet_boundary(ndim,bmask,bcs)) THEN    
    u_new = u(n)
    RETURN
  ENDIF
  
  ! ============
  ! STENCIL SUM
  ! ===========
  !
  ! Perform stencil summation at point n. At Neumann boundaries,
  ! the stencil is modified to enforce grad(u).n = 0.
  !  
  u_new = 0    
  DO i=1,ndim  
  
    ! Compute change in +/- change in n. If the point is at a boundary,
    ! then the change is reversed to enforce the homogeneous
    ! Neumann condtion 
    !
    CALL stencil_stride(ndim,bmask,nstrides,i,dn)
                                  
    ! Add to summation    
    u_new = u_new + u(n+dn(1))*wc(i) + u(n+dn(2))*wc(i)
    
  ENDDO

  u_new = (u_new - rhs(n))*wc(0)
  
END FUNCTION

! ---------------------------------------------------------------------  
!
!>@name stencil_stride
!
!>@brief Compute the linear stride for stencil 
!
!>@details
!!
!!  Compute linear stride for stencil. Enforces homogenous
!!  Neumann boundary conditions at boundaries
!!

SUBROUTINE stencil_stride(ndim,bmask,nstrides,i,dn)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,i
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: bmask,nstrides

  ! OUTPUT
  INTEGER(IT),DIMENSION(2),INTENT(OUT) :: dn

  ! LOCAL
  INTEGER(IT),DIMENSION(2),PARAMETER :: di = [-1,+1]
  
  SELECT CASE(bmask(i))
    CASE(BND_NONE)
      dn = di*nstrides(i)
    CASE(BND_LOWER)
      dn = +nstrides(i)
    CASE(BND_UPPER)
      dn = -nstrides(i)
    CASE DEFAULT
      STOP "FATAL ERROR: MASK VALUE NOT DEFINED CFV3"
  END SELECT
      
END SUBROUTINE

END MODULE
