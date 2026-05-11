

PROGRAM MAIN

  USE NDSM_ROOT
  USE NDSM_POISSON
  USE NDSM_MULTIGRID_CORE
  
  IMPLICIT NONE

  INTEGER(IT)                     ,PARAMETER :: NDIM   = 3
  INTEGER(IT)                     ,PARAMETER :: NGRIDS = 2
  INTEGER(IT)                     ,PARAMETER :: NSHAPE(NDIM) = 2*[20,20,20]
  CHARACTER(LEN=1),DIMENSION(NDIM),PARAMETER :: BCFLAG = ["P","P","N"] 

  ! Solution parameters 
  INTEGER(IT),PARAMETER :: nmax_e = 1e8
  INTEGER(IT),PARAMETER :: nmax   = 1e6
  REAL(FP)   ,PARAMETER :: vc_tol = 1e-8
  REAL(FP)   ,PARAMETER :: ex_tol = 1e-8 
  REAL(FP)   ,PARAMETER :: ms     = 5
  
  ! ARRAYS
  REAL(FP),DIMENSION(:),ALLOCATABLE :: u,ue,rhs,r
  INTEGER(IT) :: nsize

  ! ERROR/CONV METRICS
  REAL(FP) :: E_max,r_max,E_mean,du_last
  
  ! MESH
  TYPE(MG_PTR),DIMENSION(NDIM) :: mesh
  REAL(FP) :: qmin = 0 
  REAL(FP) :: dq
 
  ! MULTIGRID HANDLE
  TYPE(MG_HANDLE) :: bvp  
  CHARACTER(LEN=1),DIMENSION(NDIM,2)  :: bcs

  ! MISC
  INTEGER(IT) :: ierr
  INTEGER(IT) :: i
  
  ! Create fine mesh
  DO i=1,NDIM
    ALLOCATE(mesh(i)%val(NSHAPE(i)))
    IF (BCFLAG(i)=="P") THEN    
      dq = REAL(1,FP)/REAL(nshape(i),FP)
    ELSE
      dq = REAL(1,FP)/REAL(nshape(i)-1,FP)    
    ENDIF      
    CALL create_mesh(dq,qmin,mesh(i)%val)
  ENDDO  
  nsize = PRODUCT(nshape)
  
  ! Allocate arrays
  ALLOCATE(u(nsize))
  ALLOCATE(ue(nsize))
  ALLOCATE(rhs(nsize))
  ALLOCATE(r(nsize))
  u   = 0
  rhs = 0
  r   = 0
  
  ! Test case 
  CALL test_case(ndim,nsize,nshape,mesh,BCFLAG,ue,rhs)
  
  DEBUG=.FALSE.

  ! Setup boundary flags 
  DO i=1,NDIM
    IF (BCFLAG(i) =="P" ) THEN
      mesh(i)%periodic = .TRUE.
      bcs(i,:)         = "P"
    ELSE
      mesh(i)%periodic = .FALSE.
      bcs(i,:)         = BCFLAG(i) 
    ENDIF
  ENDDO

  ! Set flags
  DO i=1,NDIM  
    bvp%copt(i)      = bcs(i,1)
    bvp%copt(NDIM+i) = bcs(i,2)
  ENDDO
  
      
  ! Setup BVP
  CALL new_mg_handle(bvp,ndim,nshape,NGRIDS,mesh,.TRUE.,nmax_e)
  
  ! Solution parameters
  bvp%ex_tol = ex_tol
  bvp%ms     = ms

  u = 0
  CALL solve_poisson_bvp(bvp,nsize,vc_tol,nmax,u,rhs,du_last,ierr) 
  
  ! Print messages
  IF (ierr == IERR_POISSON_SUCCESS) THEN
    PRINT *,"Return success"
  ENDIF
  PRINT *,"DU FINAL: ",du_last

  ! Compute residual of solution 
  CALL residual(ndim,nsize,nshape,mesh,bcs,u,rhs,r)
  
  ! Max error 
  E_max  = MAXVAL(ABS(u-ue))
  r_max  = MAXVAL(ABS(r))
  E_mean = SUM(ABS(u-ue))/REAL(SIZE(u),FP) 
    
  ! Report max error
  PRINT *,"Max error:  ",E_max  
  PRINT *,"Mean error: ",E_mean  
  PRINT *,"MAX RES:    ",r_max
    
CONTAINS

! ---------------------------------------------------------------------
!
!>@brief Test case with periodic and one neuman boundary

SUBROUTINE test_case(ndim,nsize,nshape,qvecs,bcflag,u,rhs)
  
  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                     ,INTENT(IN) :: ndim,nsize
  TYPE(MG_PTR),DIMENSION(ndim)    ,INTENT(IN) :: qvecs
  INTEGER(IT),DIMENSION(ndim)     ,INTENT(IN) :: nshape
  CHARACTER(LEN=1),DIMENSION(ndim),INTENT(IN) :: bcflag
  
  ! OUTPUT 
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: rhs,u
  
  ! LOCAL
  INTEGER(IT) :: i,n
  REAL(FP)   ,DIMENSION(ndim) :: Ln
  INTEGER(IT),DIMENSION(ndim) :: ivec
  REAL(FP) :: phase ,qi,ki,Lx,qc,kappa2
  REAL(FP),PARAMETER :: DPI = ACOS(REAL(-1,FP))
  
  REAL(FP),PARAMETER :: phase0 = 0.4D0
  
  ! Wave number 
  ki = 2*DPI
  
  ! Get size of domain 
  ! TODO: Fix for different size
  DO i=1,ndim
    Ln(i) = 1 !MAXVAL(qvecs(i)%val) - MINVAL(qvecs(i)%val)
  ENDDO
  
  !$OMP PARALLEL DO PRIVATE(ivec,qi,i,kappa2,phase)
  DO n=1,nsize

    ! Convert to [q1,q2,q3,...,qn]
    ivec = lin2nd(ndim,nshape,n)

    ! Compute kappa^2
    kappa2 = 0
    DO i=1,ndim
      kappa2 = kappa2 + ki**2
    ENDDO

    ! Compute product of COS overe dimensions
    rhs(n) = -kappa2
    u(n)   = 1
    DO i=1,ndim
    
      ! For periodic solution set phase to something non-zero
      phase = 0
      IF (bcflag(i)=="P") THEN
        phase = phase0
      ENDIF
      
      ! Get coordinate value 
      qi     = qvecs(i)%val(ivec(i))

      ! Set SIN/COS dependin on boundary conditions      
      IF (bcflag(i) == "D") THEN
        rhs(n) = rhs(n)*SIN(ki*qi/Ln(i) )
        u(n)   = u(n)*SIN(ki*qi/Ln(i)   )
      ELSE
        rhs(n) = rhs(n)*COS(ki*qi/Ln(i) + phase)
        u(n)   = u(n)*COS(ki*qi/Ln(i)   + phase)
      ENDIF

    ENDDO
    
  ENDDO
  !$OMP END PARALLEL DO
  
END SUBROUTINE
     
END PROGRAM
