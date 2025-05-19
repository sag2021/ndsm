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
! Module for computing the vector potential in 3D
!
MODULE NDSM_VECTOR_POTENTIAL

  USE NDSM_ROOT
  USE NDSM_POISSON

  IMPLICIT NONE

  PRIVATE

  ! IOPT VECTOR
  !
  INTEGER,PUBLIC,PARAMETER :: IOPT_LEN     = 16    !< Integer options/parameters vector length
  INTEGER,PUBLIC,PARAMETER :: IOPT_MS      = 0     !< Number of smoothing sweeps
  INTEGER,PUBLIC,PARAMETER :: IOPT_NCYCLES = 1     !< Max. number of V cycles (typically set high) 
  INTEGER,PUBLIC,PARAMETER :: IOPT_FACE1   = 2     !< Approach to solving BVP
  INTEGER,PUBLIC,PARAMETER :: IOPT_IERR    = 3     !< Error solving BVP
  INTEGER,PUBLIC,PARAMETER :: IOPT_FLXCRL  = 4     !< Order in which flux correction is applied
  INTEGER,PUBLIC,PARAMETER :: IOPT_DEBUG   = 5     !< Prints more information when running
  INTEGER,PUBLIC,PARAMETER :: IOPT_DUMAX   = 6     !< If set, use max change as convergence metric
  INTEGER,PUBLIC,PARAMETER :: IOPT_TRUE    = 1
  INTEGER,PUBLIC,PARAMETER :: IOPT_FALSE   = 0

  ! ROPT VECTOR
  !
  INTEGER,PUBLIC,PARAMETER :: ROPT_VTOL = 0        !< V-cycle tolerance
  INTEGER,PUBLIC,PARAMETER :: ROPT_CTOL = 1        !< Tolerance on coarsest mesh
  INTEGER,PUBLIC,PARAMETER :: ROPT_TIM  = 2        !< Time taken
  INTEGER,PUBLIC,PARAMETER :: ROPT_LEN  = IOPT_LEN !< Real/float options/parameters vector length

  ! Base grid level
  REAL(FP),PRIVATE,PARAMETER :: BASE_GRID = 2

  ! ===============
  ! MAPPING ARRAYS:
  ! ===============
  !
  ! The six boundaries are identified with an index in 
  ! the range [1,6]. These arrays map this ID to other useful
  ! values
  !
  ! imap_cp: Maps the boundary ID to the component index that is 
  !          normal to the boundary, i.e. if b(:,:,:,i) is normal 
  !          at boundary j, then imap_cp(j) = i
  !
  ! imap_ul: Maps the boundary ID to either 1 or 2 depending on 
  !          whether the boundary is an upper or lower boundary, e.g.
  !          imap_ul(6) = 2 because ID:6 is the top boundary
  !
  ! imap_nc: 
  !
  !
  INTEGER(IT),DIMENSION(6)  ,PARAMETER :: imap_ul = [1,2,1,2,1,2]                                !< Upper-lower map
  INTEGER(IT),DIMENSION(6)  ,PARAMETER :: imap_cp = [1,1,2,2,3,3]                                !< Component map
  INTEGER(IT),DIMENSION(2,6),PARAMETER :: imap_nc = RESHAPE([2,3,2,3,1,3,1,3,1,2,1,2],[2,6])

  ! ============
  ! UNIT VECTORS
  ! ============
  !
  ! Unit vectors for the sides of the box
  ! 

  ! First tangent vector  
  !
  REAL(FP),DIMENSION(3,6) :: tvecs1 = RESHAPE([0,1,0,&
                                               0,1,0,&
                                               1,0,0,&
                                               1,0,0,&
                                               1,0,0,&
                                               1,0,0],[3,6])
  ! Second tangent vector
  ! 
  REAL(FP),DIMENSION(3,6) :: tvecs2 = RESHAPE([0,0,1,&
                                               0,0,1,&
                                               0,0,1,&
                                               0,0,1,&
                                               0,1,0,&
                                               0,1,0],[3,6])
  
  ! Normal vector
  !
  REAL(FP),DIMENSION(3,6) :: nvecs = RESHAPE([1,0,0,&
                                              1,0,0,&
                                              0,1,0,&
                                              0,1,0,&
                                              0,0,1,&
                                              0,0,1],[3,6])

  ! Public
  PUBLIC :: compute_vector_potential

CONTAINS

! ---------------------------------------------------------------------
! 
!
!>@name compute_vector_potential
!
!>@brief Computes the vector potential

SUBROUTINE compute_vector_potential(nshape,iopt,ropt,mesh,Apot,B)

  IMPLICIT NONE

  ! SUBNAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "compute_vector_potential"

  ! INPUT
  INTEGER(IT) ,DIMENSION(4),INTENT(IN) :: nshape
  TYPE(MG_PTR),DIMENSION(3),INTENT(IN) :: mesh      !< 3D mesh

  ! INPUT/OUTPUT
  INTEGER(IT),DIMENSION(0:IOPT_LEN-1)                        ,INTENT(INOUT) :: iopt !< Integer options
  REAL(FP)   ,DIMENSION(0:IOPT_LEN-1)                        ,INTENT(INOUT) :: ropt !< Real options
  REAL(FP),DIMENSION(nshape(1),nshape(2),nshape(3),nshape(4)),INTENT(INOUT) :: B

  ! OUTPUT
  REAL(FP),DIMENSION(nshape(1),nshape(2),nshape(3),nshape(4)),INTENT(OUT)   :: Apot

  ! ----------
  ! LOCAL VARS
  ! ----------

  ! 3D MESH VARIABLES
  REAL(FP)    ,DIMENSION(3)   :: dq          !< Mesh spacing
  REAL(FP)    ,DIMENSION(3)   :: Lq          !< Extent of mesh
  REAL(FP)    ,DIMENSION(6)   :: Aq
  
  ! BOUNDARY MESH
  TYPE(MG_PTR),DIMENSION(2,6) :: mesh_bn     !< Mesh on 6 boundaries

  ! BOUNDARY CONDITIONS: B FIELD
  TYPE(MG_PTR),DIMENSION(6)   :: bn              !< Boundary conditions B.n on 6 faces
  INTEGER(IT) ,DIMENSION(2,6) :: nshape_bn       !< Shape of Bn array 
  INTEGER(IT) ,DIMENSION(6)   :: nsize_bn        !< Total size of Bn arrays

  ! BOUNDARY CONDITIONS: A FIELD
  TYPE(MG_PTR),DIMENSION(2,6) :: At              !< Boundary conditions for vector potential At = -grad(chi) x n 

  ! AUXILLARY CHI FIELD  
  TYPE(MG_PTR),DIMENSION(6)   :: chi             !< Scalar field: div(grad(chi)) = bn

  ! BOUNDARY-VALUE HANDLES
  TYPE(MG_HANDLE) :: bvp
  TYPE(MG_HANDLE) :: bvp_chi

  ! MAGNETIC FLUX 
  REAL(FP),DIMENSION(6) :: phi_b        !< Flux over boundary i
  REAL(FP)              :: phi_b_total  !< Flux over all six boundaries
  REAL(FP)              :: du_last
  
  ! TEMP DATA
  REAL(FP),DIMENSION(:,:,:,:),ALLOCATABLE :: Ac

  ! LOOP INDICES
  INTEGER(IT) :: i,j
 
  ! MISC
  INTEGER(IT) :: lb,ndim,ngrids,nmin,nsize,ierr
  INTEGER(IT),PARAMETER :: i1 = 1
  LOGICAL :: USE_DU_MAX 


  ! Conv. flag 
  USE_DU_MAX = (IOPT(IOPT_DUMAX) == IOPT_TRUE)

  ! ====================
  ! GET MESH VARS
  ! ====================

  ! Get mesh extent
  DO i=1,SIZE(mesh)
    Lq(i) = MAXVAL(mesh(i)%val) - MINVAL(mesh(i)%val)
  ENDDO

  ! Get mesh spacing. Spacing is assumed to be uniform.
  ! The spacing is computed as the difference between the first
  ! and last elements. Default bounds are not assumed. The 
  ! mesh must have at least two elements.
  !  
  DO i=1,SIZE(mesh)

    ! Check that mesh is larger than two points
    IF(SIZE(mesh(i)%val) < 2) THEN
      iopt(IOPT_IERR) = 1
      RETURN
    ENDIF

    lb    = LBOUND(mesh(i)%val,1)
    dq(i) = mesh(i)%val(lb+1) - mesh(i)%val(lb)
        
  ENDDO

  ! ====================
  ! EXTRACT BCS FOR B
  ! ====================
  !
  ! Extract the normal components of B at the boundary
  !
  !
  ! The integer array imap_ic maps the linear index i
  ! to the correct vector index for nshape
  !
  ! 1 -> [2,3]  ; 4 -> [1,3]
  ! 2 -> [2,3]  ; 5 -> [1,2]
  ! 3 -> [1,3]  ; 6 -> [1,2]
  !
  ! Shape array on each boundary
  !
  ! 1: [ny,nz]
  ! 2: [ny,nz]
  ! 3: [nx,nz]
  ! 4: [nx,nz]
  ! 5: [nx,ny]
  ! 6: [nx,ny]
  !
  !
  DO i=1,SIZE(bn)
    nshape_bn(:,i) = nshape(imap_nc(:,i))
  ENDDO
  
  ! Set of each mesh
  DO i=1,SIZE(nsize_bn)
    nsize_bn(i) = PRODUCT(nshape_bn(:,i))
  ENDDO
  
  ! Allocate memory to hold boundary conditions 
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Allocate memory to hold boundary conditions...")
  DO i=1,SIZE(bn)  
    ALLOCATE(bn(i)%val(nsize_bn(i)))
  ENDDO
  
  ! Allocate mesh for boundary faces
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Allocate mesh for boundary faces...")
  DO j=1,SIZE(bn)
    DO i=1,2
      ALLOCATE(mesh_bn(i,j)%val( nshape_bn(i,j) ))
      mesh_bn(i,j)%val = mesh(imap_nc(i,j))%val
    ENDDO  
  ENDDO

  
  !
  ! Extract boundary conditions
  !
  ! The map, imap_ul, is used to select the correct layer.
  ! It returns either 1 or two depending on the value of i.
  !
  ! The map, imap_cp, maps i to the correct component for boundary
  ! i (either 1,2 or 3). The component map (imap_cp) also 
  ! maps to the dimension
  !
  !
  DO i=1,SIZE(bn)
  
    ! Pick correct layer
    IF(imap_ul(i) == 1) j = 1
    IF(imap_ul(i) == 2) j = nshape(imap_cp(i))    
  
    ! Extract boundary
    CALL extract_bn(nsize_bn(i),nshape(:3),imap_cp(i),j,     &
                    b(:,:,:,imap_cp(i)), bn(i)%val, dir=+i1 )
                    
  ENDDO
  
  ! ==============
  ! COMPUTE FLUXES
  ! ==============
 
  ! Compute fluxes 
  DO i=1,SIZE(phi_b)
    phi_b(i) = trapz_2D(nshape_bn(1,i), &
                        nshape_bn(2,i), &
                        dq(1),          &
                        dq(2),          & 
                        bn(i)%val)
  ENDDO

  ! ====================
  ! SOLVE 2D BVPS
  ! ====================
  !
  ! Solve laplace(chi) = bn on each boundary.
  !
  ! Homogenous Neumman boundary conditions are imposed 
  ! on each boundary for chi (the edges of the 3D volume).
  !

  ! Allocate memory for chi
  DO i=1,SIZE(chi)
    ALLOCATE(chi(i)%val(nsize_bn(i)))    
  ENDDO
  
  ! Compute area of each face
  
  Aq(1) = Lq(2)*Lq(3)
  Aq(2) = Lq(2)*Lq(3)
 
  Aq(3) = Lq(1)*Lq(3)
  Aq(4) = Lq(1)*Lq(3)

  Aq(5) = Lq(1)*Lq(2)
  Aq(6) = Lq(1)*Lq(2)  
  
  !
  ! Solve BVP on each boundary
  !
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Solve BVP on each boundary...")
  DO i=1,SIZE(bn)
  
    ! Compute number of grids
    nmin   = MINVAL(nshape_bn(:,i))
    ngrids = FLOOR(LOG(nmin/BASE_GRID)/LOG(REAL(2,FP)))
  
    ! Initialize chi array to zero
    chi(i)%val = 0
  
    ! Flux balance boundary
    bn(i)%val = bn(i)%val - phi_b(i)/Aq(i)
        
    ! Construct new 2D handle
    ndim = 2    
    CALL new_mg_handle(bvp,ndim,nshape_bn(:,i),ngrids,mesh_bn(:,i),USE_DU_MAX)
    
    ! Set option values
    bvp%ms                = IOPT(IOPT_MS  )  ! Smoothing sweeps
    bvp%ex_tol            = ROPT(ROPT_CTOL)  ! Exact solution tolerance
    bvp%copt(1:SIZE(bn))  = "N"              ! Boundary conditions
    
    ! Solve 
    CALL solve_poisson_bvp(bvp,nsize_bn(i),ropt(ROPT_VTOL),iopt(IOPT_NCYCLES),chi(i)%val,bn(i)%val,du_last,ierr) 
  
    ! Delete MG handle
    CALL delete_mg_handle(bvp)
    
  ENDDO

  ! =========================
  ! COMPUTE BOUNDARY CONDITIONS FOR AT
  ! ==========================

  ! Allocate memory to hold boundary conditions for vector 
  ! potential
  DO j=1,SIZE(bn)  
    DO i=1,2
      ALLOCATE(At(i,j)%val(nsize_bn(j)))
    ENDDO
  ENDDO
  
  !
  ! Compute vector potential boundary conditions
  !
  ! At = - grad(chi) x n
  !
  !
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Compute vector potential boundary conditions...")

  DO i=1,SIZE(chi)
  
    ! Initialize to zero
    At(1,i)%val = 0
    At(2,i)%val = 0
    
    ! Maps to component
    j = imap_cp(i)
    
    ! Compute At at each boundary
    CALL  compute_At_bcs(nshape_bn(:,i),chi(i)%val,dq(j)      ,tvecs1(:,i), &
                         tvecs2(:,i)   ,nvecs(:,i),At(1,i)%val,At(2,i)%val)    
  ENDDO

  ! =====================
  ! SOLVE 3D BVP
  ! =====================

  IF(DEBUG) CALL debug_msg(THIS_SUB,"Solve BVP 3D...")

  ndim  = 3
  nsize = PRODUCT(nshape(:3))

      
  !
  ! Compute with Bn nonzero on only one face at a time
  !
  SELECT CASE(IOPT_FACE1) 
  
    CASE(1)
    !
    ! Compute A one face at a time and then sum the results
    !    
    PRINT *,"FLAG SET: FACE1"
    ALLOCATE(Ac(nshape(1),nshape(2),nshape(3),3)) 

    ! Initialize A to zero
    Apot = 0
  
    ! Sum contribution from solving BVP for each of the six faces
    DO i=1,6 
      Ac = 0  
      CALL  solve_6faces(i,ROPT,IOPT,mesh,nsize_bn,nshape,At,Ac)  
      Apot = Apot + Ac
    ENDDO

    DEALLOCATE(Ac)
   
    CASE DEFAULT
    !
    ! Default: Compute A from all six faces at once
    !
    
    CALL solve(ROPT,IOPT,mesh,nsize_bn,nshape,At,Apot)
    
  END SELECT

  ! ===============
  ! COMPUTE CURL
  ! ===============
  !
  ! Compute curl and perform flux correction. The IOPT_ACBC
  ! flag determines the order in which this is performed.
  !
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Compute B = curl(B) and flux correction...")

  SELECT CASE(IOPT(IOPT_FLXCRL))
  
    CASE(1)
    !
    ! Compute curl and then perform flux correction. In this
    ! case the magnetic field is computed from the numerical
    ! curl of the ucorrected vector potential. Analytic flux
    ! correction terms are then added onto both the vector
    ! potential and the magnetic field.
    
    PRINT *,"FLAG SET: FLXCRL"
    CALL curl(dq,Apot,B)    
    CALL add_flux_balance_fields(nshape,mesh,phi_b,B,Apot)
    
    CASE DEFAULT
    !
    ! Perform flux correction and then compute the curl. In this
    ! case the magnetic field is computed from the numerical curl
    ! of the flux corrected vector potential. The correction to
    ! the vector potential is analytic.
    
    CALL add_flux_balance_fields(nshape,mesh,phi_b,B,Apot)
    CALL curl(dq,Apot,B)
    
  END SELECT
    
  ! Set error flag
  iopt(IOPT_IERR) = ierr

  ! ==================
  ! FREE MEMORY
  ! ==================

  ! Deallocate 
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Deallocate memory...")

  DO i=1,SIZE(bn)
    DEALLOCATE(bn(i)%val)
    DO j=1,2
      DEALLOCATE(mesh_bn(j,i)%val) 
      DEALLOCATE(At(j,i)%val)
    ENDDO
  ENDDO

END SUBROUTINE

! ---------------------------------------------------------------------

SUBROUTINE solve_6faces(cdim,ROPT,IOPT,mesh,nsize_bn,nshape,At,Ac)

  ! INPUT
  INTEGER(IT)                         ,INTENT(IN) :: cdim
  INTEGER(IT) ,DIMENSION(:)           ,INTENT(IN) :: nsize_bn,nshape
  REAL(FP)    ,DIMENSION(0:ROPT_LEN-1),INTENT(IN) :: ROPT
  INTEGER(IT) ,DIMENSION(0:IOPT_LEN-1),INTENT(IN) :: IOPT
  TYPE(MG_PTR),DIMENSION(:)           ,INTENT(IN) :: mesh
  
  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: Ac
  TYPE(MG_PTR),DIMENSION(:,:)            ,INTENT(INOUT) :: At
  
  ! LOCAL
  TYPE(MG_HANDLE)                       :: bvp
  REAL(FP),DIMENSION(:,:,:),ALLOCATABLE :: rhs
  INTEGER(IT)                           :: max_vcycles,ngrids,nmin,ndim,nsize,ierr
  REAL(FP)                              :: vcycle_tol,coarse_tol,du_last
  LOGICAL                               :: USE_DU_MAX

  ! LOCAL INT. PARAMS
  INTEGER(IT),PARAMETER  :: i1 = 1
  INTEGER(IT),PARAMETER  :: i2 = 2
  INTEGER(IT),PARAMETER  :: i3 = 3
  
  ! Conv. flag 
  USE_DU_MAX = (IOPT(IOPT_DUMAX) == IOPT_TRUE)

  ! Set dimensions and size parameter
  ndim = 3
  nsize = PRODUCT(nshape(:3))
  
  ! Compute number of grids to use
  nmin   = MINVAL(nshape(:3))
  ngrids = FLOOR(LOG(nmin/BASE_GRID)/LOG(2.D+0))
  
  ! Unpack
  max_vcycles = IOPT(IOPT_NCYCLES)
  vcycle_tol  = ROPT(ROPT_VTOL)
  coarse_tol  = ROPT(ROPT_CTOL)
  
  ! Allocate RHS
  ALLOCATE(rhs(nshape(1),nshape(2),nshape(3)))
  rhs = 0
  
  ! ===============
  ! AX
  ! ===============
  
  IF(cdim == 3) CALL extract_bn(nsize_bn(3),nshape(:3),i2,i1       ,Ac(:,:,:,1),At(1,3)%val,dir=-i1)
  IF(cdim == 4) CALL extract_bn(nsize_bn(4),nshape(:3),i2,nshape(2),Ac(:,:,:,1),At(1,4)%val,dir=-i1)  
  IF(cdim == 5) CALL extract_bn(nsize_bn(5),nshape(:3),i3,i1       ,Ac(:,:,:,1),At(1,5)%val,dir=-i1)
  IF(cdim == 6) CALL extract_bn(nsize_bn(6),nshape(:3),i3,nshape(3),Ac(:,:,:,1),At(1,6)%val,dir=-i1)
    
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms                = IOPT(IOPT_MS)
  bvp%ex_tol            = coarse_tol
  bvp%copt(1:6)         = ["N","D","D","N","D","D"]  
  CALL solve_poisson_bvp(bvp,nsize, vcycle_tol,max_vcycles,Ac(:,:,:,1),rhs,du_last,ierr) 
  CALL delete_mg_handle(bvp)
  
  ! ===============
  ! AY
  ! ===============

  IF(cdim == 1) CALL extract_bn(nsize_bn(1),nshape(:3),i1,i1       ,Ac(:,:,:,2),At(1,1)%val,dir=-i1)
  IF(cdim == 2) CALL extract_bn(nsize_bn(2),nshape(:3),i1,nshape(1),Ac(:,:,:,2),At(1,2)%val,dir=-i1)  
  IF(cdim == 5) CALL extract_bn(nsize_bn(5),nshape(:3),i3,i1       ,Ac(:,:,:,2),At(2,5)%val,dir=-i1)
  IF(cdim == 6) CALL extract_bn(nsize_bn(6),nshape(:3),i3,nshape(3),Ac(:,:,:,2),At(2,6)%val,dir=-i1)
  
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms         = IOPT(IOPT_MS)
  bvp%ex_tol     = coarse_tol
  bvp%copt(1:6)  = ["D","N","D","D","N","D"]  
  CALL solve_poisson_bvp(bvp,nsize,vcycle_tol,max_vcycles,Ac(:,:,:,2),rhs,du_last,ierr)   
  CALL delete_mg_handle(bvp)

  ! ===============
  ! AZ
  ! ===============

  IF(cdim == 1) CALL extract_bn(nsize_bn(1),nshape(:3),i1,i1       ,Ac(:,:,:,3),At(2,1)%val,dir=-i1)
  IF(cdim == 2) CALL extract_bn(nsize_bn(2),nshape(:3),i1,nshape(1),Ac(:,:,:,3),At(2,2)%val,dir=-i1)
  IF(cdim == 3) CALL extract_bn(nsize_bn(3),nshape(:3),i2,i1       ,Ac(:,:,:,3),At(2,3)%val,dir=-i1)
  IF(cdim == 4) CALL extract_bn(nsize_bn(4),nshape(:3),i2,nshape(2),Ac(:,:,:,3),At(2,4)%val,dir=-i1)
  
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms                = 5
  bvp%ex_tol            = coarse_tol
  bvp%copt(1:6)         = ["D","D","N","D","D","N"]  
  CALL solve_poisson_bvp(bvp,nsize,vcycle_tol,max_vcycles,Ac(:,:,:,3),rhs,du_last,ierr) 
  CALL delete_mg_handle(bvp) 

END SUBROUTINE

! ----------------------------------------------------------------------

SUBROUTINE solve(ROPT,IOPT,mesh,nsize_bn,nshape,At,Ac)

  ! INPUT
  INTEGER(IT) ,DIMENSION(:)           ,INTENT(IN) :: nsize_bn,nshape
  REAL(FP)    ,DIMENSION(0:ROPT_LEN-1),INTENT(IN) :: ROPT
  INTEGER(IT) ,DIMENSION(0:IOPT_LEN-1),INTENT(IN) :: IOPT
  TYPE(MG_PTR),DIMENSION(:)           ,INTENT(IN) :: mesh
  
  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(:,:,:,:),INTENT(INOUT) :: Ac
  TYPE(MG_PTR),DIMENSION(:,:),INTENT(INOUT) :: At
  
  ! LOCAL
  TYPE(MG_HANDLE)                       :: bvp
  REAL(FP),DIMENSION(:,:,:),ALLOCATABLE :: rhs
  INTEGER(IT)                           :: max_vcycles,ngrids,nmin,ndim,nsize,i,k,cp,ierr
  REAL(FP)                              :: vcycle_tol,coarse_tol,du_last
  CHARACTER(LEN=1),DIMENSION(6)         :: bcs
  LOGICAL                               :: USE_DU_MAX
  
  ! LOCAL INT. PARAMS
  INTEGER(IT),PARAMETER  :: i1 = 1
  INTEGER(IT),PARAMETER  :: i2 = 2
  INTEGER(IT),PARAMETER  :: i3 = 3

  ! Conv. flag 
  USE_DU_MAX = (IOPT(IOPT_DUMAX) == IOPT_TRUE)

  ! Set dimensions and size parameter
  ndim = 3
  nsize = PRODUCT(nshape(:3))
  
  ! Compute number of grids to use
  nmin   = MINVAL(nshape(:3))
  ngrids = FLOOR(LOG(nmin/BASE_GRID)/LOG(REAL(2,FP)))
  
  ! Unpack
  max_vcycles = IOPT(IOPT_NCYCLES)
  vcycle_tol  = ROPT(ROPT_VTOL)
  coarse_tol  = ROPT(ROPT_CTOL)
  
  ! Allocate RHS
  ALLOCATE(rhs(nshape(1),nshape(2),nshape(3)))
  rhs = 0
  
  ! ===============
  ! AX
  ! ===============
  
  CALL extract_bn(nsize_bn(3),nshape(:3),i2,i1       ,Ac(:,:,:,1),At(1,3)%val,dir=-i1)
  CALL extract_bn(nsize_bn(4),nshape(:3),i2,nshape(2),Ac(:,:,:,1),At(1,4)%val,dir=-i1)  
  CALL extract_bn(nsize_bn(5),nshape(:3),i3,i1       ,Ac(:,:,:,1),At(1,5)%val,dir=-i1)
  CALL extract_bn(nsize_bn(6),nshape(:3),i3,nshape(3),Ac(:,:,:,1),At(1,6)%val,dir=-i1)
    
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms                = IOPT(IOPT_MS)
  bvp%ex_tol            = coarse_tol
  bvp%copt(1:6)         = ["N","D","D","N","D","D"]  
  CALL solve_poisson_bvp(bvp,nsize, vcycle_tol,max_vcycles,Ac(:,:,:,1),rhs,du_last,ierr) 
  CALL delete_mg_handle(bvp)
  
  ! ===============
  ! AY
  ! ===============

  CALL extract_bn(nsize_bn(1),nshape(:3),i1,i1       ,Ac(:,:,:,2),At(1,1)%val,dir=-i1)
  CALL extract_bn(nsize_bn(2),nshape(:3),i1,nshape(1),Ac(:,:,:,2),At(1,2)%val,dir=-i1)  
  CALL extract_bn(nsize_bn(5),nshape(:3),i3,i1       ,Ac(:,:,:,2),At(2,5)%val,dir=-i1)
  CALL extract_bn(nsize_bn(6),nshape(:3),i3,nshape(3),Ac(:,:,:,2),At(2,6)%val,dir=-i1)
  
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms         = IOPT(IOPT_MS)
  bvp%ex_tol     = coarse_tol
  bvp%copt(1:6)  = ["D","N","D","D","N","D"]  
  CALL solve_poisson_bvp(bvp,nsize,vcycle_tol,max_vcycles,Ac(:,:,:,2),rhs,du_last,ierr)   
  CALL delete_mg_handle(bvp)

  ! ===============
  ! AZ
  ! ===============

  CALL extract_bn(nsize_bn(1),nshape(:3),i1,i1       ,Ac(:,:,:,3),At(2,1)%val,dir=-i1)
  CALL extract_bn(nsize_bn(2),nshape(:3),i1,nshape(1),Ac(:,:,:,3),At(2,2)%val,dir=-i1)
  CALL extract_bn(nsize_bn(3),nshape(:3),i2,i1       ,Ac(:,:,:,3),At(2,3)%val,dir=-i1)
  CALL extract_bn(nsize_bn(4),nshape(:3),i2,nshape(2),Ac(:,:,:,3),At(2,4)%val,dir=-i1)
  
  CALL new_mg_handle(bvp,NDIM,nshape,NGRIDS,mesh,USE_DU_MAX)
  bvp%ms                = 5
  bvp%ex_tol            = coarse_tol
  bvp%copt(1:6)         = ["D","D","N","D","D","N"]  
  CALL solve_poisson_bvp(bvp,nsize,vcycle_tol,max_vcycles,Ac(:,:,:,3),rhs,du_last,ierr) 
  CALL delete_mg_handle(bvp) 
    
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name extract_bn
!
!>@brief Extract boundary conditions

SUBROUTINE extract_bn(nsize,nshape,cdim,clay,b,bc,dir)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)             ,INTENT(IN) :: nsize,cdim,clay,dir
  INTEGER(IT),DIMENSION(:),INTENT(IN) :: nshape
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(INOUT) :: bc
  REAL(FP),DIMENSION(:,:,:),INTENT(INOUT) :: b

  ! LOCAL
  INTEGER(IT) :: lb(3),ub(3)
  INTEGER(IT) :: i,j,k,n
  
  ! Default bounds
  lb = 1
  ub = nshape
  
  ! Degen. layer
  lb(cdim) = clay
  ub(cdim) = clay
  
  ! Linear bounds
  n=1
  
  DO k=lb(3),ub(3)
    DO j=lb(2),ub(2)
      DO i=lb(1),ub(1)
      
        ! Bounds check
        IF(n>nsize) THEN
          CALL error_msg("n out of bounds","extract_bn","n")
        ENDIF
        
        IF(dir == +1) bc(n)    = b(i,j,k)
        IF(dir == -1) b(i,j,k) = bc(n)
        n = n + 1

      ENDDO
    ENDDO
  ENDDO  
  
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name curl
!
!>@brief Compute B = curl(A)
!
!>@details
!!
!! Compute curl(A) using CFITX derive module 
!!
!!
!
!>@param[in] nshape Shape of A and B
!
SUBROUTINE curl(dq,A,B)

  IMPLICIT NONE
  
  ! INPUT
  REAL(FP),DIMENSION(3)      ,INTENT(IN) :: dq
  REAL(FP),DIMENSION(:,:,:,:),INTENT(IN) :: A
  
  ! OUTPUT
  REAL(FP),DIMENSION(:,:,:,:),INTENT(OUT) :: B
  
  ! LOCAL
  INTEGER(IT),DIMENSION(3) :: ivec,dir
  INTEGER(IT),DIMENSION(4) :: nshape
  INTEGER(IT)              :: i,j,k,ndim
  REAL(FP)                 :: dAx_dy,dAx_dz
  REAL(FP)                 :: dAy_dx,dAy_dz
  REAL(FP)                 :: dAz_dx,dAz_dy
  
  ! Set dimensions/shape
  ndim   = 3
  nshape = SHAPE(b)

  ! Directions
  dir  = [1,2,3]
    
  !$OMP PARALLEL DO PRIVATE(i,j,k,ivec,dAx_dy,dAx_dz,dAy_dx,dAy_dz,dAz_dx,dAz_dy)
  DO k=1,nshape(3)
    DO j=1,nshape(2)
      DO i=1,nshape(1)
      
        ! Pack integers
        ivec = [i,j,k] 

        ! Compute derivatives
        dAx_dy = derivq(ndim,nshape,ivec,dq,dir(2),A(:,:,:,1))
        dAx_dz = derivq(ndim,nshape,ivec,dq,dir(3),A(:,:,:,1))
        dAy_dx = derivq(ndim,nshape,ivec,dq,dir(1),A(:,:,:,2))
        dAy_dz = derivq(ndim,nshape,ivec,dq,dir(3),A(:,:,:,2))
        dAz_dx = derivq(ndim,nshape,ivec,dq,dir(1),A(:,:,:,3))
        dAz_dy = derivq(ndim,nshape,ivec,dq,dir(2),A(:,:,:,3))
      
        ! Magentic field components
        B(i,j,k,1) = dAz_dy - dAy_dz
        B(i,j,k,2) = dAx_dz - dAz_dx
        B(i,j,k,3) = dAy_dx - dAx_dy
    
      ENDDO
    ENDDO
  ENDDO  
  !$OMP END PARALLEL DO 

END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name derivq
!
!>@brief Compute derivative in dimension dir
!
!>@details
!!
!! Uses second-order finite differencing to compute
!! du/dq at a point ivec = [i,j,k,..]
!!

FUNCTION derivq(ndim,nshape,ivec,dq,dir,u) 

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,dir
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nshape,ivec
  REAL(FP)   ,DIMENSION(ndim),INTENT(IN) :: dq 
  REAL(FP)   ,DIMENSION(*)   ,INTENT(IN) :: u

  ! RETURN
  REAL(FP) :: derivq 
  
  ! LOCAL
  INTEGER(IT)               :: nsize,nc,i,n
  INTEGER(IT) ,DIMENSION(3) :: stencil,nstrides
  REAL(FP)    ,DIMENSION(3) :: wc
  REAL(FP)    ,PARAMETER    :: inv2 = REAL(1,FP)/REAL(2,FP)

  ! Compute size and strides
  nsize    = PRODUCT(nshape)
  nstrides = shape_to_strides(ndim,nshape)

  ! Get linear index
  n = nd2lin(ndim,nshape,ivec)
  
  ! Compute stencil and weights
  IF(ivec(dir) == 1) THEN
    stencil = nstrides(dir)*[+0,+1,+2]
    wc      = [-3,+4,-1]*inv2/dq(dir)
    nc      = 3 
  ELSEIF(ivec(dir) == nshape(dir)) THEN
    stencil = nstrides(dir)*[-0,-1,-2]  
    wc      = [+3,-4,+1]*inv2/dq(dir)
    nc      = 3
  ELSE
    stencil(1:2) = nstrides(dir)*[-1,+1]
    wc(1:2)      = [-1,+1]*inv2/dq(dir)          ! 1/(2*dq)
    nc           = 2
  ENDIF
 
  ! Perform summation over stencil
  derivq = 0
  DO i=1,nc
    derivq = derivq + u(n+stencil(i))*wc(i)  
  ENDDO

END FUNCTION

! ----------------------------------------------------------------------
!
!>@name add_flux_balance_fields
!
!>@brief Add flux correction fields
!
SUBROUTINE add_flux_balance_fields(nshape,mesh,phi_b,B,A)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)  ,DIMENSION(3)      ,INTENT(IN) :: nshape
  TYPE(MG_PTR),DIMENSION(3),TARGET,INTENT(IN) :: mesh
  REAL(FP)    ,DIMENSION(:)       ,INTENT(IN) :: phi_b

  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(:,:,:,:),INTENT(INOUT) :: b,A
  
  ! LOCAL
  INTEGER(IT)           :: i,j,k
  REAL(FP),DIMENSION(3) :: Lq,bc,Ac,gamma,A1_l,A2_l,A3_l,A_c
  REAL(FP)              :: Vq
  REAL(FP),PARAMETER    :: o0   = REAL(0,FP)
  REAL(FP),PARAMETER    :: inv3 = REAL(1,FP)/REAL(3,FP)
    
  ! LOCAL POINTER ALIAS
  ! Makes code a bit easier to read
  REAL(FP),DIMENSION(:),POINTER :: x,y,z

  ! Compute the extent of the mesh
  DO i=1,3
    Lq(i) = MAXVAL(mesh(i)%val) - MINVAL(mesh(i)%val)
  ENDDO  
  
  ! Compute volume
  Vq = PRODUCT(Lq)
      
  ! Remap. Makes it easier to read
  x => mesh(1)%val
  y => mesh(2)%val
  z => mesh(3)%val
      
  ! Flux balances
  gamma(1) = (phi_b(2) - phi_b(1))/Vq
  gamma(2) = (phi_b(4) - phi_b(3))/Vq
  gamma(3) = (phi_b(6) - phi_b(5))/Vq 
      
  !$OMP PARALLEL DO PRIVATE(i,j,k,bc,Ac,A1_l,A2_l,A3_l,A_c)
  DO k=1,nshape(3)
    DO j=1,nshape(2)
      DO i=1,nshape(1)
      
        ! Magnetic field
        bc(1) = gamma(1)*x(i) + phi_b(1)*Lq(1)/Vq
        bc(2) = gamma(2)*y(j) + phi_b(3)*Lq(2)/Vq
        bc(3) = gamma(3)*z(k) + phi_b(5)*Lq(3)/Vq

        ! Vector potential corresponding to linear b terms
        A1_l = [-gamma(3)*y(j)*z(k),o0                       ,+gamma(1)*x(i)*y(j)]
        A2_l = [+gamma(2)*z(k)*y(j),-gamma(1)*x(i)*z(k)      ,o0                 ]        
        A3_l = [o0                 ,+gamma(3)*x(i)*z(k)      ,-gamma(2)*x(i)*y(j)]
        
        ! Vector potential corresponding to constant b terms
        A_c = -[phi_b(5)*Lq(3)*y(j), &
                phi_b(1)*Lq(1)*z(k), &
                phi_b(3)*Lq(2)*x(i)]/Vq
                
        ! Add corrections
        b(i,j,k,:) = b(i,j,k,:) + bc 
        A(i,j,k,:) = A(i,j,k,:) + A_c + inv3*(A1_l+A2_l+A3_l) 
      
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
    
END SUBROUTINE

! ----------------------------------------------------------------------
!
!>@name compute_At
!
!>@brief Computed At = -grad(chi) x n 
!
!>@details
!!
!! Computes the boundary conditions for the transverse component
!! of the vector potential on a particular boundary
!!
!! The function chi is computed from the normal component of the
!! magnetic field at the boundary as
!!     div(grad(u)) = bn
!!
!!
!
!>@param[in]  nshape  Shape of boundary. [n1,n2]
!>@param[in]  chi     Function whose perpendicualar gradient is used to compute At
!>@param[in]  dq      Mesh spacing
!>@param[in]  tvec1   First unit tangent vector in boundary
!>@param[in]  tvec2   Second unit tangent vector in boundary
!>@param[in]  nvec    Nornaml vector at boundary: nvec = tvec1 x tvec2
!
!
SUBROUTINE compute_At_bcs(nshape,chi,dq,tvec1,tvec2,nvec,At1,At2)

  IMPLICIT NONE 

  ! INPUT
  INTEGER(IT),DIMENSION(2)                  ,INTENT(IN) :: nshape
  REAL(FP)   ,DIMENSION(nshape(1),nshape(2)),INTENT(IN) :: chi
  REAL(FP)                                  ,INTENT(IN) :: dq
  REAL(FP)   ,DIMENSION(3)                  ,INTENT(IN) :: tvec1,tvec2,nvec
  
  ! OUTPUT
  REAL(FP),DIMENSION(nshape(1),nshape(2)),INTENT(OUT) :: At1,At2
  
  ! LOCAL
  REAL(FP)              :: fac,dchi_dq1,dchi_dq2
  REAL(FP),DIMENSION(3) :: grad_cross_n
  INTEGER(IT)           :: i,j,k,nq1,nq2

  ! Compute spacing factor 
  fac = REAL(1,DP)/(REAL(2,DP)*dq)
  
  ! Get size
  nq1 = SIZE(At1,1)
  nq2 = SIZE(At1,2)
  
  !$OMP PARALLEL DO PRIVATE(i,j,k,dchi_dq1,dchi_dq2,grad_cross_n)
  DO j=1,nq2
    DO i=1,nq1
    
      ! Regular gradient 
      IF((i==1).OR. (i==nq1)) THEN
        dchi_dq1 = 0
      ELSE
        dchi_dq1 = fac*(chi(i+1,j)-chi(i-1,j))
      ENDIF
      
      IF((j==1).OR. (j==nq2)) THEN
        dchi_dq2 = 0
      ELSE
        dchi_dq2 = fac*(chi(i,j+1)-chi(i,j-1))        
      ENDIF
    
      ! Compute grad(chi) x n
      grad_cross_n = 0
      grad_cross_n = dchi_dq1*CROSS(tvec1,nvec) + dchi_dq2*CROSS(tvec2,nvec)
       
      ! Project onto unit vectors
      At1(i,j) = -DOT_PRODUCT(tvec1,grad_cross_n)
      At2(i,j) = -DOT_PRODUCT(tvec2,grad_cross_n)       
             
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 

END SUBROUTINE

! --------------------------------------------------------------------  
!
!>@name cross
!
!>@brief Compute the cross product of vectors A and B, i.e AxB
!
!>@details
!! Computes
!!\f[
!!  \mathbf C = (a_2 b_3 - a_3 b_2)\hat{x}
!!            + (a_3 b_1 - a_1 b_3)\hat{y}
!!            + (a_1 b_2 - a_2 b_1)\hat{z}.
!!\f]
!
!>@param[in] A Cartesian vector
!>@param[in] B Cartesian vector 
!>@return    Cross product of A and B
!
FUNCTION cross(A,B) RESULT(C)

  IMPLICIT NONE

  ! INPUT
  REAL(FP),DIMENSION(3),INTENT(IN) :: A,B
  
  ! RETURN VALUE
  REAL(FP),DIMENSION(3) :: C 
  
  ! Components of cross product
  C(1) = A(2)*B(3) - A(3)*B(2)
  C(2) = A(3)*B(1) - A(1)*B(3)
  C(3) = A(1)*B(2) - A(2)*B(1)
    
END FUNCTION

! ---------------------------------------------------------------------

FUNCTION trapz_2D(n1,n2,dq1,dq2,f)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)              ,INTENT(IN) :: n1,n2
  REAL(FP)                 ,INTENT(IN) :: dq1,dq2
  REAL(FP),DIMENSION(n1,n2),INTENT(IN) :: f
  
  ! RETURN
  REAL(FP) :: trapz_2D
  
  ! LOCAL
  REAL(FP),DIMENSION(:,:),ALLOCATABLE :: w
  REAL(FP),PARAMETER                  :: inv4 = REAL(1,FP)/REAL(4,FP)
  REAL(FP),PARAMETER                  :: inv2 = REAL(1,FP)/REAL(2,FP)
  
  ALLOCATE(w(n1,n2))
   
  ! Volume
  w = 1
  
  ! Edge: Must be set before corners 
  w(1,: ) = inv2
  w(:,1 ) = inv2
  w(n1,:) = inv2
  w(:,n2) = inv2
  
  ! Corners
  w(1,1 )  = inv4
  w(1,n2)  = inv4
  w(n1,1)  = inv4
  w(n1,n2) = inv4
  
  trapz_2D = SUM(w*f)*dq1*dq2

END FUNCTION
  
END MODULE
