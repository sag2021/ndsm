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

! --------------------------------------------------------------------
!
! Test that the restriction/interpolation operators 
! have the Galerkin property. 
!
! Restriction/interpolation operators R[]/P[] satisfy the 
! Galerkin property if for any function, u, 
!
! < u_c | R[u_f] > = < P[u_c] | u_f > .
!  
! Where u_c is function defined on the coarse mesh, u_f is the
! function defined on the fine grid, and < | > are inner products
! defined on the coarse and fine grids. The inner product is defined
! as 
!           _
! <u,v> =  \
!          /_ u(i)*v(i) * dV.
!
! This program computes an analytic function on a coarse and fine
! mesh. The function is then interpolated/restricted to the 
! fine/coarse mesh respectively and the LHS and RHS inner products
! are computed.
!
!
! USAGE NOTES:
!
! The size of the mesh is in terms of length (not number of points)
! is random. The max value is set by LFACTOR. 
! 
!

PROGRAM MAIN

  USE NDSM_INTERP
  
  IMPLICIT NONE

  ! USER DEFINED  
  INTEGER(IT),PARAMETER  :: NDIM    = 4   ! Dimensions
  REAL(FP)   ,PARAMETER :: LFACTOR = 16  ! Max extent of mesh
  INTEGER(IT),PARAMETER  :: LBC     = 5   ! Array bound
  
  ! MAIN ARRAY  
  REAL(FP)   ,DIMENSION(:),ALLOCATABLE :: u_c,u_rc,u_f,u_if
  INTEGER(IT),DIMENSION(NDIM)          :: nvec_c,nvec_f
  INTEGER(IT)                          :: nsize_c,nsize_f
  
  ! MESH
  TYPE(MG_PTR),DIMENSION(NDIM)     :: mesh_c,mesh_f
  REAL(FP)    ,DIMENSION(NDIM,2)   :: ext        !
  REAL(FP)    ,DIMENSION(NDIM)     :: dq_c,dq_f         
  INTEGER(IT)                      :: nt_c,nt_f    
  
  ! COEFFICENT TERMS
  REAL(FP),DIMENSION(NDIM) :: B,M
  
  ! LOOP
  INTEGER(IT),DIMENSION(NDIM) :: ivec
  INTEGER(IT)                 :: n,i,lb
  
  ! MISC
  REAL(FP),DIMENSION(NDIM) :: q0
  REAL(FP)                 :: u_ana,u_num,tmp,mean
  
  ! ERRORS
  REAL(FP) :: inner_c,inner_f
    
  ! Construct extent
  !
  DO n=1,NDIM
    nvec_f(n) = 32
    nvec_c(n) = 15
  ENDDO

  ! Compute the total number of points
  nt_c = PRODUCT(nvec_c)
  nt_f = PRODUCT(nvec_f)

  ! Welcome message
  PRINT *,"UNIT TEST: Galerkin"  
  PRINT *,"Dimensions:                         ",NDIM
  PRINT *,"Total number of coarse mesh points: ",nt_c
  PRINT *,"Total number of fine mesh points:   ",nt_f
  
  !
  ! Construct mesh that spans random extent
  !
  DO n=1,NDIM
  
    ! Random extent
    CALL RANDOM_NUMBER(ext(n,:))
    ext(n,:) = ext(n,:)*LFACTOR
    
    ! Ensure that ext(i,2) is larger than ext(i,1)
    IF(ext(n,1) > ext(n,2)) ext(n,1:2) = ext(n,2:1:-1)
  
    ! Get spacing
    dq_c(n) = (ext(n,2) - ext(n,1))/REAL(nvec_c(n)-1,FP)
    dq_f(n) = (ext(n,2) - ext(n,1))/REAL(nvec_f(n)-1,FP)
  
    ! Generate mesh
    ALLOCATE(mesh_c(n)%val(1-LBC:nvec_c(n)-LBC))
    ALLOCATE(mesh_f(n)%val(1-LBC:nvec_f(n)-LBC))

    lb = LBOUND(mesh_c(n)%val,1)
    DO i=1,nvec_c(n)    
      mesh_c(n)%val(i+lb-1) = (i-1)*dq_c(n) + ext(n,1)
    ENDDO

    lb = LBOUND(mesh_f(n)%val,1)    
    DO i=1,nvec_f(n)
      mesh_f(n)%val(i+lb-1) = (i-1)*dq_f(n) + ext(n,1)
    ENDDO

  ENDDO

  ! Allocate array to hold data at mesh points
  ALLOCATE(u_c(nt_c),u_rc(nt_c))
  ALLOCATE(u_f(nt_f),u_if(nt_f))
  
  ! Populate mesh with original functions
  CALL populate_mesh(ndim,nvec_c,mesh_c,u_c)
  CALL populate_mesh(ndim,nvec_f,mesh_f,u_f)

  ! Create interpolated fine grid
  PRINT *,"Interpolating..."
  DO n=1,nt_f

    ivec = lin2nd(ndim,nvec_f,n)
    
    DO i=1,NDIM
      lb    = LBOUND(mesh_f(i)%val,1)
      q0(i) = mesh_f(i)%val(ivec(i)+lb-1)
    ENDDO
    u_if(n) = ninterp(ndim,nt_c,nvec_c,mesh_c,q0,u_c)
  
  ENDDO
  
  ! Create restricted coarse mesh
  PRINT *,"Restricting..."  
  DO n=1,nt_c
  
    ivec = lin2nd(ndim,nvec_c,n)
    
    DO i=1,NDIM
      lb    = LBOUND(mesh_c(i)%val,1)
      q0(i) = mesh_c(i)%val(ivec(i)+lb-1)
    ENDDO
    
    u_rc(n) = nrestrict(ndim,nt_c,nt_f,nvec_c,nvec_f,mesh_c,mesh_f,q0,u_f)
  
  ENDDO
 
  ! Compute inner products
  inner_c = inner_product(ndim,nvec_c,mesh_c,u_c,u_rc)
  inner_f = inner_product(ndim,nvec_f,mesh_f,u_f,u_if)
  
  ! Print results
  PRINT *,"Difference should be of order rounding error"
  PRINT *,"LHS inner product   < u_c    | R[u_f]> :",inner_c
  PRINT *,"RHS inner product   < P[u_c] | u_f   > :",inner_f
  PRINT *,"R[.]: Restriction operator"
  PRINT *,"P[.]: Interpolation operator"
  
  STOP "DONE"
    
 CONTAINS

! ---------------------------------------------------------------------
!
!>@name populate_mesh
!
!>@brief Compute the function u on the given mesh

SUBROUTINE populate_mesh(ndim,nvec,mesh,u)

  ! INPUT
  INTEGER(IT)                 ,INTENT(IN) :: ndim
  INTEGER(IT) ,DIMENSION(ndim),INTENT(IN) :: nvec
  TYPE(MG_PTR),DIMENSION(ndim),INTENT(IN) :: mesh
  
  ! OUTPUT
  REAL(FP),DIMENSION(*),INTENT(OUT) :: u
  
  ! LOCAL
  INTEGER(IT) :: n,i,nt,ivec(ndim),lb
  
  ! Total number of points
  nt = PRODUCT(nvec)
  
  DO n = 1,nt
    
    ! Construct NDIM index vector for linear index n
    ivec = lin2nd(ndim,nvec,n)
    
    DO i=1,NDIM
      lb    = LBOUND(mesh_f(i)%val,1)
      q0(i) = mesh(i)%val(ivec(i)+lb-1)
    ENDDO
    
    ! Construct point at q0
    u(n) = test_function(ndim,q0)
  
  ENDDO

END SUBROUTINE
  
! ---------------------------------------------------------------------

!>@name test_function
!
!>@brief A nonlinear test function
!
PURE FUNCTION test_function(ndim,q0) RESULT(f)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)             ,INTENT(IN) :: ndim
  REAL(FP),DIMENSION(ndim),INTENT(IN) :: q0
  
  ! RETURN 
  REAL(FP) :: f
  
  ! LOCAL
  INTEGER :: i
  
  f = 1
  DO i=1,ndim
    IF(MOD(i,2)== 0) THEN
      f = f*SIN(q0(i))
    ELSE
      f = f*COS(q0(i))
    ENDIF
  ENDDO
  
END FUNCTION
     
END PROGRAM
