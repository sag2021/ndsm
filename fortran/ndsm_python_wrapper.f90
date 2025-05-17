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

MODULE NDSM_PYTHON_WRAPPER

  USE,INTRINSIC :: ISO_C_BINDING,ONLY:C_INT,C_SIZE_T,C_DOUBLE,C_PTR
  USE NDSM_ROOT
  USE NDSM_VECTOR_POTENTIAL

  IMPLICIT NONE

  CONTAINS

! ---------------------------------------------------------------------
!
!>@name ndsm_vector_solve
!
!>@brief Compute vector potential A and magnetic field B = curl(A)
!>@details
!!  Built to be called from Python via ctypes
!!
!
!>@param[in]     nsize:   Total number of points in B: nx*ny*nz*3 
!>@param[in]     nshape4: Shape of B: (nx,ny,nx,3)
!>@param[inout]  ioptc:   Integer-value options 
!>@parama[inout] ropt:    Real-valued options
!>@param[in]     x: Mesh vector x, size nx 
!>@param[in]     y: Mesh vector y, size ny 
!>@param[in]     z: Mesh vector z, size nz  
!>@param[out]    A: Vector potential, (nx,ny,nz,3)
!>@param[inout]  B: Magnetic field, (nx,ny,nz,3)
!
!
FUNCTION ndsm_vector_solve(nsize,nshape4,ioptc,ropt,x,y,z,A,B) BIND(C) RESULT(ierr)

  IMPLICIT NONE

  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "ndsm_vector_solve"

  ! INPUT
  INTEGER(C_SIZE_T),VALUE                         :: nsize
  INTEGER(C_INT),DIMENSION(4)         ,INTENT(IN) :: nshape4
  REAL(C_DOUBLE),DIMENSION(nshape4(1)),INTENT(IN) :: x
  REAL(C_DOUBLE),DIMENSION(nshape4(2)),INTENT(IN) :: y
  REAL(C_DOUBLE),DIMENSION(nshape4(3)),INTENT(IN) :: z

  ! INPUT/OUTPUT
  INTEGER(C_INT),DIMENSION(0:IOPT_LEN-1),INTENT(INOUT) :: ioptc
  REAL(C_DOUBLE),DIMENSION(0:IOPT_LEN-1),INTENT(INOUT) :: ropt
  REAL(C_DOUBLE),DIMENSION(nsize)       ,INTENT(INOUT) :: b
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: A
  
  ! RETURN
  INTEGER(C_INT) :: ierr
  
  ! LOCAL 
  REAL(FP)                  :: tend,tstart
  TYPE(MG_PTR),DIMENSION(3) :: mesh
  !REAL(FP)                  :: dq
  INTEGER(IT)               :: iopt(0:IOPT_LEN-1)
  
  ! SHAPE
  INTEGER(IT),DIMENSION(4) :: nshape

  ! LOOP
  INTEGER :: i,j
  
  ! Copy in 
  nshape = nshape4
  iopt   = ioptc

  ! Set debug flag 
  DEBUG = (iopt(IOPT_DEBUG) == 1)

  ! Start
  tstart = get_cpu_time()
        
  ! ====================
  ! CONSTRUCT 3D MESH
  ! ====================

  ! Mesh spacing: Uniform in all dimensions
  !
  !dq = REAL(1,FP)/(nshape(1)-REAL(1,FP))
    
  ! Allocate memory
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Allocating memory for mesh...")
  DO i=1,SIZE(mesh)
    ALLOCATE(mesh(i)%val(nshape(i)))
  ENDDO
  
  ! Copy mesh
  mesh(1)%val = x
  mesh(2)%val = y
  mesh(3)%val = z
  
  ! =============
  ! SOLVE 
  ! =============
  !
  ! Call NDSM vector potential subroutine
  !
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Calling compute_vector_potential...")
  CALL compute_vector_potential(nshape,iopt,ropt,mesh,A,B)
    
  ! Check for generic error,
  !ierr = iopt(IOPT_IERR)
  
  ! =============
  ! MEMORY FREE
  ! ============
    
  ! Free mesh
  IF(DEBUG) CALL debug_msg(THIS_SUB,"Deallocating memory for mesh...")
  DO i=1,SIZE(mesh)
    DEALLOCATE(mesh(i)%val)
  ENDDO

  ! End time
  tend = get_cpu_time()

  ! Save time
  ROPT(ROPT_TIM) = tend-tstart  

  ! Copy out
  ioptc = iopt

  IF(DEBUG) CALL debug_msg(THIS_SUB,"Exiting Fortran lib...")
    
END FUNCTION

! ---------------------------------------------------------------------

!>@name get_iopt_len
!>@brief Returns IOPT_LEN parameter
FUNCTION get_iopt_len() BIND(C)
  IMPLICIT NONE
  INTEGER(C_INT) :: get_iopt_len
  get_iopt_len   = IOPT_LEN
END FUNCTION

FUNCTION get_iopt_ierr() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val   = IOPT_LEN
END FUNCTION

FUNCTION get_iopt_ms() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = IOPT_MS
END FUNCTION

FUNCTION get_iopt_ncycles() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = IOPT_NCYCLES
END FUNCTION

FUNCTION get_iopt_debug() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = IOPT_DEBUG
END FUNCTION

FUNCTION get_ropt_tim() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = ROPT_TIM
END FUNCTION

FUNCTION get_ropt_vtol() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = ROPT_VTOL
END FUNCTION

FUNCTION get_ropt_ctol() BIND(C) RESULT(val)
  IMPLICIT NONE
  INTEGER(C_INT) :: val
  val = ROPT_CTOL
END FUNCTION

END MODULE
