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
!>@name ndsm_vector
!
!>@brief Computes vector potential and B = curl(A)
!
!>@details
!!
!! Wrapper function designed to be called like a C function from
!! IDL
!!
!>@param[in] argc         Number of arguments
!>@param[in] argv_cptr    Void C pointer to arguments 
!
!
FUNCTION ndsm_vector_solve(nsize,nshape4,ioptc,ropt,x,y,z,A,B) BIND(C) RESULT(ierr)

  IMPLICIT NONE

  ! INPUT
  INTEGER(C_SIZE_T),VALUE                         :: nsize
  INTEGER(C_INT),DIMENSION(4)         ,INTENT(IN) :: nshape4
  REAL(C_DOUBLE),DIMENSION(nshape4(1)),INTENT(IN) :: x
  REAL(C_DOUBLE),DIMENSION(nshape4(2)),INTENT(IN) :: y
  REAL(C_DOUBLE),DIMENSION(nshape4(3)),INTENT(IN) :: z

  ! INPUT/OUTPUT
  INTEGER(C_INT),DIMENSION(IOPT_LEN),INTENT(INOUT) :: ioptc
  REAL(C_DOUBLE),DIMENSION(IOPT_LEN),INTENT(INOUT) :: ropt
  REAL(C_DOUBLE),DIMENSION(nsize)   ,INTENT(INOUT) :: b
  
  ! OUTPUT
  REAL(FP),DIMENSION(nsize),INTENT(OUT) :: A
  
  ! RETURN
  INTEGER(C_INT) :: ierr
  
  ! LOCAL 
  REAL(FP)                  :: tend,tstart
  TYPE(MG_PTR),DIMENSION(3) :: mesh
  !REAL(FP)                  :: dq
  INTEGER(IT)               :: iopt(IOPT_LEN)
  
  ! SHAPE
  INTEGER(IT),DIMENSION(4) :: nshape

  ! LOOP
  INTEGER :: i,j
  
  ! Copy in 
  nshape = nshape4
  iopt   = ioptc

  ! Start
  tstart = get_cpu_time()
        
  ! ====================
  ! CONSTRUCT 3D MESH
  ! ====================

  ! Mesh spacing: Uniform in all dimensions
  !
  !dq = REAL(1,FP)/(nshape(1)-REAL(1,FP))
    
  ! Allocate memory
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
  CALL compute_vector_potential(nshape,iopt,ropt,mesh,A,B)
    
  ! Check for generic error
  ierr = iopt(IOPT_IERR)
  
  ! =============
  ! MEMORY FREE
  ! ============
    
  ! Free mesh
  DO i=1,SIZE(mesh)
    DEALLOCATE(mesh(i)%val)
  ENDDO

  ! End time
  tend = get_cpu_time()

  ! Save time
  ROPT(ROPT_TIM) = tend-tstart  

  ! Copy out
  ioptc = iopt
    
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
