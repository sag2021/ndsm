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
!>@name ndsm_root
!
!>@brief Root NDSM code
!
!>@details
!!
!! The root module has these functions:
!!
!! 1) Define real/integer types
!!
!! 2) Define error message subroutines
!!
!! 3) Define array index manipulation subroutines, i.e. 
!!    the subroutines that map 1D to ND indices etc.
!!
!! 4) Define time keeping subroutine
!!

MODULE NDSM_ROOT
 
  USE, INTRINSIC :: ISO_C_BINDING,ONLY:C_DOUBLE,C_INT64_T

  IMPLICIT NONE
  
  ! =================
  ! DEFINE REAL TYPES
  ! =================
  !
  INTEGER,PARAMETER :: DP = KIND(1.D0)  !< Double
  INTEGER,PARAMETER :: SP = KIND(1.)    !< Single
  INTEGER,PARAMETER :: FP = DP          !< Default float precision

  ! Default INT
  INTEGER,PARAMETER :: IT = C_INT64_T

  ! For debug
  LOGICAL :: debug

  ! ====================
  ! DEFINE: PTR_TO_ARRAY
  ! ====================
  !
  ! Derived type used to create arrays of arrays
  !
  TYPE :: MG_PTR
    REAL(FP),DIMENSION(:),ALLOCATABLE :: val
  END TYPE

  !>@name assert_n_inbounds
  !>@brief Generic interface for array bounds assertions
  INTERFACE assert_n_inbounds
    MODULE PROCEDURE  assert_n_inbounds_ND
    MODULE PROCEDURE  assert_n_inbounds_1D
  END INTERFACE
    
 CONTAINS

! ---------------------------------------------------------------------
!
!>@name nd2lin
!
!>@brief Convert N-dimensional indices to linear index n
!
!>@details
!! An N-D array can be indexed either by a vector of length N
!! or a linear index n.
!!
!! In accordance with ISO/IEC 1539:1991, an N-D Fortran 
!! array maps to a 1D array with "linear" index:
!!
!! A(s1,s2,s3,...) = A(n)
!!
!! Where, from Table 6.1 of the ISO standard,
!!
!! n = 1 + (s1-j1) + (s2-j2)*d1 + (s3-j3)*d1*d2 + ...
!!
!! Here j1,j2,j3, ... are the starting indices in each dimension
!! (by default 1, unless the array is declared with non-standard 
!! bounds). Here also, d1 = k1 - j1 + 1 ,d2 = k2 -j2 + 1, 
!! where k1,k2,k3 ... is the max index in dimension s1,s2,s3 etc.
!! 
!! In the following function, it is assumed that the data has 
!! the standard bounds, i.e. j1 = 1, j2 = 1,j3 = 1 ... etc. 
!! and k1 = d1, k2 = d2, ... etc.
!!
!!

FUNCTION nd2lin(ndim,nshape,ivec,tag)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: ivec
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nshape
  
  ! OPTIONAL
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: tag
  
  ! RETURN
  INTEGER(IT) :: nd2lin
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "nd2lin"
  
  ! LOCAL
  INTEGER(IT) :: i
  INTEGER(IT) :: dd,nsize
  CHARACTER(LEN=16) :: tag_
  
  ! Set 
  IF(PRESENT(tag)) THEN
    tag_ = tag
  ELSE
    tag_ = "AC5H"
  ENDIF
  
  dd     = 1
  nd2lin = ivec(1)
  DO i=1,ndim-1  
    dd     = dd*nshape(i)   
    nd2lin = nd2lin + (ivec(i+1)-1)*dd
  ENDDO
  
  ! ASSERT: Assert that nd2lin is in the range [1,NSIZE]
  nsize = PRODUCT(nshape)
  CALL assert_n_inbounds(ndim,nsize,nd2lin,ivec,THIS_SUB,TAG)
  
END FUNCTION

! ---------------------------------------------------------------------
!
!>@name shape_to_strides
!
!>@brief Convert shape to strides

FUNCTION shape_to_strides(ndim,nshape) RESULT(nstrides)

  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nshape

  ! RETURN
  INTEGER(IT),DIMENSION(ndim) :: nstrides

  ! LOCAL
  INTEGER(IT) :: i

  nstrides(1) = 1
  DO i=2,ndim
    nstrides(i) = nshape(i-1)*nstrides(i-1)
  ENDDO

END FUNCTION

! ---------------------------------------------------------------------
!
!>@name lin2ND
!
!>@brief Convert 1D index to N-D index vector
!
!>@details
!!
!! This is the inverse of ND2lin and has the same 
!! restrictions
!!

FUNCTION lin2nd(ndim,nshape,n)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,n
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nshape
  
  ! RETURN 
  INTEGER(IT),DIMENSION(ndim) :: lin2ND
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "lin2nd"
  
  ! LOCAL
  INTEGER(IT) :: i,nloc,nsize
  
  ! Copy 
  nloc = n-1
  
  ! Product: K = d_1 * d_2 * ... * d_(NDIM-1)
  nsize = PRODUCT(nshape(1:ndim-1))
    
  DO i=ndim,2,-1
    
    ! ASSERT: nsize is nonzero
    CALL assert_nonzero(nsize,THIS_SUB,"AC31")
        
    ! Integer divide 
    !
    lin2nd(i) = nloc/nsize 
    
    ! Subtract term
    !
    nloc  = nloc - nsize*(lin2nd(i))
    
    ! d_1*d_2*d_3*...*d_k -> d_1*d_2*d_3*...*d_(k-1)
    !
    nsize = nsize/nshape(i-1)
    
  ENDDO
  lin2nd(1) = nloc
  
  ! Default Fortran arrays start at 1
  lin2nd = lin2nd + 1
    
  ! ASSERT: lin2nd is in bounds
  CALL assert_n_inbounds_ND(ndim,nshape,lin2nd,THIS_SUB,"AC6")
   
END FUNCTION

! ---------------------------------------------------------------------
!
!>@name lin2nd_s
!
!>@brief Same as lin2nd, but excepts strides

FUNCTION lin2nd_s(ndim,nshape,nstrides,n)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,n
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nstrides,nshape
  
  ! RETURN 
  INTEGER(IT),DIMENSION(ndim) :: lin2ND_s
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "lin2nd_s"
  
  ! LOCAL
  INTEGER(IT) :: i,nloc
  
  ! Copy 
  nloc = n-1
  
  DO i=ndim,2,-1
    
    ! ASSERT: nsize is nonzero
    CALL assert_nonzero(nstrides(i),THIS_SUB,"AC319")
        
    ! Integer divide 
    !
    lin2ND_s(i) = nloc/nstrides(i) 
    
    ! Subtract term
    !
    nloc  = nloc - nstrides(i)*lin2nd_s(i)
        
  ENDDO
  lin2ND_s(1) = nloc
    
  ! ASSERT: lin2nd is in bounds
  CALL assert_n_inbounds_ND(ndim,nshape,lin2ND_s,THIS_SUB,"AC57")
   
END FUNCTION

! ---------------------------------------------------------------------
!
! ASSERTIONS
!
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
!
!>@name assert_nonzero
!
!>@brief Assert that n is not zero
!
!>@param[in] n       Integer
!>@param[in] subname Name of subroutine that calls the assert. Used in error message
!>@param[in] eid     ID string to identify the location of the call in source code

SUBROUTINE assert_nonzero(n,subname,eid)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)     ,INTENT(IN) :: n
  CHARACTER(LEN=*),INTENT(IN) :: subname,eid
  
  IF(n == 0) THEN
     CALL error_msg("Value must be non-zero",subname,eid)
     STOP "FATAL ERROR"  
  ENDIF
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name assert_n_inbounds_1D
!
!>@brief Assert that n is in range [1,NSIZE]
!
!>@details
!! Asserts that the integer n is in the range 1 <= n <= NSIZE.
!! 
!! A fatal error is raised if this condition is not met.
!! A fatal error causes a hard crash with a call to STOP.
!!
!
!>@param[in] ndim    Number of dimension 
!>@param[in] nsize   Integer: Array size
!>@param[in] n       Integer: Array index
!>@param[in] subname Name of subroutine that calls the assert. Used in error message
!>@param[in] eid     ID string to identify the location of the call in source code

SUBROUTINE assert_n_inbounds_1D(ndim,nsize,n,ivec,subname,eid)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,nsize,n
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: ivec
  CHARACTER(LEN=*)           ,INTENT(IN) :: subname,eid
  
  ! LOGICAL
  LOGICAL :: bounds_error
  
  ! Bound error check
  bounds_error = (n < 0) .OR. (n > nsize)

  ! Error in self-consistency
  IF(bounds_error) THEN
     CALL error_msg("Array bounds error 1D",subname,eid)
     PRINT *,"N:    ",N
     PRINT *,"NSIZE:",NSIZE
     PRINT *,"NDIM: ",NDIM
     !PRINT *,"IVEC: ",IVEC
     STOP "FATAL ERROR"
  ENDIF
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name assert_n_inbounds_ND
!
!>@brief Assert that a vector of indices is consistent with shape
!
!>@details
!! Asserts that  ALL(ivec >= 1) .AND. ALL(ivec <= nshape)
!! 
!! A fatal error is raised if this condition is not met.
!! A fatal error causes a hard crash with a call to STOP.
!!
!
!>@param[in] ndim    Number of dimension 
!>@param[in] nshape  Vector of integers
!>@param[in] ivec    Vector of integers
!>@param[in] subname Name of subroutine that calls the assert. Used in error message
!>@param[in] eid     ID string to identify the location of the call in source code

SUBROUTINE assert_n_inbounds_ND(ndim,nshape,ivec,subname,eid)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: ivec,nshape
  CHARACTER(LEN=*)           ,INTENT(IN) :: subname,eid
  
  ! LOGICAL
  LOGICAL :: bounds_error
  
  ! Bound error check
  bounds_error = ANY(ivec < 1) .OR. ANY(ivec > nshape)

  ! Error in self-consistency
  IF(bounds_error) THEN
     CALL error_msg("Array bounds error ND",subname,eid)
     STOP "FATAL ERROR"
  ENDIF
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name assert_size_shape_consistent
!
!>@brief Assert that PRODUCT(nshape) = nsize. 
!
!>@detail
!! Checks that PRODUCT(nshape) = nsize
!!
!! A fatal error is raised if this condition is not met.
!! A fatal error causes a hard crash with a call to STOP.
!!
!!
!
!>@param[in] ndim    Number of dimension 
!>@param[in] nsize   Integer 
!>@param[in] nshape  Vector of integers of length NDIM
!>@param[in] subname Name of subroutine that calls the assert. Used in error message.
!>@param[in] eid      ID string to identify the location of the call in source code
!
SUBROUTINE assert_size_shape_consistent(ndim,nsize,nshape,subname,eid)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                ,INTENT(IN) :: ndim,nsize
  INTEGER(IT),DIMENSION(ndim),INTENT(IN) :: nshape
  CHARACTER(LEN=*)           ,INTENT(IN) :: subname,eid
  
  ! LOGICAL
  LOGICAL :: self_consistency 
  INTEGER(IT) :: nt
  
  ! Self-consistency check
  nt                = PRODUCT(nshape)
  self_consistency  = (nt .ne. nsize)
  
  ! Error in self-consistency
  IF(self_consistency) THEN
     CALL error_msg("Shape and total number of points not consistent",subname,eid)
     STOP "FATAL ERROR"
  ENDIF
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name error_msg
!
!>@brief Print error message in a standard format
!
!>@details
!!  
!!  Print an error message in a standard format
!!
!!  ERROR(#sub):#msg:#eid
!!
!!  #sub: Name of subroutine that raises the error
!!  #msg: Details of error
!!  #eid: ID string to identify the location of the call in source code
!!
!! The message is written to ERROR_UNIT, defined in 
!! ISO_FORTRAN_ENV.
!!
SUBROUTINE error_msg(msg,sub,eid)

  USE, INTRINSIC :: ISO_FORTRAN_ENV,ONLY: ERROR_UNIT

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: msg,sub,eid

  WRITE(ERROR_UNIT,'(A)') "ERROR("//TRIM(sub)//"):"//TRIM(msg)//":"//TRIM(eid)

END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name debug_msg
!>@brief Standard format for debug message

SUBROUTINE debug_msg(sub,msg)

  USE, INTRINSIC :: ISO_FORTRAN_ENV,ONLY: ERROR_UNIT

  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: sub,msg

  WRITE(ERROR_UNIT,'(A)') "DEBUG("//TRIM(sub)//"):"//TRIM(msg)

END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@get_cpu_time
!
!>@brief Get CPU time from some arbitrary start point in seconds
!
!>@details
!! Uses Fortran CPU_TIME subroutine by default. 
!!
!! However, if compiled with OpenMP, the OMP wall time subroutine 
!! is used. Lines with !$ are treated as comments when OpenMP is not
!! used.
!!
!
!>@param[return] Wall time in seconds
!
FUNCTION get_cpu_time() RESULT(wtime)

!$ USE OMP_LIB,ONLY: OMP_GET_WTIME

  IMPLICIT NONE
  
  ! RETURN
  REAL(DP) :: wtime
  
  ! Fortran CPU time
  CALL CPU_TIME(wtime)
  
  ! Replace with OpenMP 
!$  wtime = OMP_GET_WTIME()
  
END FUNCTION

END MODULE
