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
!>@name ndsm_interp
!
!>@brief Performs interpolation/restriction in N dimensions
!
!>@details
!!
!! RESTRICTIONS:
!! -------------
!!  Mesh must be uniform in a given dimension, although it can 
!!  vary between dimensions
!!
!!  Note: The functions for converting between N-D and 1D indices
!!        assume standard Fortran bounds. However, the nrestrict/ninterp
!!        routines perform internal remappings, so the arrays passed
!!        to these routines may have non-default bounds (including
!!        the mesh vectors).
!! 
!!  IMPORTANT: Non-default bounds for the mesh vectors has not 
!!             been tested. [HAS BEEN TESTED NOW 30 JUNE 2019]
!!
!!
!! Mesh derived types:
!! -------------------
!!
!! A derived type (MG_PTR) is used to define an "vector of vectors".
!! The derived type has one field corresponding to a 1D allocatable
!! array. This makes it possible to store vectors of different 
!! sizes in one object.
!!
!!
!! Generic Name conventions:
!! -------------------------
!!
!! qvec: Vector of derived types that holds the mesh vector for each
!!       dimension, i.e. the mesh coordinates in dimension k are 
!!       qvec(k)%val(:). 
!!

MODULE NDSM_INTERP

  USE NDSM_ROOT

  IMPLICIT NONE
    
 CONTAINS

! ---------------------------------------------------------------------
!
!>@name ninterp
!
!>@brief Peform N-linear interpolation in N dimensions
!
!>@param[in] ndim    Number of dimensions
!>@param[in] nsize   Total number of points 
!>@param[in] nshape  Vector containing array shape information
!>@paran[in] qvecs   Mesh vectors
!>@param[in] q0      Point at which interpolation is performed
 
FUNCTION ninterp(ndim,nsize,nshape,qvecs,q0,f)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                        ,INTENT(IN) :: ndim,nsize
  INTEGER(IT)     ,DIMENSION(ndim)   ,INTENT(IN) :: nshape
  TYPE(MG_PTR),DIMENSION(ndim),TARGET,INTENT(IN) :: qvecs
  REAL(FP)    ,DIMENSION(ndim)       ,INTENT(IN) :: q0
  REAL(FP)    ,DIMENSION(nsize)      ,INTENT(IN) :: f
  
  ! RETURN
  REAL(FP) :: ninterp
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "ninterp"
  
  ! LOCAL
  REAL(FP)    ,DIMENSION(2**ndim) :: fs      !< Number of samples
  INTEGER(IT) ,DIMENSION(ndim,2 ) :: bpts    !< Bracket points
  REAL(FP)                        :: wh,wl,ql,qh,dq
  INTEGER(IT)                     :: i,j,nq,NC,ierr,lb

  ! ASSERT: Check that nsize is consistent with total size implied by nshape
  !
  CALL assert_size_shape_consistent(ndim,nsize,nshape,THIS_SUB,"AC1")
  
  ! Compute bracket points
  !
  ! Note: bracket points are returned in the default range [1,NQ].
  !       The subroutine get_interpolation_values expects bpts in
  !       this range. The vector mesh(i)%val may have non-standard
  !       bounds and must be addressed with bpts + LBOUND(mesh(i)%val,1)-1
  !       instead of bpts directly
  !
  DO i=1,ndim  
    nq               = SIZE(qvecs(i)%val)
    CALL find_bracket_points_uniform(qvecs(i)%val,nq,q0(i),bpts(i,1),bpts(i,2),ierr)  
  ENDDO
  
  ! Compute samples
  CALL get_interpolation_values(ndim,nsize,nshape,bpts,f,fs)
  
  DO i=ndim,1,-1
  
    ! Get lower bound
    lb = LBOUND(qvecs(i)%val,1)
  
    ! Coordinates of bracket points
    ql = qvecs(i)%val(bpts(i,1)+lb-1)
    qh = qvecs(i)%val(bpts(i,2)+lb-1)
  
    ! Distance between bracket points
    dq = qh - ql
    
    ! Interpolation weights
    wl = +(q0(i) - ql)/dq
    wh = -(q0(i) - qh)/dq

    ! Perform interpolation in dimension i. This involves
    ! a weighted sum over 2^(i-1) points. The data are packed into 
    ! the fs vector such that fs(j) and fs(j+NC) are at ql and qh
    ! respectively
    !
    NC = 2**(i-1)    
    DO j=1,NC
      fs(j) = wh*fs(j) + wl*fs(j+NC)          
    ENDDO
  
  ENDDO
  
  ninterp = fs(1)
  
END FUNCTION

! ---------------------------------------------------------------------
!
!>@name nrestrict
!
!>@brief Perform full-weighted restriction in N dimensions
!
!>@details
!! This restriction operator is the adjoin of the N-dimensional
!! interpolation operator defined by ninterp
!!
!! The mesh must have constant spacing in each dimension 
!! (although the spacing can differ between dimensions)
!!
!
!>@param[in] ndim      Number of dimensions
!>@param[in] nsize_c   Total number points in coarse mesh
!>@param[in] nsize_f   Total number of points in fine mesh
!>@param[in] nshape_c  Number of pts. in each dimension (coarse mesh)
!>@param[in] nshape_f  Number of pts. in each dimension (fine   mesh)
!>@param[in] qvec_c    Mesh vectors for coarse mesh
!>@param[in] qvec_f    Mesh vectors for fine mesh
!>@param[in] q0        Point where restriction is performed
!>@param[in] f         Function to be restricted. Defined on fine mesh. Must have size PRODUCT(nvec_f)
!
!>@return              Function, f, restricted to coarse mesh point q0
!
FUNCTION nrestrict(ndim,nsize_c,nsize_f,nshape_c,nshape_f,qvec_c,qvec_f,q0,f) RESULT(fc)

  ! INPUT
  INTEGER(IT)                        ,INTENT(IN) :: ndim,nsize_c,nsize_f
  INTEGER(IT)     ,DIMENSION(ndim)   ,INTENT(IN) :: nshape_c,nshape_f
  TYPE(MG_PTR),DIMENSION(ndim),TARGET,INTENT(IN) :: qvec_c,qvec_f
  REAL(FP)    ,DIMENSION(ndim)       ,INTENT(IN) :: q0
  REAL(FP)    ,DIMENSION(nsize_f)    ,INTENT(IN) :: f
  
  ! RETURN
  REAL(FP) :: fc
  
  ! THIS SUBROUTINE
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "nrestrict"
  
  ! LOCAL
  REAL(FP)                      :: w,c1,c2
  INTEGER(IT)                   :: qil,qih,i,j,n,nq,ierr,nsize_s,lb_f,lb_c
  INTEGER(IT),DIMENSION(ndim,2) :: bpts
  INTEGER(IT),DIMENSION(ndim)   :: nshape_s,ivec,qic
  REAL(FP)   ,DIMENSION(ndim)   :: qc,dq_c,dq_f,w2
    
  ! ASSERT: Assert that nsize_c/nsize_f is consistent with the array size 
  !         implied by the shape vector nshape_c/nshape_f
  !
  CALL assert_size_shape_consistent(ndim,nsize_c,nshape_c,THIS_SUB,'AC2')
  CALL assert_size_shape_consistent(ndim,nsize_f,nshape_f,THIS_SUB,'AC3')
    
  !
  ! Compute the range of indices on the fine mesh that contribute
  ! to the restricted point on the coarse mesh
  !
  DO i=1,ndim
  
    ! Get lower bounds
    lb_c = LBOUND(qvec_c(i)%val,1)
    lb_f = LBOUND(qvec_f(i)%val,1)
        
    ! Spacing
    dq_c(i) = qvec_c(i)%val(lb_c+1) - qvec_c(i)%val(lb_c)
    dq_f(i) = qvec_f(i)%val(lb_f+1) - qvec_f(i)%val(lb_f)

    ! Weighting
    w2(i) = dq_f(i)/dq_c(i)**2   
    
    ! Bracket points
    nq = SIZE(qvec_f(i)%val)
    qil = 0 ; qih = 0
    CALL find_bracket_points_uniform(qvec_f(i)%val,nq,q0(i)-dq_c(i),qil,qih,ierr) 
    IF(ierr < 0) THEN
      bpts(i,1) = qil  
    ELSE
      bpts(i,1) = qih    
    ENDIF
    
    qil = 0 ; qih = 0
    CALL find_bracket_points_uniform(qvec_f(i)%val,nq,q0(i)+dq_c(i),qil,qih,ierr) 
    IF(ierr > 0) THEN
      bpts(i,2) = qih  
    ELSE
      bpts(i,2) = qil    
    ENDIF
      
    ! Stencil shape
    nshape_s(i) = bpts(i,2)- bpts(i,1) + 1
        
  ENDDO

  ! Compute total size of stencil
  nsize_s = PRODUCT(nshape_s)
  
  ! Initialize variable to zero
  fc = 0
  
  !
  ! Perform restriction
  !
  DO j=1,nsize_s
    
    ! Convert local linear index "j" into N-Dimensional
    ivec = 0
    ivec = lin2nd(ndim,nshape_s,j)
        
    ! Get qc
    DO i=1,NDIM
      lb_f    = LBOUND(qvec_f(i)%val,1)
      qic(i)  = bpts(i,1) + ivec(i) - 1
      qc(i)   = qvec_f(i)%val(qic(i)+lb_f-1) 
    ENDDO
    
    ! Compute weight
    w = 1
    DO i=1,NDIM
      c1 = ABS(qc(i) - q0(i))
      c2 = ABS(dq_c(i)  - c1)
      w  = w * c2 *w2(i)
    ENDDO
    
    ! Compute 1D index
    n = 0
    n = nd2lin(ndim,nshape_f,qic)

    fc = fc + w*f(n)
    
  ENDDO

END FUNCTION

! ---------------------------------------------------------------------
!
!>@name get_interpolation_values
!
!>@brief Takes data points used for interpolation and packs them into a vector of length 2^N
!
!>@details
!! For a function f, defined on a mesh. The subroutine takes the 
!! values of f for a given cell defined by bpts and packs them
!! into a vector fs. 
!!
!!
!
!>@param[in] ndim Number of dimensions
!>@param[in] nvec Defines shape on f array
!>@param[in] bpts The coordinates of the cell boundary
!>@param[in] f    A function defined on a mesh in N dimensions
!
!>@param[out] fs  Values of f at the 2^N cell points

SUBROUTINE get_interpolation_values(ndim,nsize,nshape,bpts,f,fs)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                     ,INTENT(IN) :: ndim,nsize
  INTEGER(IT) ,DIMENSION(ndim)    ,INTENT(IN) :: nshape  
  INTEGER(IT) ,DIMENSION(ndim,0:1),INTENT(IN) :: bpts
  REAL(FP),DIMENSION(nsize)       ,INTENT(IN) :: f
  
  ! OUTPUT
  REAL(FP),DIMENSION(2**ndim),INTENT(OUT) :: fs
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "get_interpolation_values"
  
  ! LOCAL
  INTEGER(IT),DIMENSION(ndim)    :: svec
  INTEGER(IT),DIMENSION(ndim)    :: ivec
  INTEGER(IT)                    :: i,j,n
  INTEGER(IT),PARAMETER          :: i1 = 1
  
  ! ASSERT: PRODUCT(nshape)= nsize
  !
  CALL assert_size_shape_consistent(ndim,nsize,nshape,THIS_SUB,"AC4")
  
  DO n = 1,2**ndim
  
    ! Convert the binary representation of n-1 into an array
    ! of ones and zeros
    svec = 0
    CALL bit2array(n-i1,svec)
    
    ! N-Dimensional array index vector 
    DO i=1,NDIM
      ivec(i) = 0
      ivec(i) = bpts(i,svec(i)) 
    ENDDO
    
    ! Convert to linear index 
    j = 0
    j = nd2lin(ndim,nshape,ivec,"KCA")
 
    ! ASSERT: Assert that j is in range [1,NSIZE]
    CALL assert_n_inbounds(ndim,nsize,j,ivec,THIS_SUB,"ACB5")
     
    ! Add element n
    fs(n) = f(j)
  
  ENDDO
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name find_bracket_points_uniform
!
!>@brief Find mesh points that bracket point q0

SUBROUTINE find_bracket_points_uniform(qvec,nq,q0,qil,qih,ierr)

  IMPLICIT NONE
    
  ! INPUT 
  INTEGER(IT)           ,INTENT(IN) :: nq
  REAL(FP),DIMENSION(nq),INTENT(IN) :: qvec
  REAL(FP)              ,INTENT(IN) :: q0

  ! OUTPUT 
  INTEGER(IT),INTENT(OUT) :: qil,qih
  INTEGER(IT),INTENT(OUT) :: ierr
  
  ! NAME
  CHARACTER(LEN=*),PARAMETER :: THIS_SUB = "find_bracket_points_uniform"
  
  ! LOCAL
  REAL(FP) :: dq
  
  ! Special case of nq = 1
  IF(nq == 1) THEN
    CALL error_msg("Vector has length 1","THIS_SUB","ERC4")  
    STOP "FATAL ERROR"
  ENDIF
  
  ! Point below bottom of vector
  IF(q0 .le. qvec(1)) THEN
      qil  = 1 
      qih  = 2
      ierr = -1
    RETURN
  ENDIF

  ! Point above top of vector
  IF(q0 .ge. qvec(nq)) THEN
      qil  = nq - 1
      qih  = nq 
      ierr = +1
    RETURN
  ENDIF
  
  ! Determine grid spacing
  dq = qvec(2)-qvec(1)
    
  ! Compute lower index. Need to add 1 since Fortran vectors start
  ! at 1 by default
  qil = FLOOR((q0-qvec(1))/dq) + 1

  !
  ! Rounding error in the quotiont (q0-qvec(1))/dq can 
  ! lead to qil == nq without q0 < qvec(nq)
  !
  IF(qil .ge. nq) THEN
    qil = nq - 1
    qih = nq
  ELSE
    qih = qil + 1
  ENDIF
  
  ! No error
  ierr = 0
  
END SUBROUTINE                              

! ---------------------------------------------------------------------
!
!>@name bit2array
!
!>@brief Convert bit field to array of ones and zeros

SUBROUTINE bit2array(bit_field,bit_array)

  ! INPUT
  INTEGER(IT),INTENT(IN) :: bit_field

  ! OUTPUT
  INTEGER(IT),DIMENSION(:),INTENT(OUT) :: bit_array

  ! LOCAL
  INTEGER(IT) :: i,bit_field_copy
  INTEGER(IT),PARAMETER :: i1 = 1

  
  ! Check that integer has enough bits to match array.
  IF(BIT_SIZE(bit_field)-1 < SIZE(bit_array)) THEN
    CALL error_msg("Integer size is too small","bit2array","ERC$@")
    STOP "FATAL ERROR"
  ENDIF
  
  ! Copy bit field
  bit_field_copy = bit_field
  
  DO i=1,SIZE(bit_array)
  
    ! Set element i to the bit i
    bit_array(i) = IAND(bit_field_copy,i1)
      
    ! Shift bits right 
    bit_field_copy = ISHFT(bit_field_copy,-1)
  
  ENDDO
  
END SUBROUTINE

! ---------------------------------------------------------------------

PURE FUNCTION nlinear_function(ndim,q0,M,B) RESULT(f)

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)             ,INTENT(IN) :: ndim
  REAL(FP),DIMENSION(ndim),INTENT(IN) :: q0
  REAL(FP),DIMENSION(ndim),INTENT(IN) :: M,B
  
  ! RETURN VALUE
  REAL(FP) :: f
  
  ! LOCAL
  INTEGER :: i
    
  f = 1
  DO i=1,ndim
    f = (M(i)*q0(i) + B(i))*f 
  ENDDO
  
END FUNCTION

! ---------------------------------------------------------------------

FUNCTION inner_product(ndim,nvec,qvec,u,v) RESULT(outp)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                 ,INTENT(IN) :: ndim
  INTEGER(IT) ,DIMENSION(ndim),INTENT(IN) :: nvec
  TYPE(MG_PTR),DIMENSION(ndim),INTENT(IN) :: qvec
  REAL(FP),DIMENSION(*)       ,INTENT(IN) :: u,v
  
  ! RETURN VALUE
  REAL(FP) :: outp
  
  ! LOCAL
  REAL(FP),DIMENSION(ndim) :: dq
  REAL(FP)                 :: dV
  INTEGER(IT)              :: i,ntotal
  
  ! Initialize output
  outp = 0
  
  ! Spacing
  DO i=1,ndim
    dq(i) = qvec(i)%val(2) - qvec(i)%val(1)
  ENDDO
  dV = PRODUCT(dq)
  
  ! Total size
  ntotal = PRODUCT(nvec)
  
  ! Product
  !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:outp) ORDERED
  DO i =1,ntotal 
    outp = outp + u(i)*v(i)      
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Scale by cell size
  outp = outp*dV
  
END FUNCTION

END MODULE
