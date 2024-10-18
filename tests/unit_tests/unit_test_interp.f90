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
! Unit test for interpolation N-Dimensional interpolation
!
! Tests interpolation module by performing interpolation on
! the interpolation function that defines the method. The 
! particular function can always be reproduced without truncation
! error. The test is a success if the error is consistent with 
! rounding error. The number of dimensions is fixed, but most 
! other parameters are chosen using pseudo-random numbers.
! 
! This program constructs a random mesh in NDIM dimensions. 
! The mesh has constant spacing in any particular dimension, but
! may differ between dimensions. An N-linear function is constructed
! on the mesh. Random points selected within extent of the mesh 
! (but not necessarily mesh points) are selected
!

PROGRAM UT_INTERP

  USE NDSM_INTERP
  
  IMPLICIT NONE
  
  ! USER-DEFINED PARAMETERS
  INTEGER(IT),PARAMETER :: NDIM     = 5   !< Number of dimensions
  REAL(FP)   ,PARAMETER :: LFACTOR  = 10  !< Max. mesh extent 
  INTEGER(IT),PARAMETER :: NSAMPLES = 16  !< Number of test points
  
  ! MAIN ARRAY  
  REAL(FP)   ,DIMENSION(:),ALLOCATABLE :: u     !< Data
  INTEGER(IT),DIMENSION(NDIM)          :: nvec  !< Shape
  
  ! MESH
  TYPE(MG_PTR),DIMENSION(NDIM)     :: mesh       !< Mesh vectors
  REAL(FP)    ,DIMENSION(NDIM,2)   :: ext        !< Mesh extent in each dimension
  REAL(FP)    ,DIMENSION(NDIM)     :: dq         !< Mesh spacing in each dimension
  INTEGER(IT)                      :: ntotal     !< Total number of points
  
  ! COEFFICENT TERMS
  REAL(FP),DIMENSION(NDIM) :: B,M
  
  ! LOOP
  INTEGER(IT),DIMENSION(NDIM) :: ivec
  INTEGER(IT)                 :: n,i
  
  ! MISC
  REAL(FP),DIMENSION(NDIM) :: q0
  REAL(FP)                 :: u_ana,u_num,tmp,mean
  
  ! ERRORS
  REAL(FP) :: max_a_error,max_r_error,err
  
  ! Seed RNG with time
  CALL seed_time()
  
  ! Construct extent
  !
  DO n=1,NDIM
    CALL RANDOM_NUMBER(tmp)
    nvec(n) = CEILING(31*tmp) + 1
  ENDDO

  ! Compute the total number of points
  ntotal = PRODUCT(nvec)

  ! Welcome message
  PRINT *,"UNIT TEST: Interpolation"
  
  PRINT *,"Dimensions:                  ",NDIM
  PRINT *,"Total number of mesh points: ",NTOTAL
  

  !
  ! Construct mesh
  !
  DO n=1,NDIM
  
    ! Random extent
    CALL RANDOM_NUMBER(ext(n,:))
    ext(n,:) = ext(n,:)*LFACTOR
    
    ! Ensure that ext(i,2) is larger than ext(i,1)
    IF(ext(n,1) > ext(n,2)) ext(n,1:2) = ext(n,2:1:-1)
  
    ! Get spacing
    dq(n) = (ext(n,2) - ext(n,1))/REAL(nvec(n)-1,FP)
  
    ! Generate mesh
    ALLOCATE(mesh(n)%val(nvec(n)))

    DO i=1,nvec(n)
      mesh(n)%val(i) = (i-1)*dq(n) + ext(n,1)
    ENDDO
      
  ENDDO
  
  ! Allocate array to hold data at mesh points
  ALLOCATE(u(ntotal))
  
  CALL RANDOM_NUMBER(M)
  CALL RANDOM_NUMBER(B)
  
  !
  ! Construct function on mesh
  !
  DO n = 1,ntotal
  
    ! Construct NDIM index vector for linear index n
    ivec = lin2nd(ndim,nvec,n)
    
    DO i=1,NDIM
      q0(i) = mesh(i)%val(ivec(i))
    ENDDO
    
    ! Construct point at q0
    u(n) = nlinear_function(ndim,q0,M,B)
  
  ENDDO

  ! Initialize 
  max_a_error = 0
  max_r_error = 0
  
  !
  ! Perform interpolation at random points and measure
  ! the error
  !
  DO n=1,NSAMPLES
  
    ! Construct random point inside domain
    CALL RANDOM_NUMBER(q0)
    DO i=1,NDIM
      q0(i) = q0(i)*(ext(i,2)-ext(i,1)) + ext(i,1)
    ENDDO
  
    ! Interpolate
    u_num = ninterp(ndim,ntotal,nvec,mesh,q0,u)
    
    ! Exact value
    u_ana = nlinear_function(ndim,q0,M,B)
    
    ! Compute error
    err = ABS(u_ana - u_num)
    
    ! Absolute Mean value
    mean = ABS(u_ana + u_num)/REAL(2,FP)
    
    ! Set errors
    max_a_error = MAX(max_a_error,err)
    max_r_error = MAX(max_r_error,err/mean)
        
  ENDDO
    
  PRINT *,"Max abs. error: ",max_a_error  
  PRINT *,"Max rel. error: ",max_r_error  

  STOP "DONE"    

CONTAINS

! ---------------------------------------------------------------------
!
!>@name seed_time
!>@brief Seed Fortran RNG with time/date. 
!>@details
!! Used to produce unpredictable random variates

SUBROUTINE seed_time()

  IMPLICIT NONE
  
  INTEGER :: seed_size,i,stat
  INTEGER,DIMENSION(:),ALLOCATABLE :: seed
  CHARACTER(LEN=10) :: time
  
  ! Get size of seed array
  CALL RANDOM_SEED(SIZE=seed_size)
  
  ! Allocate seed array
  ALLOCATE(seed(seed_size))
    
  ! Set seed
  CALL DATE_AND_TIME(time=time)
  
  ! Build seed vector from time/date string. Characters can't be
  ! converted to int, so the date parts are replaced by i. If seed_size > 10,
  ! then parts of the seed are non initialized
  DO i=1,MIN( LEN(time),seed_size )  
    READ(time(i:i),*,IOSTAT=stat) seed(i)
    IF(STAT .ne. 0) seed(i) = i            ! (for non integer characters)
  ENDDO
  
  CALL RANDOM_SEED(PUT=seed)

END SUBROUTINE
    
END PROGRAM
