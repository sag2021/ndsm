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

MODULE NDSM_OPTIMIZED

  USE NDSM_ROOT

  IMPLICIT NONE

  CONTAINS

! ---------------------------------------------------------------------  
!
!>@name red_black_gauss_3D
!
!>@brief Performs one sweep of red/black Gauss-Seidel relaxation in 3D
  
SUBROUTINE red_black_gauss_3D(bcs,nx,ny,nz,x,y,z,rhs,u)

  USE NDSM_MULTIGRID_CORE,ONLY: mean

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                    ,INTENT(IN)  :: nx,ny,nz
  REAL(FP),DIMENSION(nx)         ,INTENT(IN)  :: x
  REAL(FP),DIMENSION(ny)         ,INTENT(IN)  :: y
  REAL(FP),DIMENSION(nz)         ,INTENT(IN)  :: z
  REAL(FP),DIMENSION(nx,ny,nz)   ,INTENT(IN)  :: rhs
  CHARACTER(LEN=1),DIMENSION(3,2),INTENT(IN)  :: bcs

  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(nx,ny,nz),INTENT(INOUT) :: u

  ! LOCAL
  REAL(FP)    :: unew,umean
  REAL(FP)    :: hx,hy,hz,w1,wx,wy,wz
  INTEGER(IT) :: i,j,k
  INTEGER(IT) :: xl,yl,zl,xh,yh,zh
  INTEGER(IT) :: lb(3),ub(3),nsize

  ! LOCAL LITS.
  INTEGER(IT),PARAMETER :: i2 = 2
    
  ! Get array bounds
  lb = [1 ,1 ,1]
  ub = [nx,ny,nz]
  
  ! Total size
  nsize = SIZE(u)

  ! Modify bounds based on boundary conditions
  WHERE(bcs(:,1) == "D") lb = lb + 1
  WHERE(bcs(:,2) == "D") ub = ub - 1

  ! Grid has uniform spacing
  hx = x(2)-x(1)
  hy = y(2)-y(1)
  hz = z(2)-z(1)

  ! Second-order centred difference weights for off-diagonal terms,
  ! e.g. u(i+1,j,k) etc.
  !
  !
  wx = REAL(1,FP)/hx**2
  wy = REAL(1,FP)/hy**2
  wz = REAL(1,FP)/hz**2
 
  ! Second-order centred difference weights for diagonal term
  !
  w1 = 2*(wx+wy+wz)
  w1 = REAL(1,FP)/w1     
       
  ! ========================
  ! RELAXATION 
  ! ========================
  !
  ! Red-Black Gauss-Sidel
  !
  
  !$OMP  PARALLEL DO PRIVATE(i,j,k,xl,xh,yl,yh,zl,zh,unew) 
  DO k=lb(3),ub(3)
    DO j=lb(2),ub(2)
      DO i=lb(1)+MOD(j+MOD(k,i2),i2),ub(1),2

        ! Get points surrounding (i,j,k)
        xl = i-1; xh = i+1
        yl = j-1; yh = j+1    
        zl = k-1; zh = k+1
          
        IF(xl < 1 ) xl = 2
        IF(xh > nx) xh = nx-1
      
        IF(yl < 1 ) yl = 2
        IF(yh > ny) yh = ny-1
            
        IF(zl < 1 ) zl = 2
        IF(zh > nz) zh = nz-1
          
        ! Perform relaxation 
        unew =  (u(xh,j,k) +  u(xl,j,k))*wx &
             +  (u(i,yh,k) +  u(i,yl,k))*wy &
             +  (u(i,j,zh) +  u(i,j,zl))*wz &
             -  rhs(i,j,k)
               
        ! Update 
        u(i,j,k) = w1*unew
        
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  !$OMP  PARALLEL DO PRIVATE(i,j,k,xl,xh,yl,yh,zl,zh,unew) 
  DO k=lb(3),ub(3)
    DO j=lb(2),ub(2)
      DO i=lb(1)+MOD(j+1+MOD(k,i2),i2),ub(1),2 

        ! Get points surrounding (i,j,k)
        xl = i-1; xh = i+1
        yl = j-1; yh = j+1    
        zl = k-1; zh = k+1
          
        IF(xl < 1 ) xl = 2
        IF(xh > nx) xh = nx-1
      
        IF(yl < 1 ) yl = 2
        IF(yh > ny) yh = ny-1
            
        IF(zl < 1 ) zl = 2
        IF(zh > nz) zh = nz-1
          
        ! Perform relaxation 
        unew =  (u(xh,j,k) +  u(xl,j,k))*wx &
             +  (u(i,yh,k) +  u(i,yl,k))*wy &
             +  (u(i,j,zh) +  u(i,j,zl))*wz &
             -  rhs(i,j,k)
                       
        ! Update 
        u(i,j,k) = w1*unew
        
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Compute mean. Used to pick unique Neumann 
  ! solution
  !
  !
  IF(ALL(bcs== "N")) THEN  

    ! Compute mean
    umean = 0
    umean = mean(nsize,u)
    
    !$OMP  PARALLEL DO PRIVATE(i,j,k) 
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)        
          u(i,j,k) = u(i,j,k) - umean        
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO 

  ENDIF
  
END SUBROUTINE

! ---------------------------------------------------------------------  
!
!>@name red_black_gauss_2D
!
!>@brief Performs one sweep of red/black Gauss-Seidel relaxation in 2D
  
SUBROUTINE red_black_gauss_2D(bcs,nx,ny,x,y,rhs,u)

  USE NDSM_MULTIGRID_CORE,ONLY: mean

  IMPLICIT NONE
  
  ! INPUT
  INTEGER(IT)                    ,INTENT(IN)  :: nx,ny
  REAL(FP),DIMENSION(nx)         ,INTENT(IN)  :: x
  REAL(FP),DIMENSION(ny)         ,INTENT(IN)  :: y
  REAL(FP),DIMENSION(nx,ny)      ,INTENT(IN)  :: rhs
  CHARACTER(LEN=1),DIMENSION(2,2),INTENT(IN)  :: bcs

  ! INPUT/OUTPUT
  REAL(FP),DIMENSION(nx,ny),INTENT(INOUT) :: u

  ! LOCAL
  REAL(FP) :: unew,umean
  REAL(FP) :: hx,hy,w1,wx,wy
  INTEGER(IT)  :: i,j
  INTEGER(IT)  :: xl,yl,xh,yh
  INTEGER(IT)  :: lb(2),ub(2),nsize

  ! LOCAL LITS.
  INTEGER(IT),PARAMETER :: i2 = 2
    
  ! Get array bounds
  lb = LBOUND(u)
  ub = UBOUND(u)
  
  ! Total size
  nsize = SIZE(u)

  ! Modify bounds based on boundary conditions
  WHERE(bcs(:,1) == "D") lb = lb + 1
  WHERE(bcs(:,2) == "D") ub = ub - 1

  ! Grid has uniform spacing
  hx = x(2)-x(1)
  hy = y(2)-y(1)

  ! Second-order centred difference weights for off-diagonal terms,
  ! e.g. u(i+1,j,k) etc.
  !
  !
  wx = REAL(1,FP)/hx**2
  wy = REAL(1,FP)/hy**2
 
  ! Second-order centred difference weights for diagonal term
  !
  w1 = REAL(2,FP)/hx**2 &  
     + REAL(2,FP)/hy**2  
    
  w1 = REAL(1,FP)/w1
            
  ! ========================
  ! RELAXATION 
  ! ========================
  !
  ! Relax using weighted Jacobi method
  !
  
  !$OMP  PARALLEL DO PRIVATE(i,j,xl,xh,yl,yh,unew) 
  DO j=lb(2),ub(2)
    DO i=lb(1)+MOD(j,i2),ub(1),2

        ! Get points surrounding (i,j)
        xl = i-1; xh = i+1
        yl = j-1; yh = j+1    
          
        IF(xl < 1 ) xl = 2
        IF(xh > nx) xh = nx-1
      
        IF(yl < 1 ) yl = 2
        IF(yh > ny) yh = ny-1
            
        ! Perform relaxation 
        unew =  (u(xh,j) +  u(xl,j))*wx &
             +  (u(i,yh) +  u(i,yl))*wy &
             -  rhs(i,j)
               
        ! Update 
        u(i,j) = w1*unew
        
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  !$OMP  PARALLEL DO PRIVATE(i,j,xl,xh,yl,yh,unew) 
  DO j=lb(2),ub(2)
    DO i=lb(1)+MOD(j+1,i2),ub(1),2 

        ! Get points surrounding (i,j,k)
        xl = i-1; xh = i+1
        yl = j-1; yh = j+1    
          
        IF(xl < 1 ) xl = 2
        IF(xh > nx) xh = nx-1
      
        IF(yl < 1 ) yl = 2
        IF(yh > ny) yh = ny-1
                      
        ! Perform relaxation 
        unew =  (u(xh,j) +  u(xl,j))*wx &
             +  (u(i,yh) +  u(i,yl))*wy &
             -  rhs(i,j)
                       
        ! Update 
        u(i,j) = w1*unew
        
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
  ! Compute mean. Used to pick unique Neumann 
  ! solution
  !
  !
  umean = 0
  
  IF(ALL(bcs== "N")) THEN  
    umean = mean(nsize,u)
  ENDIF
    
  ! 
  ! Copy back to ut and subtract mean
  !
  !$OMP  PARALLEL DO PRIVATE(i,j) 
  DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)        
        u(i,j) = u(i,j) - umean        
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 
    
END SUBROUTINE

! ---------------------------------------------------------------------
!
!>@name poisson_residual_3D
!
!>@brief Compute residual for 3D Poisson equation 
!
!>@details
!! A second-order stencil is used to compute the derivatives
!!

SUBROUTINE poisson_residual_3D(bcs,nx,ny,nz,x,y,z,rhs,u,r)

  IMPLICIT NONE

  ! INPUT
  INTEGER(IT)                    ,INTENT(IN)  :: nx,ny,nz
  REAL(FP),DIMENSION(nx)         ,INTENT(IN)  :: x
  REAL(FP),DIMENSION(ny)         ,INTENT(IN)  :: y
  REAL(FP),DIMENSION(nz)         ,INTENT(IN)  :: z
  REAL(FP),DIMENSION(nx,ny,nz)   ,INTENT(IN)  :: rhs,u
  CHARACTER(LEN=1),DIMENSION(3,2),INTENT(IN)  :: bcs

  ! OUTPUT
  REAL(FP),DIMENSION(nx,ny,nz),INTENT(OUT) :: r

  ! LOCAL
  INTEGER(IT) :: i,j,k
  INTEGER(IT) :: xl,yl,zl,xh,yh,zh
  INTEGER(IT) :: lb(3),ub(3)
  REAL(FP)    :: hx,hy,hz,wx,wy,wz,wc
    
  ! Get array bounds
  lb = LBOUND(u)
  ub = UBOUND(u)

  ! Modify bounds based on boundary conditions
  WHERE(bcs(:,1) == "D") lb = lb + 1
  WHERE(bcs(:,2) == "D") ub = ub - 1
  
  ! Grid has uniform spacing
  hx = x(2)-x(1)
  hy = y(2)-y(1)
  hz = z(2)-z(1)
  
  ! Weights
  wx = REAL(1,FP)/hx**2
  wy = REAL(1,FP)/hy**2
  wz = REAL(1,FP)/hz**2
  wc = 2*(wx+wy+wz)
  
  ! ==================
  ! Calculate residue
  ! ==================
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        r(i,j,k) = 0
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO 
  
  !$OMP PARALLEL DO PRIVATE(i,j,k,xl,xh,yl,yh,zl,zh)
  DO k=lb(3),ub(3)
    DO j=lb(2),ub(2)
      DO i=lb(1),ub(1)

      ! Get points surrounding (i,j,k)
      xl = i-1
      xh = i+1
      yl = j-1
      yh = j+1
      zl = k-1
      zh = k+1
            
      IF(xl < 1 ) xl = 2
      IF(xh > nx) xh = nx-1
      
      IF(yl < 1 ) yl = 2
      IF(yh > ny) yh = ny-1
            
      IF(zl < 1 ) zl = 2
      IF(zh > nz) zh = nz-1
       
      ! 
      ! Comptute R = L[u] - rhs
      !
      r(i,j,k) =  (u(xl,j ,k ) + u(xh,j ,k ))*wx  &
               +  (u(i ,yl,k ) + u(i ,yh,k ))*wy  &
               +  (u(i ,j ,zl) + u(i ,j ,zh))*wz  &
               -  rhs(i,j,k)   - u(i,j,k)*wc
              
      ! Invert       
      r(i,j,k) = -r(i,j,k)
                            
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
   
  ! Homogeneou Dirichlet boundary conditions are enforced
  ! exactly and so have zero residue
  IF(bcs(1,1) == "D") r(1,:,:)    = 0
  IF(bcs(2,1) == "D") r(:,1,:)    = 0
  IF(bcs(3,1) == "D") r(:,:,1)    = 0

  IF(bcs(1,2) == "D") r(nx,: ,: ) = 0
  IF(bcs(2,2) == "D") r(: ,ny,: ) = 0
  IF(bcs(3,2) == "D") r(: ,: ,nz) = 0
  
END SUBROUTINE

END MODULE
