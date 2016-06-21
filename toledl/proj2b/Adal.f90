   PROGRAM Adal
!------------------------------------------------------------
! This program calculates the expectation values <x> and 
! <x*x> for the system H = H0 + Ht (from project 2a)
!------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
   INTEGER                           ::  i, j 
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, rho, rhon1, rhon2, xmat, x2mat
   REAL                              ::  hbaromega, delt, alpha
!------- Output ---------------------------------------------
   OPEN(UNIT=11,FILE=   'expx.dtx'        ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'expxx.dtx'       ,STATUS='NEW')
!------------------------------------------------------------
! Define H, x, x*x and the initial values of rho. 
   H = Czero                                                         ! Define H
   xmat = Czero
   DO j = 1, Nf 
     DO i = 1, Nf
       IF(j == i) H(i,j) = i-0.5
       IF(ABS(j-i) == 1) H(i,j) = SQRT(FLOAT((i+j-1)/2))
       IF(ABS(j-i) == 1) xmat(i,j) = SQRT(FLOAT((i+j-1)/2))
     END DO
   END DO
   xmat = xmat/SQRT(FLOAT(2))
   x2mat = matmul(xmat,xmat)
   rho = Czero
   rho(1,1) = 1                                                      ! Initial value of rho
   rhon1 = rho
! Print the initial expectation values
   WRITE(11,FMT='(E15.8,2X,E15.8)') 0, REAL(tr(matmul(rho,xmat)))    ! The expectation of x/a at time t=0 
   WRITE(12,FMT='(E15.8,2X,E15.8)') 0, REAL(tr(matmul(rho,x2mat)))   ! The expectation of (x/a)*(x/a) at time t=0
! Define the value of hbar*omega in our system, the timestep
! used and the constant that appears in our equation. 
   hbaromega  = 1E-3                                                 ! hbar*omega in eV
   delt       = 1E-3                                                 ! delta t in ps 
   alpha      = hbaromega*delt/(2*hbar)                              ! The constant in our equation below
!------------------------------------------------------------
! Calculate and print the expectation values for t>0 
   DO j=1,Nt/delt                                                    ! Time grid
      DO                                                             ! Iterations
        rhon2 = rho + alpha*(lambda(rho) + lambda(rhon1))
        IF(err(rhon1,rhon2)<1E-1) EXIT
        rhon1 = rhon2
      END DO
      rho=rhon2
      WRITE(11,FMT='(E15.8,2X,E15.8)') FLOAT(j)*delt, REAL(tr(matmul(rho,xmat)))
      WRITE(12,FMT='(E15.8,2X,E15.8)') FLOAT(j)*delt, REAL(tr(matmul(rho,x2mat)))
   END DO
!-----------------------------------------------------------
! Functions used 
   CONTAINS
   FUNCTION lambda(mat)                                       ! Calculates the commutator in L-vN 
      COMPLEX, DIMENSION(Nf,Nf)                :: lambda
      COMPLEX, DIMENSION(:,:), INTENT(IN)      :: mat(:,:)
      lambda = -ci*(matmul(H,mat)-matmul(mat,H))
   END FUNCTION lambda

   FUNCTION err(x1,x2)                                        ! Calculates the error
      REAL                                :: err 
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: x1,x2
      err = SQRT(ABS(SUM(x2-x1)))
   END FUNCTION err

   FUNCTION tr(matrix)                                        ! Calculates the trace of a matrix
      INTEGER                             :: m
      COMPLEX                             :: tr
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: matrix 
      tr = Czero 
      DO m = 1, Nf
         tr = tr+matrix(m,m)
      END DO
   END FUNCTION tr
!-----------------------------------------------------------
END PROGRAM Adal
