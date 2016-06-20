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
   OPEN(UNIT=11,FILE=   'exp(x).dtx'        ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'exp(x*x).dtx'      ,STATUS='NEW')
!------------------------------------------------------------
! Define the hamiltonian and the initial values of rho. 
   H = Czero                                                  ! Define H
   DO j = 1, Nf 
     DO i = 1, Nf
       IF(j == i) H(i,j) = i-0.5
       IF(ABS(j-i) == 1) H(i,j) = SQRT(FLOAT((j+i-1)/2))
     END DO
   END DO
   rho = Czero
   rho(1,1) = 1                                               ! Initial value of rho
   rhon1 = rho

   xmat = (1/SQRT(FLOAT(2)))*Ht
   x2mat = matmul(xmat,xmat)
! We begin by printing the initial expectation values
   WRITE(11,FMT='(E15.8,2X,E15.8)') 0, REAL(tr(matmul(rho,xmat)))    ! The expectation of x/a at time t=0 
   WRITE(12,FMT='(E15.8,2X,E15.8)') 0, REAL(tr(matmul(rho,x2mat)))   ! The expectation of (x/a)*(x/a) at time t=0
! Define the value of hbar*omega in our system, the timestep
! used and a constant that appears in our equation. 
   hbaromega  = 1E-3                                          ! hbar*omega in eV
   delt       = 1E-3                                          ! delta t in ps 
   alpha      = hbaromega*delt/(2*hbar)                       ! The constant in our equation below
!------------------------------------------------------------
! Then we calculate and print the entries of rho for other 

   DO j=1,Nt                        ! Time grid
      DO                            ! Iterations
        rhon2(:,:) = rho(:,:) + alpha*(lambda(rho(:,:)) + lambda(rhon1(:,:)))
        IF(err(rhon1,rhon2)<1E-7) EXIT
        rhon1 = rhon2
      END DO
      rho=rhon2
      WRITE(11,FMT='(E15.8,2X,E15.8)') FLOAT(j)*0.0001, REAL(tr(matmul(rho,xmat)))    ! The expectation of x/a at time t=0 
      WRITE(12,FMT='(E15.8,2X,E15.8)') FLOAT(j)*0.0001, REAL(tr(matmul(rho,x2mat)))   ! The expectation of (x/a)*(x/a) at time t=0
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
