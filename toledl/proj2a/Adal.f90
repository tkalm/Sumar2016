   PROGRAM Adal
!------------------------------------------------------------
! This program calculates the occupation over time  of chosen 
! states of a system with a hamiltonian H0 + V(t) where H0 is 
! the harmonic oscillator and V(t) is 0 for t<0 but (a+a*) 
! for t>=0 where a is the ladder operator. 
! This is done by a Crank-Nicolson approximation of the 
! Liouville-von Neumann equation. 
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
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, rho, rhon1, rhon2
   REAL                              ::  hbaromega, delt, alpha
!------- Output ---------------------------------------------
   OPEN(UNIT=10,FILE=   'trace.dtx'      ,STATUS='NEW')
   OPEN(UNIT=11,FILE=   'State0.dtx'      ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'State1.dtx'      ,STATUS='NEW')
   OPEN(UNIT=13,FILE=   'State2.dtx'      ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'State3.dtx'      ,STATUS='NEW')
   OPEN(UNIT=15,FILE=   'State4.dtx'      ,STATUS='NEW')
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

! Print the initial state
   DO i=1,5                                                   ! The number of states printed
      WRITE(10+i,FMT='(E15.8,2X,E15.8)') 0, REAL(rho(i,i))    ! The occupation of states at time t=0
   END DO
   WRITE(10,FMT='(E15.8,2X,E15.8)') 0, REAL(tr(rho))          ! The trace of rho at time t=0

! Define the value of hbar*omega in our system, the timestep
! used and a constant that appears in our equation. 
   hbaromega  = 1E-3                                          ! hbar*omega in eV
   delt       = 1E-3                                          ! delta t in ps 
   alpha      = hbaromega*delt/(2*hbar)                       ! The constant in our equation below
!------------------------------------------------------------
! Calculate and print entries of rho for t>0 by iteration 
! of the Crank-Nicolson approx. the of  L-vN equation. 
   DO j=1,Nt                                                                  ! Time grid
      DO                                                                      ! Iterations
        rhon2(:,:) = rho(:,:) + alpha*(lambda(rho(:,:)) + lambda(rhon1(:,:)))
        IF(err(rhon1,rhon2)<1E-4) EXIT
        rhon1 = rhon2
      END DO
      rho=rhon2
      DO i=1,5                                                                ! The number of states printed
        WRITE(10+i,FMT='(E15.8,2X,E15.8)') FLOAT(j)*delt, REAL(rho(i,i))   ! The occupation of states at time t=j
      END DO
      WRITE(10,FMT='(E15.8,2X,E15.8)') FLOAT(j)*delt, REAL(tr(rho))         ! The trace of rho at time t=j
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
