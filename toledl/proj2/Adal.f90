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
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, Ht, rho, rhon
!------- Output ---------------------------------------------
   OPEN(UNIT=11,FILE=   'State0.dtx'      ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'State1.dtx'      ,STATUS='NEW')
   OPEN(UNIT=13,FILE=   'State2.dtx'      ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'State3.dtx'      ,STATUS='NEW')
   OPEN(UNIT=15,FILE=   'State4.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
! We begin be defining the hamiltonian and the initial values
! of the density matrix rho. 

   H = Czero                    ! We define H as H0 and define Ht and then take the sum
   Ht = Czero
   DO j = 1, Nf 
     DO i = 1, Nf
       IF(j == i) H(i,j) = i-0.5
       IF(ABS(j-i) == 1) Ht(i,j) = SQRT(FLOAT((j+i-1)/2))
     END DO
   END DO
   H=H+Ht

   rho(:,:) = Czero
   rho(1,1) = 1                 ! Initial value of rho
   rhon = rho
! We begin by printing the initial state
   DO i=1,5                                                  ! The number of states printed
      WRITE(10+i,FMT='(E15.8,2X,E15.8)') 0, REAL(rho(i,i))   ! The occupation of the state at time t=0
   END DO

!------------------------------------------------------------
! Then we calculate and print the entries of rho for other 
! values of time by iteration of our approximation 
! of the L-vN equation. 

   DO j=1,Nt                        ! Time grid
      DO i=1,10                     ! Iterations
        rhon(:,:) = rho(:,:) + 0.006044485*(lambda(rho(:,:)) + lambda(rhon(:,:)))
      END DO
      rho=rhon
      DO i=1,5                                                             ! The number of states printed
         WRITE(10+i,FMT='(E15.8,2X,E15.8)') FLOAT(j)*0.05, REAL(rho(i,i))   ! The occupation of the state at time t=0
      END DO
   END DO

!-----------------------------------------------------------
! The function lambda is the superoperator that calculates
! the commutator of H and rho in the L-vN eq. 

   CONTAINS
   FUNCTION lambda(mat) 
      COMPLEX, DIMENSION(Nf,Nf)                :: lambda
      COMPLEX, DIMENSION(:,:), INTENT(IN)      :: mat(:,:)
      lambda = matmul(H,mat)-matmul(mat,H)
      lambda = -ci*lambda 
   END FUNCTION

!-----------------------------------------------------------
END PROGRAM Adal
