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
   COMPLEX, DIMENSION(Nf,Nf,0:Nt)    ::  rho
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, Ht
!------- Output ---------------------------------------------
   OPEN(UNIT=12,FILE=   'OccupationOfStates.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
! We begin be defining the hamiltonian and the initial values
! of the density matrix rho. 

   rho(:,:,:) = Czero
   rho(1,1,0) = 1               ! Initial value of rho
   
   H = Czero                    ! We define H as H0 and define Ht and then take the sum
   Ht = Czero
   DO j = 1, Nf 
     DO i = 1, Nf
       IF(j == i) H(i,j) = i-0.5
       IF(ABS(j-i) == 1) Ht(i,j) = SQRT(FLOAT((j+i-1)/2))
     END DO
   END DO
   H=H+Ht

!------------------------------------------------------------
! Then we calculate the entries of rho for other values of time
! by iteration of our approximation of the L-vN equation. 

   DO j=1,Nt                        ! Time grid
      DO i=1,10                     ! Iterations
        rho(:,:,j) = rho(:,:,j-1)
        rho(:,:,j) = rho(:,:,j-1) + 0.01208897*(lambda(rho(:,:,j-1)) + lambda(rho(:,:,j)))
      END DO
   END DO

!------------------------------------------------------------
! And finally we write our results to a file

   DO i=1,5                                                  ! The number of states printed
      DO j=0,Nt
         WRITE(12,FMT='(E15.8,2X,E15.8)') FLOAT(j)*0.1, REAL(rho(i,i,j))  ! The occupation of the state at time t
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
