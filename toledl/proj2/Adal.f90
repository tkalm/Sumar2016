   PROGRAM Adal
!------------------------------------------------------------
! This program calculates the occupation of chosen states
! of the harmonic oscillator given a constant potential 
! is turned on at time t=0. 
!------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
   INTEGER         :: t, i, ierr
   COMPLEX, DIMENSION(Nf,Nf,Nt)  ::  rho1
   COMPLEX, DIMENSION(Nf,Nf)  ::  rho0, H0, V
!------- Output ---------------------------------------------
   OPEN(UNIT=12,FILE=   'OccupationOfStates.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
   ierr = 0
!----Initial matrix definitions------------------------------
   rho0 = Czero
   rho0(2,2) = 1
   rho1 = Czero
   rho1(:,:,1) = rho0
   
   V = Czero                                 ! Set V as the zero matrix
   DO t=1, Nf                                ! Set V as a* + a in energy basis 
      DO i=1, Nf
         IF(ABS(t-i) == 1) V(i,t) = SQRT(FLOAT((t+i-1)/2))
      END DO
   END DO
   H0 = Czero                                ! Set H0 as the zero matrix
   DO t = 1, Nf                              ! Define the entries of H0 
     DO i = 1, Nf
       IF(t == i) H0(i,t) = i-0.5
     END DO
   END DO
   H0=H0+V
!----First we define rho1(:,:,1), the first time step--------
   DO i=1,4                        ! Iterations
      rho1(:,:,1) = rho0(:,:) + 0.1208897*(lambda(rho0(:,:)) + lambda(rho1(:,:,1)))
   END DO
!----Then we define rho1(:,:,t) for the rest of time---------
   DO t=2,Nt                        ! Time grid
      DO i=1,10                     ! Iterations
        rho1(:,:,t) = rho1(:,:,t-1)
        rho1(:,:,t) = rho1(:,:,t-1) + 0.1208897*(lambda(rho1(:,:,t-1)) + lambda(rho1(:,:,t)))
      END DO
   END DO
!----And finally we write our result to the output----------
   DO i=1,10                                                 ! The number of states printed
      WRITE(12,FMT='(I2,2X,E15.8)') 0, REAL(rho0(i,i))       ! The initial occupation of the state 
      DO t=1,Nt
         WRITE(12,FMT='(I2,2X,E15.8)') t, REAL(rho1(i,i,t))  ! The occupation of the state at time t
      END DO
   END DO
!-----------------------------------------------------------
   CONTAINS
   FUNCTION lambda(rho) 
      COMPLEX, DIMENSION(Nf,Nf)                :: comm
      COMPLEX, DIMENSION(Nf,Nf)                :: lambda
      COMPLEX, DIMENSION(:,:), INTENT(IN)      :: rho(:,:)
      comm = matmul(H0,rho)-matmul(rho,H0)
      lambda = -ci*comm 
   END FUNCTION 
!-----------------------------------------------------------
END PROGRAM Adal
