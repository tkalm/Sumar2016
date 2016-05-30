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
   COMPLEX, DIMENSION(Nf,Nf)  ::  rho0, H0
!------- Output ---------------------------------------------
   OPEN(UNIT=12,FILE=   'OccupationOfStates.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
   ierr = 0
 !  ALLOCATE(H0(Nf,Nf),rho0(Nf,Nf),rho1(Nf,Nf,Nt), STAT=ierr)
!----Initial matrix definitions------------------------------

   H0 = Czero
   DO i=1,Nf
      H0(i,i)=i-0.5
   END DO
   rho0 = Czero
   rho0(1,1) = 1
   rho1 = Czero
   rho1(:,:,1) = rho0

!----First we define rho1(:,:,1), the first time step--------
  
   DO i=1,10                        ! Iterations
      rho1(:,:,1) = rho0(:,:) + (lambda(rho0(:,:)) + lambda(rho1(:,:,1)))
   END DO

!----Then we define rho1(:,:,t) for the rest of time---------

   DO t=2,Nt                        ! Time grid
      DO i=1,10                     ! Iterations
        rho1(:,:,t) = rho1(:,:,t-1) + (lambda(rho1(:,:,t-1)) + lambda(rho1(:,:,t)))
      END DO
   END DO

!----And finally we write our result to the output----------
   DO i=1,4                                               ! The number of states printed
      WRITE(12,FMT='(I2,2X,(E15.8,E15.8))') 0, rho0(i,i)       ! The initial occupation of the state 
      DO t=1,Nt
         WRITE(12,FMT='(I2,2X,(E15.8,E15.8))') t, rho1(i,i,t)  ! The occupation of the state at time t
      END DO
   END DO

!   DEALLOCATE(rho0,rho1, STAT=ierr) 

!-----------------------------------------------------------
CONTAINS
FUNCTION lambda(rho) 
   COMPLEX, DIMENSION(Nf,Nf)         :: comm
   COMPLEX, DIMENSION(Nf,Nf)         :: lambda
   COMPLEX, DIMENSION(:,:), INTENT(IN)             :: rho(:,:)
   comm = matmul(H0,rho)-matmul(rho,H0)
   lambda = -ci*comm 
END FUNCTION 
!-----------------------------------------------------------
END PROGRAM Adal
