   PROGRAM Adal
!------------------------------------------------------------
!  Let H0 be the harmonic oscillator and V=h_bar*omega*(x/a)^4
!  where 'a' is the natural length of the system.
!  This program calculates how the energy of each steady state 
!  of the system H = H0 + lambda*V increases with lambda.
!  This is achieved by evaluating H in the basis of H0
!  and then calculating the eigenvalues of H for
!  different values of lambda. 
!------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
   INTEGER         :: i, j, ierr, n
   REAL (KIND=dp)  :: lambda
!------- Output ---------------------------------------------
   OPEN(UNIT=11,FILE=   'State0.dtx'      ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'State1.dtx'      ,STATUS='NEW')
   OPEN(UNIT=13,FILE=   'State2.dtx'      ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'State3.dtx'      ,STATUS='NEW')
   OPEN(UNIT=15,FILE=   'State4.dtx'      ,STATUS='NEW')
   OPEN(UNIT=16,FILE=   'State5.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
   ierr = 0
   ALLOCATE(H0(Nf,Nf),V(Nf,Nf),Hmat(Nf,Nf), STAT=ierr)
!------------------------------------------------------------
! First we define H0 and V.

   V = Czero
   DO j = 1, Nf 
     DO i = 1, Nf
       IF(ABS(j-i) == 1) V(i,j) = SQRT(FLOAT((j+i-1)/2))
     END DO
   END DO
   V = (1/SQRT(FLOAT(2)))*V
   V = matmul(V,V)
   V = matmul(V,V)

   H0 = Czero                                ! Set H0 as the zero matrix
   DO j = 1, Nf                              ! Define the entries of H0 
     DO i = 1, Nf
       IF(j == i) H0(i,j) = i-0.5
     END DO
   END DO

!------------------------------------------------------------
! Here we calculate the eigenvalues for different
! values of lambda.

   DO n=1, 6                                                      ! n is the number of eigenvalues we print 
      DO lambda=0, 1.2, 0.01
         Hmat = (H0 + V*lambda)                                   ! Define H for this lambda 
         ALLOCATE(Eigval(Nf),Eigvect(Nf,Nf), STAT=ierr)
         CALL HEEVR(Hmat,Eigval,UPLO,Eigvect)                     ! Finds the eigenvalues (and vectors) 
         WRITE(10+n,FMT='(E15.8,2X,E15.8)') lambda, Eigval(n)       ! Eigenvalues printed (along with lambda)
         DEALLOCATE(Eigval,Eigvect, STAT=ierr)
      END DO
   END DO
   
   DEALLOCATE(H0,V,Hmat, STAT=ierr)

!-----------------------------------------------------------
END PROGRAM Adal
