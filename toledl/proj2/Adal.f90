   PROGRAM Adal
!------------------------------------------------------------
!  
!------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
   INTEGER         :: i, j, ierr, lambda, n
!------- Output ---------------------------------------------
   OPEN(UNIT=12,FILE=   'Eigenval.dtx'      ,STATUS='NEW')
!------------------------------------------------------------
   ierr = 0
   ALLOCATE(H0(Nf,Nf),V(Nf,Nf),Hmat(Nf,Nf), STAT=ierr)
!------------------------------------------------------------
! First we define H0 and V.

   V = Czero                                 ! Set V as the zero matrix
   DO j=1, Nf                                ! Set V as x/a in the basis of H0
      DO i=1, Nf
         IF(ABS(j-i) == 1) V(i,j) = SQRT(FLOAT(i+j+1))*0.5
      END DO
   END DO
   V = matmul(V,V)                           ! Take V to the second power
   V = matmul(V,V)                           ! Take V to the fourth power

   H0 = Czero                                ! Set H0 as the zero matrix
   DO j = 1, Nf                              ! Define the entries of H0 
     DO i = 1, Nf
       IF(j == i) H0(i,j) = i-0.5
     END DO
   END DO

!------------------------------------------------------------
! Here we calculate the eigenvalues for different
! values of lambda.

   DO n=1, 6                                                      ! n will be the number of 
                                                                  ! eigenvalues we print for each lambda
      DO lambda=0, 120
         Hmat = H0 + V*0.01*lambda                                ! Define H for this lambda 

         ALLOCATE(Eigval(Nf),Eigvect(Nf,Nf), STAT=ierr)
         CALL HEEVR(Hmat,Eigval,UPLO,Eigvect)                     ! Subroutine from MKL for finding
                                                                  ! eigenvalues and vectors, input is a 
                                                                  ! lower or upper triangular
                                                                  
         WRITE(12,FMT='(E15.8,2X,E15.8)') 0.01*lambda, Eigval(n)  ! Eigenvalues printed (along with
                                                                  ! the value of lambda)
         DEALLOCATE(Eigval,Eigvect, STAT=ierr)
      END DO
   END DO
   
   DEALLOCATE(H0,V,Hmat, STAT=ierr)

!-----------------------------------------------------------
END PROGRAM Adal
