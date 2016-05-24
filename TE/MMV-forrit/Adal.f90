   PROGRAM Adal
!------------------------------------------------------------
!  Program to show the use of the mkl-soubroutine for 
!  diagonalizing a Hermitian complex matrix
!------------------------------------------------------------
!
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
!
   INTEGER         :: i, j, ierr
!------- Úttak ----------------------------------------------
!
   OPEN(UNIT=11,FILE=      'Hmat.dtx'       ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'Eigenval.dtx'      ,STATUS='NEW')
   OPEN(UNIT=13,FILE=  'Eigenvect.dtx'      ,STATUS='NEW')

!  STATUS='NEW' the program will only run if these output files
!  do not exist
!------------------------------------------------------------

   ierr = 0
   ALLOCATE(Hmat(Nf,Nf), STAT=ierr) ! 4x4 matrix

   Hmat = Czero              ! Zero-ing, Czero is defined in Mod_Init.f90

   DO j = 1, Nf              ! Matrix made
     DO i = 1, Nf
       IF(j .EQ. i) Hmat(i,j) = FLOAT(i)*pid2
       IF(j .GT. i) Hmat(i,j) = ci*pid2*0.3_dp*SQRT(FLOAT(i*j))
       WRITE(11,FMT='(2(I4,1X),2(E15.8,1X))') i, j, Hmat(i,j)
     END DO
   END DO

   ALLOCATE(Eigval(Nf),Eigvect(Nf,Nf), STAT=ierr)

   CALL HEEVR(Hmat,Eigval,UPLO,Eigvect)   ! Subroutine from MKL for finding
                                 ! eigenvalues and vectors, input is a 
                                 ! lower or upper triangular

   DO i = 1, Nf
     WRITE(12,FMT='(I4,2X,E15.8)') i, Eigval(i)  ! Eigenvalues printed
   END DO

   DO j = 1, Nf
     DO i = 1, Nf 
       WRITE(13,FMT='(I4,2X,I4,2X,E15.8,1X,E15.8)') i, j, Eigvect(i,j)
     END DO
     WRITE(13,FMT='')   ! To simplify a 3D plot with gnuplot
   END DO
 
   DEALLOCATE(Hmat,Eigval,Eigvect, STAT=ierr)

!------------------------------------------------------------
!
END PROGRAM Adal







