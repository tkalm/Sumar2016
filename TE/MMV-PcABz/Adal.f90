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
   USE blas95                ! Subroutines

   USE PAB
   USE PcAB
   USE PABz

   IMPLICIT NONE
!------- Local variables ------------------------------------
!
   INTEGER         :: i, j, ierr, info
!------- Úttak ----------------------------------------------
!
   OPEN(UNIT=11,FILE=      'Hmat.dtx'       ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'Eigenval.dtx'      ,STATUS='NEW')
   OPEN(UNIT=13,FILE=  'Eigenvect.dtx'      ,STATUS='NEW')
   OPEN(UNIT=14,FILE=         'Hn.dtx'      ,STATUS='NEW')

!  STATUS='NEW' the program will only run if these output files
!  do not exist
!------------------------------------------------------------

   ierr = 0
   ALLOCATE(Hmat(Nf,Nf), Hn(Nf,Nf), VD(Nf,Nf), Henn(Nf,Nf), STAT=ierr)
   ALLOCATE(xmat(Nf,Nf), STAT=ierr)

   xmat = Czero             ! Set x-matrix as zero
   DO j=1, Nf               ! Define entries of x-matrix
      DO i=1, Nf
         IF(ABS(j-i) == 1) xmat(i,j) = SQRT(FLOAT(i+j+1))*0.5
      END DO
   END DO
   xmat = matmul(xmat,xmat)
   xmat = matmul(xmat,xmat) ! Take x to the fourth power

   Hmat = Czero              ! Set H matrix as zero
   DO j = 1, Nf              ! Matrix made
     DO i = 1, Nf
       IF(j .EQ. i) Hmat(i,j) = i-0.5
     END DO
   END DO
   Hmat = Hmat + xmat      ! Add the xmatrix
   DO j = 1, Nf              ! Print matrix 
     DO i = 1, Nf
       WRITE(11,FMT='(I4,2X,I4,2X,E15.8,1X,E15.8)') i, j, Hmat(i,j)
     END DO
   END DO

  Henn = Hmat

   ALLOCATE(Eigval(Nf), STAT=ierr)

   info = 0
   CALL HEEVD(Hmat,Eigval,JOB,UPLO,info)   ! Subroutine from MKL for finding
                                       ! eigenvalues and vectors, input is a 
   VD = Hmat                           ! lower or upper triangular

   DO i = 1, Nf
     WRITE(12,FMT='(I4,2X,E15.8)') i, Eigval(i)  ! Eigenvalues printed
   END DO

   DO j = 1, Nf
     DO i = 1, Nf 
       WRITE(13,FMT='(I4,2X,I4,2X,E15.8,1X,E15.8)') i, j, VD(i,j)
     END DO
     WRITE(13,FMT='')   ! To simplify a 3D plot with gnuplot
   END DO

   Hn = MATMULVG(MATMULVGc(VD,Henn),VD)
 

   DO j = 1, Nf
     DO i = 1, Nf 
       WRITE(14,FMT='(I4,2X,I4,2X,E15.8,1X,E15.8)') i, j, Hn(i,j)
     END DO
     WRITE(14,FMT='')   ! To simplify a 3D plot with gnuplot
   END DO

   DEALLOCATE(Hmat, Eigval, Hn, Henn, VD, STAT=ierr)

!------------------------------------------------------------

END PROGRAM Adal







