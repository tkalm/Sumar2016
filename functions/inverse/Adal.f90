PROGRAM Adal
!----------------------------------------------------------------------
! This program writes out a matrix, calculates its inverse and
! writes the inverse.
!----------------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ----------------------------------------------
   INTEGER                           ::  i, j 
   REAL, DIMENSION(Nf,Nf)            ::  mat, inmat 
!------- Output -------------------------------------------------------
   OPEN(UNIT=10,FILE=   'origmat.dtx'         ,STATUS='NEW')
   OPEN(UNIT=11,FILE=   'invmat.dtx'          ,STATUS='NEW')
!------- Definitions --------------------------------------------------
   mat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
!         IF(j == i+1) amat(i,j) = SQRT(FLOAT(i))
      END DO
   END DO
   mat(1,1)=4
   mat(1,2)=1
   mat(2,1)=3
   mat(2,2)=2
   inmat = inv(mat)
   DO j = 1, Nf
      DO i = 1, Nf
         WRITE(10,FMT='(E15.8)') REAL(mat)
         WRITE(11,FMT='(E15.8)') REAL(inmat)
      END DO
   END DO
!------ Functions used ------------------------------------------------
   CONTAINS
   FUNCTION inv(A) RESULT(Ainv)
      REAL,     DIMENSION(Nf,Nf), intent(in)    :: A
      REAL,     DIMENSION(SIZE(A,1),SIZE(A,2))  :: Ainv
      REAL,     DIMENSION(size(A,1))            :: work                ! Work array for LAPACK
      INTEGER,  DIMENSION(size(A,1))            :: ipiv                ! Pivot indices
      INTEGER                                   :: n, info
   ! External procedures defined in LAPACK
      EXTERNAL DGETRF
      EXTERNAL DGETRI
   ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = SIZE(A,1)
   ! DGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.
      CALL DGETRF(n, n, Ainv, n, ipiv, info)
      IF (info /= 0) STOP 'Matrix is numerically singular!'
   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.
      CALL DGETRI(n, Ainv, n, ipiv, work, n, info)
      IF (info /= 0) STOP 'Matrix inversion failed!'
   end function inv
!----------------------------------------------------------------------
END PROGRAM Adal
