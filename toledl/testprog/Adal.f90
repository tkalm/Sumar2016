PROGRAM Adal
!----------------------------------------------------------------------
! This program finds the eigenvalues and vectors of the density matrix
! rho in the Liouville-von Neumann equation in Liouville space. 
!----------------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Precision_Mod         ! Module for setting double precision
   USE Init_Mod              ! Initial values
   USE Fields_Mod            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   USE KronAB                ! Function KronM:  Kronecker product
   USE KronATB               ! Function KronTM: First matrix transposed
   USE PAB                   ! Function MATMULVG: Product of matrices
   USE PcAB                  ! Function MATMULVGc: First matrix transp 
   USE PABz                  ! Function MATMULVGz: Latter matrix transp

   IMPLICIT NONE
!------- Local variables ----------------------------------------------
   INTEGER                           ::  i, j, l, ierr, ssloc
   REAL(KIND=dp)                     ::  Nsum, Eave, time
   COMPLEX(KIND=dp)                  ::  cnorm
!------- Output -------------------------------------------------------
   OPEN(UNIT=12,FILE=   'eigvallog.dtx'     ,STATUS='NEW')
!------- Definitions --------------------------------------------------
   ierr = 0
   ALLOCATE(vec(Nf), mat(Nf,Nf), locvec(Nf), STAT=ierr)
   
   mat = Czero
   
!   DO j = 1, Nf
      DO i = 1, Nf
         vec(i) = i
      END DO
!   END DO
   
   vec(4) = -5

   locvec = MINLOC(vec,1)

!   DO j = 1, Nf
      DO i = 1, Nf
         WRITE(12,FMT="(I6,1X,E15.8)") i, vec(i) 
      END DO
      WRITE(12,FMT="(1X)") 
      DO i = 1, Nf
         WRITE(12,FMT="(I6,1X,E15.8)") i, locvec(i) 
      END DO
!   END DO

   DEALLOCATE(mat,  locvec, STAT=ierr)
!------------------------------------------------------------------------------
END PROGRAM Adal
