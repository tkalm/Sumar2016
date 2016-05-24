   MODULE Kronvec
   CONTAINS
   FUNCTION KronV(Af)

!  Kronecker vector constructed for a complex square matrix

   USE Mod_Precision

   USE mkl95_lapack
   USE blas95
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: Nd1, Nd2, ierr, NdK, i, j
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: Af, Bf
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)   :: KronV

   Nd1 = SIZE(Af,1)
   Nd2 = Nd1**2
   ierr = 0

 
   ALLOCATE(KronV(Nd2), STAT=ierr)

   KronV = (0.0_dp, 0.0_dp)

!$OMP PARALLEL DO SHARED(Af, KronV) SCHEDULE(DYNAMIC) 
   DO j = 1, Nd1
     DO i = 1, Nd1
       KronV(Nd1*(j-1)+i) =  Af(i,j)
     END DO
   END DO
!$OMP END PARALLEL DO

!----------------------------------------
   END FUNCTION KronV
   END MODULE Kronvec
