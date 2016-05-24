   MODULE PABz
   CONTAINS
   FUNCTION MATMULVGz(Af,Bf)

   USE Mod_Precision

   USE mkl95_lapack
   USE blas95
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: Nd1, Nd2, ierr
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: MATMULVGz, Af, Bf

   Nd1 = SIZE(Af,1)
   Nd2 = SIZE(Bf,1)
   ierr = 0

   ALLOCATE(MATMULVGz(Nd1,Nd2), STAT=ierr)

   MATMULVGz = (0.0_dp, 0.0_dp)
   CALL GEMM3M(Af,Bf,MATMULVGz,'N','C')

!----------------------------------------
   END FUNCTION MATMULVGz
   END MODULE PABz
