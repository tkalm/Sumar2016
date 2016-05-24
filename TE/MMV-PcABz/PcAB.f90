   MODULE PcAB
   CONTAINS
   FUNCTION MATMULVGc(Af,Bf)

   USE Mod_Precision

   USE mkl95_lapack
   USE blas95
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: Nd1, Nd2, ierr
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: MATMULVGc, Af, Bf

   Nd1 = SIZE(Af,2)
   Nd2 = SIZE(Bf,2)
   ierr = 0

   ALLOCATE(MATMULVGc(Nd1,Nd2), STAT=ierr)

   MATMULVGc = (0.0_dp, 0.0_dp)
   CALL GEMM3M(Af,Bf,MATMULVGc,'C','N')

!----------------------------------------
   END FUNCTION MATMULVGc
   END MODULE PcAB
