MODULE PcAB
   CONTAINS
   FUNCTION MATMULVGc(Af,Bf)
      USE Precision_Mod
      USE mkl95_lapack
      USE blas95
      USE omp_lib
      IMPLICIT NONE
!-------------------------------------------------------------------------
      INTEGER                                        ::  Nd1, Nd2, ierr
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  ::  MATMULVGc, Af, Bf
!-------------------------------------------------------------------------
      ierr = 0
      Nd1  = SIZE(Af,2)
      Nd2  = SIZE(Bf,2)
      ALLOCATE(MATMULVGc(Nd1,Nd2), STAT=ierr)
!-------------------------------------------------------------------------
      MATMULVGc = (0.0_dp, 0.0_dp)
      CALL GEMM3M(Af,Bf,MATMULVGc, "C", "N")
!-------------------------------------------------------------------------
   END FUNCTION MATMULVGc
END MODULE PcAB
