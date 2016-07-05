MODULE PABz
   CONTAINS
   FUNCTION MATMULVGz(Af,Bf)
      USE Precision_Mod
      USE mkl95_lapack
      USE blas95
      USE omp_lib
      IMPLICIT NONE
!-------------------------------------------------------------------------
      INTEGER                                        ::  Nd1, Nd2, ierr
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  ::  MATMULVGz, Af, Bf
!-------------------------------------------------------------------------
      ierr = 0
      Nd1  = SIZE(Af,1)
      Nd2  = SIZE(Bf,1)
      ALLOCATE(MATMULVGz(Nd1,Nd2), STAT=ierr)
!-------------------------------------------------------------------------
      MATMULVGz = (0.0_dp, 0.0_dp)
      CALL GEMM3M(Af,Bf,MATMULVGz,"N","C")
!-------------------------------------------------------------------------
   END FUNCTION
END MODULE PABz
