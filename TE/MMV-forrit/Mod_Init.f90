   MODULE Mod_Init

!-----------------------------------
! Global parameters, parameters can not be changed in a program
!-----------------------------------

    USE Mod_Precision

    IMPLICIT NONE

    INTEGER,          PARAMETER     :: NumThreads     = 2 !4

    INTEGER,          PARAMETER     :: Nf  = 3 !128    ! 128


!----------------------------------------------------------------------------

    REAL(KIND=dp),     PARAMETER    :: pi   = 3.14159265358979324_dp
    REAL(KIND=dp),     PARAMETER    :: pi2i = 1.0_dp/(2.0_dp*pi)
    REAL(KIND=dp),     PARAMETER    :: pid2 = pi/2.0_dp

    COMPLEX(KIND=dp),  PARAMETER    :: ci    = CMPLX(0.0_dp, 1.0_dp)
    COMPLEX(KIND=dp),  PARAMETER    :: CUnit = CMPLX(1.0_dp, 0.0_dp)
    COMPLEX(KIND=dp),  PARAMETER    :: Czero = CMPLX(0.0_dp, 0.0_dp)

    CHARACTER(LEN=1), PARAMETER     :: JOB   = 'V', UPLO  = 'U'
    CHARACTER(LEN=1), PARAMETER     :: TRANS = 'N', RANGO = 'I'


!---------------------------------------------------------------------------
   END MODULE Mod_Init
