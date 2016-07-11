MODULE Init_Mod
!---------------------------------------------------------------------
    USE Precision_Mod
    IMPLICIT NONE
!---------------------------------------------------------------------
! Calculation constants
    INTEGER,          PARAMETER     :: NumThreads = 2
    INTEGER,          PARAMETER     :: Nf      = 16
    INTEGER,          PARAMETER     :: Nf2     = Nf*Nf
    INTEGER,          PARAMETER     :: mu0     = 1                       ! Initial state of rho
! 
    REAL(KIND=dp),     PARAMETER    :: lambda    = 0.4_dp                ! Omega (capital) divided by omega, the ratio of ang. freq.
    REAL(KIND=dp),     PARAMETER    :: hbom      = 1.0_dp                ! h_bar*omega in units of [meV] 
    REAL(KIND=dp),     PARAMETER    :: kappa     = 0.2_dp*0.5_dp         ! The strength of dissipation 
    REAL(KIND=dp),     PARAMETER    :: hbar      = 6.582119514E-4        ! In units of eV*ps 
    REAL(KIND=dp),     PARAMETER    :: hbarinv   = 1.5193_dp             ! 1/[h_bar] in units of 1/[ps*meV]
    REAL(KIND=dp),     PARAMETER    :: timetotal = 100                  ! Length of simulation (in ps)  
    REAL(KIND=dp),     PARAMETER    :: delt      = 0.1                   ! Size of timestep (in ps)
! Mathematical constants
    REAL(KIND=dp),     PARAMETER    :: pi      = 3.14159265358979324_dp
    REAL(KIND=dp),     PARAMETER    :: pi2i    = 1.0_dp/(2.0_dp*pi)
    REAL(KIND=dp),     PARAMETER    :: pid2    = pi/2.0_dp
    REAL(KIND=dp),     PARAMETER    :: sq2     = SQRT(REAL(2.0,dp))

    COMPLEX(KIND=dp),  PARAMETER    :: ci      = CMPLX(0.0_dp, 1.0_dp)
    COMPLEX(KIND=dp),  PARAMETER    :: CUnit   = CMPLX(1.0_dp, 0.0_dp)
    COMPLEX(KIND=dp),  PARAMETER    :: Czero   = CMPLX(0.0_dp, 0.0_dp)

    CHARACTER(LEN=1), PARAMETER     :: JOB     = 'V', UPLO  = 'U'
    CHARACTER(LEN=1), PARAMETER     :: TRANS   = 'N', RANGO = 'I'

!---------------------------------------------------------------------
END MODULE Init_Mod
