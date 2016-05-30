   MODULE Mod_Fields

!-----------------------------------

    USE Mod_Precision
    USE Mod_Init

    IMPLICIT NONE

    INTEGER                     :: Nmax

    REAL(KIND=dp)               :: al 

    COMPLEX(KIND=dp)            :: cl

    CHARACTER(LEN=1)            :: TransA, TransB



!----------------------------------------------------------------

    REAL(KIND=dp),        ALLOCATABLE, DIMENSION(:)      :: Eigval
    
!    COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:,:)    :: H0, rho0
    
!    COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:,:,:)    :: rho1
   

!--------------------------------------------------------------------------

   END MODULE Mod_Fields
