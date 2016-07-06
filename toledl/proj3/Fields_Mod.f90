MODULE Fields_Mod
   USE Precision_Mod
   USE Init_Mod
   IMPLICIT NONE
!-- Scalars -----------------------------------------------------------------------------------------
   INTEGER                     :: Nmax
   REAL(KIND=dp)               :: al 
   COMPLEX(KIND=dp)            :: cl
   CHARACTER(LEN=1)            :: TransA, TransB
!-- Vectors -----------------------------------------------------------------------------------------
   COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:)      :: Eigval, rho0v, rhotv
!-- Matrices ----------------------------------------------------------------------------------------
   COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:,:)    :: Imat, Hmat, amat, admat, &
                         Nmat, Limat, vl, vr, ssExpEnMat, rhoss, einnT, expiLt, vlU, vrV
!----------------------------------------------------------------------------------------------------
END MODULE Fields_Mod
