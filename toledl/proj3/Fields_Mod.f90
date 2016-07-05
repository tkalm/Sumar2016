MODULE Fields_Mod
!--------------------------------------------------------------------------
   USE Precision_Mod
   USE Init_Mod
   IMPLICIT NONE
!--------------------------------------------------------------------------
!  Scalars
   INTEGER                     :: Nmax
   REAL(KIND=dp)               :: al 
   COMPLEX(KIND=dp)            :: cl
   CHARACTER(LEN=1)            :: TransA, TransB
!--------------------------------------------------------------------------
!  Vectors
   COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:)      :: Eigval, minIm,&
                         minABS, rho0v, rhotv
!--------------------------------------------------------------------------
!  Matrices 
   COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:,:)    :: Hmat, Nrho, &
                         amat, admat, Nmat, Imat, Limat, vl, vr, einnT, &
                         rhoss, EaveM, expiLt, vlU, vrV 
!--------------------------------------------------------------------------
END MODULE Fields_Mod
