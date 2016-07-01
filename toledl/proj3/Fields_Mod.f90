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
                         minABS
!--------------------------------------------------------------------------
!  Matrices 
   COMPLEX(KIND=dp),     ALLOCATABLE, DIMENSION(:,:)    :: Hmat, rhomat, &
                         amat, admat, Nmat, Imat, Limat, Eigvect, vl, vr,&
                         rhoss
!--------------------------------------------------------------------------
END MODULE Fields_Mod
