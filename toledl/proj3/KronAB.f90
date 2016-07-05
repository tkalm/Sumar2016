MODULE KronAB
   CONTAINS
   FUNCTION KronM(mat1,mat2)
      USE Precision_Mod
      USE mkl95_lapack
      USE blas95
      USE omp_lib
      IMPLICIT NONE
!-------------------------------------------------------------------------
      INTEGER                                        ::  r, s, v, w, m, n, &
                                                         p, q, ierr
      COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN)   ::  mat1(:,:), mat2(:,:)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  ::  KronM 
!----------------------------------------------------------------------------
      ierr = 0
      m    = SIZE(mat1,1)
      n    = SIZE(mat1,2)
      p    = SIZE(mat2,1)
      q    = SIZE(mat2,2)
      ALLOCATE(KronM(m*p,n*q), STAT=ierr)
!----------------------------------------------------------------------------
      DO r=1,m
         DO s=1,n
            DO v=1,p
               DO w=1,q
                  KronM(p*(r-1)+v,q*(s-1)+w) = mat1(r,s)*mat2(v,w)
               END DO
            END DO
         END DO
      END DO
!----------------------------------------------------------------------------
   END FUNCTION KronM
END MODULE KronAB 
