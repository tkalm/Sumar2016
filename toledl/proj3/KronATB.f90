MODULE KronATB
   CONTAINS
   FUNCTION KronTM(mat1,mat2)
      USE Precision_Mod
      USE mkl95_lapack
      USE blas95
      USE omp_lib
      IMPLICIT NONE
!-------------------------------------------------------------------------
      INTEGER                                        ::  r, s, v, w, m, n, &
                                                         p, q, ierr
      COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN)   ::  mat1(:,:), mat2(:,:)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  ::  mat1t(:,:)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:)  ::  KronTM 
!----------------------------------------------------------------------------
      ierr = 0
      n     = SIZE(mat1,1)
      m     = SIZE(mat1,2)
      p     = SIZE(mat2,1)
      q     = SIZE(mat2,2)
      ALLOCATE(mat1t(m,n), STAT=ierr) 
      mat1t = TRANSPOSE(mat1)
      ALLOCATE(KronTM(m*p,n*q), STAT=ierr)
!----------------------------------------------------------------------------
      DO r=1,m
         DO s=1,n
            DO v=1,p
               DO w=1,q
                  KronTM(p*(r-1)+v,q*(s-1)+w) = mat1t(r,s)*mat2(v,w)
               END DO
            END DO
         END DO
      END DO
!----------------------------------------------------------------------------
      DEALLOCATE(mat1t, STAT=ierr) 
   END FUNCTION KronTM
END MODULE KronATB 
