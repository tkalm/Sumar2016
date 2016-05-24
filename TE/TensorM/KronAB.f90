   MODULE KronAB
   CONTAINS
   FUNCTION KronM(Af,Bf)

! Kronecker product of two complex square matrices of the
! same size 

   USE Mod_Precision

   USE mkl95_lapack
   USE blas95
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: Nd1, Nd2, ierr, i, j, k, l
   COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: KronM, Af, Bf

   Nd1 = SIZE(Af,1)
   Nd2 = Nd1**2
   ierr = 0

 
   ALLOCATE(KronM(Nd2,Nd2), STAT=ierr)

   KronM = (0.0_dp, 0.0_dp)

!$OMP PARALLEL DO SHARED(Af, Bf, KronM) SCHEDULE(DYNAMIC) 
   DO j = 1, Nd1
     DO i = 1, Nd1
       DO l = 1, Nd1
         DO k = 1, Nd1  
           KronM(Nd1*(i-1)+k,Nd1*(j-1)+l) =  Af(i,j)*Bf(k,l)
         END DO
       END DO
     END DO
   END DO
!$OMP END PARALLEL DO
   

!----------------------------------------
   END FUNCTION KronM
   END MODULE KronAB
