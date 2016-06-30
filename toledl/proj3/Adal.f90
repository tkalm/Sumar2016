PROGRAM Adal
!----------------------------------------------------------------------
! This program finds the eigenvalues and vectors of the density matrix
! rho in the Liouville-von Neumann equation in Liouville space. 
!----------------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Precision_Mod         ! Module for setting double precision
   USE Init_Mod              ! Initial values
   USE Fields_Mod            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ----------------------------------------------
   INTEGER                           ::  i, j, ierr
!------- Output -------------------------------------------------------
   OPEN(UNIT=13,FILE=   'eigval.dtx'        ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'eigvec.dtx'        ,STATUS='NEW')
!------- Definitions --------------------------------------------------
   ierr = 0

! Unit matrix
   ALLOCATE(Imat(Nf,Nf), STAT=ierr)
   Imat = Czero
   DO j = 1, Nf
      DO i = 1 , Nf
         IF(j == i) Imat(i,j) = 1
      END DO
   END DO

! Ladder operators
   ALLOCATE(amat(Nf,Nf),admat(Nf,Nf), STAT=ierr)
   amat = Czero
   admat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
         IF(j == i+1)  admat(i,j)  = SQRT(FLOAT(i))                    ! Lowering operator
         IF(i == j+1)  amat(i,j) = SQRT(FLOAT(j))                      ! Raising operator
      END DO
   END DO
   Nmat       = matmul(admat,amat)                                     ! Defined for convenience 

! Hamiltonian
   ALLOCATE(Hmat(Nf,Nf), STAT=ierr)
   Hmat          = Czero
   DO i = 1, Nf
      Hmat(i,i)  = CMPLX((REAL(i-1,dp))+0.5_dp,0.0_dp,dp)
   END DO
   Hmat          = Hmat + Odo*(amat + admat)/SQRT(FLOAT(2))                  ! Add the external static electric field 

! Initial state of density matrix 
   ALLOCATE(rhomat(Nf,Nf), STAT=ierr)
   rhomat        = Czero
   rhomat(1,1)   = 1

! The Liouville operator 
   ALLOCATE(Limat(Nf2,Nf2), STAT=ierr)
   Limat      = kron(Imat,Hmat) - kron(TRANSPOSE(Hmat),Imat) + ci*0.5*kappa*&
     (2.0_dp*kron(amat,amat) - kron(Imat,Nmat) - kron(Nmat,Imat)) 
!------- Calculation & Output -----------------------------------------

   ALLOCATE(Eigval(Nf2),             STAT=ierr)
   ALLOCATE(vl(Nf2,Nf2),vr(Nf2,Nf2), STAT=ierr)
   CALL GEEV(Limat,Eigval,vl,vr)                                       ! Finds eigenvectors and values 

   DO i = 1, Nf2
      WRITE(13,FMT="(I6,1X,2(E15.8,1X))") i, Eigval(i)
   END DO

   DO j=1, Nf2
      DO i=1, Nf2
         WRITE(11,FMT="(E15.8,1X)") LOG(ABS(REAL(vl(i,1))))
      END DO
   END DO


   DEALLOCATE(vl, vr, Eigval, Limat, rhomat, Hmat, admat, amat, Imat)
!------ Functions used ------------------------------------------------
   CONTAINS
   FUNCTION kron(mat1,mat2)                                             ! Calculates the Kronecker product of matrices
      INTEGER                              :: r, s, v, w, m, n, p, q
      COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN)  :: mat1(:,:), mat2(:,:)
      COMPLEX(KIND=dp), DIMENSION(SIZE(mat1,1)*SIZE &
      (mat2,1),SIZE(mat1,2)*SIZE(mat2,2))  :: kron
      m = SIZE(mat1,1)
      n = SIZE(mat1,2)
      p = SIZE(mat2,1)
      q = SIZE(mat2,2)
      DO r=1,m
         DO s=1,n
            DO v=1,p
               DO w=1,q
                  kron(p*(r-1)+v,q*(s-1)+w) = mat1(r,s)*mat2(v,w)
               END DO
            END DO
         END DO
      END DO
   END FUNCTION kron
!----------------------------------------------------------------------
END PROGRAM Adal
