PROGRAM Adal
!----------------------------------------------------------------------
! This program finds the eigenvalues and vectors of the density matrix
! rho in the Liouville-von Neumann equation in Liouville space. 
!----------------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ----------------------------------------------
   INTEGER                           ::  i, j 
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, rho, amat, admat, Nmat, Imat, CapA
   COMPLEX, DIMENSION(Nf*Nf)         ::  vecrho
   REAL, DIMENSION(Nf*Nf)            ::  Eigenval
   COMPLEX, DIMENSION(Nf*Nf,Nf*Nf)   ::  CompMat, Eigenvect
   REAL                              ::  hbaromega, kappa, Odo
!------- Output -------------------------------------------------------
   OPEN(UNIT=10,FILE=   'eigval.dtx'       ,STATUS='NEW')
   OPEN(UNIT=11,FILE=   'reigvec.dtx'       ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'ieigvec.dtx'       ,STATUS='NEW')
!------- Definitions --------------------------------------------------
! Define the unit matrix
   Imat = Czero
   DO j = 1, Nf
      DO i = 1 , Nf
         IF(j == i) Imat(i,j) = 1
      END DO
   END DO

! Define hamiltonian of system and ladder operators used in lambda 
   H = Czero
   amat = Czero
   admat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
         IF(j == i+1)  admat(i,j)  = SQRT(FLOAT(i))                    ! Lowering operator
         IF(i == j+1)  amat(i,j) = SQRT(FLOAT(j))                      ! Raising operator
         IF(i == j)    H(i,j)     = i-0.5                              ! Hamiltonian of harmonic oscillator
      END DO
   END DO
   Odo        = 0.4                                                    ! Omega (capital) divided by omega
   H          = H + Odo*(amat + admat)/SQRT(FLOAT(2))                  ! Add the external static electric field 
   Nmat       = matmul(admat,amat)                                     ! Defined for convenience 

! Initial state of rho   
   rho        = Czero
   rho(1,1)   = 1
   vecrho     = vec(rho)

! Define constants used 
   hbaromega  = 1E-3                                                   ! Our energy scale, hbar*omega in eV
   kappa      = 1E-1                                                   ! The strength of dissipation

! Define A and the complete matrix in the linear system
   CapA = H - ci*kappa*0.5*Nmat
   CompMat = kron(Imat,CapA) + kron(TRANSPOSE(CapA),Imat) + ci*kappa*kron(amat,amat)

!------- Calculation --------------------------------------------------

   CALL HEEVR(CompMat,Eigenval,UPLO,Eigenvect)   ! Subroutine from MKL for finding

   DO j=1,Nf*Nf
      WRITE(10,FMT="(E15.8)") Eigenval(j)
         DO i=1,Nf*Nf
            WRITE(11,FMT="(E15.8)") LOG(ABS(REAL(Eigenvect(i,j))))
            WRITE(12,FMT="(E15.8,2X,E15.8)") CompMat(i,j)
         END DO
   END DO 
!------ Functions used ------------------------------------------------
   CONTAINS
   FUNCTION vec(mat)
      INTEGER                                     :: k, l
      COMPLEX, DIMENSION(:,:), INTENT(IN)         :: mat(:,:)
      COMPLEX, DIMENSION(SIZE(mat,1)*SIZE(mat,2)) :: vec
      DO k=1,SIZE(mat,1)
         DO l=1,SIZE(mat,2)
            vec((k-1)*SIZE(mat,2)+l) = mat(k,l)
         END DO
      END DO
   END FUNCTION vec

   FUNCTION kron(mat1,mat2)
      INTEGER                                      :: r, s, v, w, m, n, p, q
      COMPLEX, DIMENSION(:,:), INTENT(IN)          :: mat1(:,:), mat2(:,:)
      COMPLEX, DIMENSION(SIZE(mat1,1)*SIZE(mat2,1),SIZE(mat1,2)*SIZE(mat2,2))  :: kron
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

   FUNCTION tr(matrix)                                                 ! Calculates the trace of a matrix
      INTEGER                             :: m
      COMPLEX                             :: tr
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: matrix 
      tr = Czero 
      DO m = 1, Nf
         tr = tr+matrix(m,m)
      END DO
   END FUNCTION tr
!----------------------------------------------------------------------
END PROGRAM Adal
