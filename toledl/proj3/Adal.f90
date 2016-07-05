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

   USE KronAB                ! Function KronM:  Kronecker product
   USE KronATB               ! Function KronTM: First matrix transposed
   USE PAB                   ! Function MATMULVG: Product of matrices
   USE PcAB                  ! Function MATMULVGc: First matrix transp 
   USE PABz                  ! Function MATMULVGz: Latter matrix transp

   IMPLICIT NONE
!------- Local variables ----------------------------------------------
   INTEGER                           ::  i, j, l, ierr, ssloc
   REAL(KIND=dp)                     ::  Nsum, Eave, time
   COMPLEX(KIND=dp)                  ::  cnorm
!------- Output -------------------------------------------------------
   OPEN(UNIT=12,FILE=   'test.dtx'          ,STATUS='NEW')
   OPEN(UNIT=13,FILE=   'eigval.dtx'        ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'eigvec.dtx'        ,STATUS='NEW')
   OPEN(UNIT=15,FILE=   'nit.dtx'           ,STATUS='NEW')
   OPEN(UNIT=16,FILE=   'einnT.dtx'         ,STATUS='NEW')
   OPEN(UNIT=18,FILE=   'NEave.dtx'         ,STATUS='NEW')
!------- Definitions --------------------------------------------------
   ierr = 0
! Unit matrix
   ALLOCATE(Imat(Nf,Nf), STAT=ierr)
   Imat = Czero
   DO j = 1, Nf
      Imat(j,j) = CUnit
   END DO

! Ladder operators
   ALLOCATE(amat(Nf,Nf),admat(Nf,Nf), Nmat(Nf,Nf), STAT=ierr)
   amat = Czero
   admat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
         IF(j == i+1)  admat(i,j) = SQRT(FLOAT(i))                     ! Lowering operator
         IF(i == j+1)  amat(i,j)  = SQRT(FLOAT(j))                     ! Raising operator
      END DO
   END DO
   Nmat          = MATMULVG(admat,amat)                                ! Defined for convenience 

! Hamiltonian
   ALLOCATE(Hmat(Nf,Nf), STAT=ierr)
   Hmat          = Czero
   DO i = 1, Nf
      Hmat(i,i)  = CMPLX((REAL(i-1,dp))+0.5_dp,0.0_dp,dp)
   END DO
   Hmat          = Hmat + Odo*(amat + admat)/SQRT(FLOAT(2))                  ! Add the external static electric field 

! The Liouville operator 
   ALLOCATE(Limat(Nf2,Nf2), STAT=ierr)
   Limat      = KronM(Imat,Hmat) - KronTM(Hmat,Imat) - ci*kappa*KronTM(admat,amat)
   Limat = Limat - ci*0.5*kappa*(KronM(Imat,Nmat) - KronTM(Nmat,Imat)) 

!------------------------------------------------------------------------------
   ALLOCATE(Eigval(Nf2),             STAT=ierr)
   ALLOCATE(vl(Nf2,Nf2),vr(Nf2,Nf2), STAT=ierr)
   CALL GEEV(Limat,Eigval,vl,vr)                                       ! Finds eigenvectors and values 

   DO i = 1, Nf2
      WRITE(13,FMT="(I6,1X,2(E15.8,1X))") i, Eigval(i)
   END DO

   DO j=1, Nf2
      DO i=1, Nf2
         WRITE(14,FMT="(I4,1X,I4,1X,1024(E15.8,1X))") i, j, vr(i,j)
      END DO
      WRITE(14,FMT="()")
      WRITE(14,FMT="()")
   END DO

   DEALLOCATE(Limat, STAT=ierr)

!------------------------------------------------------------------------------
   ALLOCATE(rhoss(Nf,Nf), minIm(Nf2), minABS(Nf2), STAT=ierr)
   WRITE(12,FMT="(I6)") 1 
   minIm  = 100
   minABS = 100
   minIm  = MINLOC(ABS(AIMAG(Eigval)))
   minABS = MINLOC(ABS(      Eigval))
   IF(minIm(1) == minABS(1)) ssloc = minABS(1)                         ! Steady state location
   DO j = 1, Nf
      DO i = 1, Nf
         l = Nf*(j-1) + i 
         rhoss(i,j) = vr(l,ssloc)                                      ! Rho steady state 
      END DO
   END DO
   
   ALLOCATE(Nrho(Nf,Nf),EaveM(Nf,Nf), STAT=ierr)                       ! Energy average matrix
   Nrho  = MATMULVG(Nmat,rhoss)
   EaveM = MATMULVG(Hmat,rhoss)
   Nsum  = 0.0_dp
   Eave  = 0.0_dp
   DO i = 1, Nf
      Nsum = Nsum + REAL(Nrho(i,i), dp)
      Eave = Eave + REAL(EaveM(i,i),dp)
   END DO

   WRITE(18,*)                          '# Nsum, Eave'
   WRITE(18,FMT='(A1,1X,2(E15.8,1X))')  '#', Nsum, Eave
   WRITE(18,*)                          '#'
   WRITE(18,*)                          '# ssloc'
   WRITE(18,FMT='(A1,1X,I5)')           '#', ssloc
   WRITE(18,*)                          '# i occ'
   DO i = 1, Nf
      WRITE(18,FMT='(I5,1X,E15.8,1X)')   i, REAL(rhoss(i,i),dp)
   END DO
   DEALLOCATE(minIm, minABS, STAT=ierr)
   DEALLOCATE(rhoss, Nrho, EaveM, STAT=ierr)

!------------------------------------------------------------------------------
   ALLOCATE(expiLt(Nf2,Nf2), rho0v(Nf2), rhotv(Nf2), einnT(Nf2,Nf2), STAT=ierr)
   expiLt               =  Czero
   rho0v                =  Czero
   einnT                =  Czero
   rho0v(Nf*(mu0-1)+mu0) =  CUnit
   einnT                =  MATMULVGc(vr,vl)
   ALLOCATE(vlU(Nf2,Nf2), vrV(Nf2,Nf2), STAT=ierr)
   DO j = 1, Nf2
      cnorm             =  1.0_dp/SQRT(einnT(j,j))
      DO i = 1, Nf2
         vlU(i,j)       =  vl(i,j)*cnorm
         vrV(i,j)       =  vr(i,j)*CONJG(cnorm)
      END DO
   END DO
   DEALLOCATE(vl, vr, STAT=ierr)
   einnT                =  MATMULVGc(vrV,vlU)
   DO j = 1, Nf2
      DO i = 1, Nf2
         WRITE(16,FMT="(I5,1X,I5,1X,2(E15.8,1X))") i, j, (einnT(i,j))
      END DO
      WRITE(16,FMT="()")
   END DO

   DO j = 1, timeps
      time              =  REAL((j-1),dp)*theta
      DO i = 1, Nf2
         expiLt(i,i)    =  EXP(-ci*time*Eigval(i))
      END DO
      rhotv             =  MATMUL(MATMULVGz(MATMULVG(vrV,expiLt),vlU),rho0v) 
      DO i = 1, Nf
         WRITE(15,FMT="()") time*0.65820_dp, (REAL(rhotv(Nf*(i-1)+i)))
      END DO
   END DO

!------------------------------------------------------------------------------
   DEALLOCATE(Eigval,                                   STAT=ierr)
   DEALLOCATE(expiLt, rho0v, rhotv, einnT,              STAT=ierr)
   DEALLOCATE(Hmat, admat, amat, Imat, Nmat, vlU, vrV,  STAT=ierr)
!------------------------------------------------------------------------------
END PROGRAM Adal
