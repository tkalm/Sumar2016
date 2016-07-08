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
   REAL(KIND=dp)                     ::  ssExpEn, time
   COMPLEX(KIND=dp)                  ::  cnorm
!------- Output -------------------------------------------------------
   OPEN(UNIT=12,FILE=   'eigvallog.dtx'     ,STATUS='NEW')
   OPEN(UNIT=13,FILE=   'eigval.dtx'        ,STATUS='NEW')
   OPEN(UNIT=14,FILE=   'eigvec.dtx'        ,STATUS='NEW')
   OPEN(UNIT=15,FILE=   'nit.dtx'           ,STATUS='NEW')
   OPEN(UNIT=16,FILE=   'einnT.dtx'         ,STATUS='NEW')
   OPEN(UNIT=18,FILE=   'steadystate.dtx'   ,STATUS='NEW')
!------- Definitions --------------------------------------------------
   ierr = 0
! Unit matrix
   ALLOCATE(Imat(Nf,Nf), STAT=ierr)
   Imat = Czero
   DO j = 1, Nf
      Imat(j,j) = CUnit
   END DO

! Ladder operators
   ALLOCATE(amat(Nf,Nf),admat(Nf,Nf), Nmat(Nf,Nf), xmat(Nf,Nf), STAT=ierr)
   amat  = Czero
   admat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
         IF(j == i+1)  amat(i,j)   = SQRT(FLOAT(i))                             ! Lowering operator
         IF(i == j+1)  admat(i,j)  = SQRT(FLOAT(j))                             ! Raising  operator
      END DO
   END DO
   Nmat          = MATMULVG(admat,amat)                                         ! Number   operator 
   xmat          = (amat + admat)/sq2                                           ! Position operator

! Hamiltonian
   ALLOCATE(Hmat(Nf,Nf), STAT=ierr)
   Hmat          = Czero
   DO i = 1, Nf
      Hmat(i,i)  = CMPLX((REAL(i-1,dp))+0.5_dp,0.0_dp,dp)
   END DO
   Hmat          = Hmat + lambda*xmat

! The Liouville operator 
   ALLOCATE(Limat(Nf2,Nf2), STAT=ierr)
   Limat      = KronM(Imat,Hmat) - KronTM(Hmat,Imat) + &
   ci*kappa*0.5_dp*(2.0_dp*KronTM(admat,amat) - KronM(Imat,Nmat) - KronTM(Nmat,Imat)) 

!------- Finding the eigenvalues & -vectors -----------------------------------
   ALLOCATE(Eigval(Nf2), vl(Nf2,Nf2),vr(Nf2,Nf2), STAT=ierr)
   CALL GEEV(Limat,Eigval,vl,vr)                                                ! Find eigenvalues and left
                                                                                ! and right eigenvectors.
   DO i = 1, Nf2
      WRITE(13,FMT="(I6,1X,2(E15.8,1X))") i, Eigval(i)                          ! Print eigenvalues
      WRITE(12,FMT="(I6,1X,2(E15.8,1X))") i, LOG(ABS(REAL(Eigval(i)))), &       ! Print for log plot
                                             LOG(ABS(AIMAG(Eigval(i))))
   END DO

   DO j=1, Nf2
      DO i=1, Nf2
         WRITE(14,FMT="(I4,1X,I4,1X,1024(E15.8,1X))") i, j, vr(i,j)             ! Print r(ight)/l(eft) eigenvectors
      END DO
      WRITE(14,FMT="()")                                                        ! Blank line seperators between 
      WRITE(14,FMT="()")                                                        ! the vectors
   END DO

   DEALLOCATE(Limat, STAT=ierr)
!------- Finding the steady state ---------------------------------------------
   ALLOCATE(rhoss(Nf,Nf), STAT=ierr)

   ssloc = MINLOC(ABS(Eigval),1)
   DO j = 1, Nf
      DO i = 1, Nf
         l = Nf*(j-1) + i 
         rhoss(i,j) = vr(l,ssloc)                                               ! Rho steady state 
      END DO
   END DO
   
   ALLOCATE(ssExpEnMat(Nf,Nf), STAT=ierr)
   ssExpEnMat = MATMULVG(Hmat,rhoss)                                            ! Expectation of the energy of the steady state
   ssExpEn    = 0.0_dp
   DO i = 1, Nf
      ssExpEn = ssExpEn + REAL(ssExpEnMat(i,i),dp)
   END DO

   WRITE(18,*)                          '# Expectation of the energy of the steady state [meV]'
   WRITE(18,FMT='(1X,A1,1X,E15.8,1X)')  '#', ssExpEn
   WRITE(18,*)                          '#'
   WRITE(18,*)                          '# Location of the steady state'
   WRITE(18,FMT='(1X,A1,1X,I5)')        '#', ssloc
   WRITE(18,*)                          '#'
   WRITE(18,*)                          '# Occupation of harmonic oscillator states by the steady state'
   WRITE(18,*)                          '#  i  Occupation'
   DO i = 1, Nf
      WRITE(18,FMT='(I5,2X,E15.8,1X)')   i, REAL(rhoss(i,i),dp)
   END DO
   DEALLOCATE(rhoss, ssExpEnMat, STAT=ierr)

!------- Finding the time evolution -------------------------------------------
   ALLOCATE(expiLt(Nf2,Nf2), rho0v(Nf2), rhotv(Nf2), einnT(Nf2,Nf2), STAT=ierr)
   expiLt                =  Czero
   rho0v                 =  Czero
   einnT                 =  Czero
   rho0v(Nf*(mu0-1)+mu0) =  CUnit
   einnT                 =  MATMULVGc(vr,vl)
   ALLOCATE(vlU(Nf2,Nf2), vrV(Nf2,Nf2), STAT=ierr)

   DO j = 1, Nf2
      cnorm              =  1.0_dp/SQRT(einnT(j,j))
      DO i = 1, Nf2
         vlU(i,j)        =  vl(i,j)*cnorm
         vrV(i,j)        =  vr(i,j)*CONJG(cnorm)
      END DO
   END DO

   DEALLOCATE(vl, vr, STAT=ierr)

   einnT                 =  MATMULVGc(vrV,vlU)
   DO j = 1, Nf2
      DO i = 1, Nf2
         WRITE(16,FMT="(I5,1X,I5,1X,2(E15.8,1X))") i, j, (einnT(i,j))
      END DO
      WRITE(16,FMT="()")
   END DO

   DO j = 1, timetotal/delt
      time               =  REAL((j-1),dp)*delt
      DO i = 1, Nf2
         expiLt(i,i)     =  EXP(-ci*time*hbarinv*Eigval(i))
      END DO
      rhotv              =  MATMUL(MATMULVGz(MATMULVG(vrV,expiLt),vlU),rho0v) 
      DO i = 1, Nf
         WRITE(15,FMT="(1000(E15.8,1X))") time, (REAL(rhotv(Nf*(i-1)+i)))
      END DO
   END DO

!------------------------------------------------------------------------------
   DEALLOCATE(Eigval, expiLt, rho0v, rhotv, einnT,      STAT=ierr)
   DEALLOCATE(Hmat, admat, amat, Imat, Nmat, vlU, vrV,  STAT=ierr)
!------------------------------------------------------------------------------
END PROGRAM Adal
