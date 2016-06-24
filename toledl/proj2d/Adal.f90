   PROGRAM Adal
!----------------------------------------------------------------------
! This program does the same thing as project 2b (calculates occupation
! of states of the system H=H0+H'(t) for t>0 given the state at t=0) 
! except this time dissipation is added to the model (so if kappa=0 
! this will give the same result).
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
   COMPLEX, DIMENSION(Nf,Nf)         ::  H, rho, rhon1, rhon2, amat, admat, Nmat, xmat, xxmat
   REAL                              ::  hbaromega, delt, alpha, kappa, Odo
!------- Output -------------------------------------------------------
   OPEN(UNIT=11,FILE=   'expx.dtx'       ,STATUS='NEW')
   OPEN(UNIT=12,FILE=   'expxx.dtx'      ,STATUS='NEW')
!------- Definitions --------------------------------------------------
! Define hamiltonian of system and ladder operators used in lambda 
   H = Czero
   amat = Czero
   admat = Czero
   DO j = 1, Nf
      DO i = 1, Nf
         IF(j == i+1)  amat(i,j)  = SQRT(FLOAT(i))                     ! Lowering operator
         IF(i == j+1)  admat(i,j) = SQRT(FLOAT(j))                     ! Raising operator
         IF(i == j)    H(i,j)     = i-0.5                              ! Hamiltonian of harmonic oscillator
      END DO
   END DO
   xmat       = (amat + admat)/SQRT(FLOAT(2))
   xxmat      = matmul(xmat,xmat)
   Odo        = 0.4                                                    ! Omega (capital) divided by omega
   H          = H + Odo*xmat                                           ! Add the external static electric field 
   Nmat       = matmul(admat,amat)                                     ! Defined for convenience (used in lambda)

! Initial state of rho   
   rho        = Czero
   rho(2,2)   = 1
   rhon1      = rho                                                    ! Setup for the iteration

! Print the initial expectation values 
   WRITE(11,FMT="(E15.8,2X,E15.8)") 0, REAL(tr(matmul(rho,xmat)))      ! The trace of rho times x at time t=0
   WRITE(12,FMT="(E15.8,2X,E15.8)") 0, REAL(tr(matmul(rho,xxmat)))     ! The trace of rho times x^2 at time t=0

! Define constants used 
   hbaromega  = 1E-3                                                   ! Our energy scale, hbar*omega in eV
   delt       = 1E-2                                                   ! The timestep of our approximation in ps 
   alpha      = hbaromega*delt/(2*hbar)                                ! The constant in our equation below (hbar is in eV*ps)
   kappa      = (5E-2)/2                                               ! The strength of dissipation (kappa/2) (appears in lambda)
!------- Calculation --------------------------------------------------
! Calculate rho for t>0 by iteration of the Crank-Nicolson 
! approximation of the of  L-vN equation.  
   DO j=1, Nt/delt                                                     ! Time grid
      DO                                                               ! Iterations
        rhon2 = rho + alpha*(lambda(rho) + lambda(rhon1))              ! Our approximation
        IF(   err(rhon1,rhon2) < 1E-3   )  EXIT                        ! Exit loop when desired accuracy is achieved
        rhon1 = rhon2                                                  ! Setup for next iteration
      END DO
      rho     = rhon2                                                  ! Setup for next timestep
      rhon1   = rho                                                    ! Setup for next iteration

! Print the expectation values at t = j*delt. 
      WRITE(11,FMT="(E15.8,2X,E15.8)") j*delt, REAL(tr(matmul(rho,xmat)))      ! The trace of rho times x at time t=0
      WRITE(12,FMT="(E15.8,2X,E15.8)") j*delt, REAL(tr(matmul(rho,xxmat)))     ! The trace of rho times x^2 at time t=0
   END DO
!------ Functions used ------------------------------------------------
   CONTAINS
   FUNCTION lambda(mat)                                                ! The superoperator in our equation 
      COMPLEX, DIMENSION(Nf,Nf)                :: lambda, matad
      COMPLEX, DIMENSION(:,:), INTENT(IN)      :: mat(:,:)
      lambda = -ci*(matmul(H,mat)-matmul(mat,H))+ &                    ! Commutator of hamiltonian and rho 
      kappa*(2*matmul(amat,matmul(mat,admat)) - &                      ! The dissipation terms
      matmul(Nmat,mat) - matmul(mat,Nmat))                             ! The dissipation terms
   END FUNCTION lambda

   FUNCTION err(x1,x2)                                                 ! Calculates the error
      REAL                                :: err 
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: x1,x2
      err = SQRT(ABS(SUM(x2-x1)))
   END FUNCTION err

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
