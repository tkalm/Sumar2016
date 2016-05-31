FUNCTION lambda(rho) 
!------------------------------------------------------------
! This function takes in a matrix and calculates its 
! commutator with the hamiltonian of the harmonic oscillator
!------------------------------------------------------------
   USE omp_lib               ! For OpenMP parallel processing
   USE Mod_Precision         ! Module for setting double precision
   USE Mod_Init              ! Initial values
   USE Mod_Fields            ! Global variables 

   USE mkl95_lapack          ! Subroutines
   USE mkl95_blas            ! Subroutines

   IMPLICIT NONE
!------- Local variables ------------------------------------
   INTEGER         :: t
   COMPLEX         :: rho(:,:), comm(Nf,Nf)
!------------------------------------------------------------
   ierr = 0
!------------------------------------------------------------

   comm = matmul(H0,rho)-matmul(rho,H0)
   lambda = -ci*comm 

!-----------------------------------------------------------
END FUNCTION 
