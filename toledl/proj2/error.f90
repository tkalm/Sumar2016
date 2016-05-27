SUBROUTINE Trace_and_error(Trace,error)

USE Mod_Precision
USE Mod_Init
USE Mod_Fields
USE Mod_Matrices_Master_Equation
!-----------------------------------------
IMPLICIT NONE
!-----------------------------------------
INTEGER          :: mu,nu
COMPLEX(KIND=dp) :: error,Trace
!-----------------------------------------
Trace = (0.0_dp, 0.0_dp)
DO mu = 1, N_mes
Trace = Trace+rho1new(mu,mu)
END DO
error = (0.0_dp, 0.0_dp)
error = SQRT(ABS(SUM(rho1new-rho1)))
END SUBROUTINE Trace_and_error
