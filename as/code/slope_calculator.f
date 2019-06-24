! source file: /net/mare/home1/eby/as/ism/slope_calculator.F
      subroutine slope_calculator(c)

      return
      end
!          ****************************************************************************
      FUNCTION ANGVEC (V1, V2) RESULT (COSTHETA)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: V1, V2
      DOUBLE PRECISION :: COSTHETA
      COSTHETA = DOT_PRODUCT(V1,V2)/
     &(SQRT(DOT_PRODUCT(V1,V1))*SQRT(DOT_PRODUCT(V2,V2)))
      RETURN
      END FUNCTION ANGVEC

