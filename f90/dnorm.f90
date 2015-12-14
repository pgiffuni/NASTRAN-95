SUBROUTINE dnorm(x,mag)
     
!     DOUBLE PRECISION NORMALIZATION
 
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: x(3)
 DOUBLE PRECISION, INTENT(OUT)            :: mag
 DOUBLE PRECISION :: a
 
 mag= 0.d0
 a= x(1)*x(1) + x(2)*x(2) +x(3)*x(3)
 IF (a > 0.d0) mag= DSQRT(a)
 IF(mag == 0.0D0) RETURN
 x(1) = x(1) / mag
 x(2) = x(2) / mag
 x(3) = x(3) / mag
 RETURN
END SUBROUTINE dnorm
