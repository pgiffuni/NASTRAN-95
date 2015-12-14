DOUBLE PRECISION FUNCTION dvmag (v1,eps)
     
!     RETURNS DOUBLE PRECISION MAGNITUDE OF VECTOR V1
!        DVMAG= 0.D0 WHEN .LE. EPS
 
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: v1(3)
 DOUBLE PRECISION, INTENT(IN OUT)         :: eps
 DOUBLE PRECISION :: a, dadotb
 
 
 dvmag= 0.d0
 a= dadotb(v1,v1)
 IF (a > 0.d0) dvmag= DSQRT(a)
 IF ( dvmag <= eps)             dvmag= 0.d0
 RETURN
END FUNCTION dvmag
