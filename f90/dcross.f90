SUBROUTINE dcross(x,y,z)
     
!     DOUBLE PRECISION CROSS PRODUCT
 
 
 
 DOUBLE PRECISION, INTENT(IN)             :: x(3)
 DOUBLE PRECISION, INTENT(IN)             :: y(3)
 DOUBLE PRECISION, INTENT(OUT)            :: z(3)
 
 
 z(1) = x(2)*y(3) - x(3)*y(2)
 z(2) = y(1)*x(3) - y(3)*x(1)
 z(3) = x(1)*y(2) - x(2)*y(1)
 RETURN
END SUBROUTINE dcross
