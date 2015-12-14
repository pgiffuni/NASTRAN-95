SUBROUTINE sub(x,y,a,b)
!*******
!     SUB WILL FORM Y = A*X - B*Y  WHERE A AND B ARE SCALAR MULTIPLIERS
!     FOR THE VECTORS X AND Y
!*******
 
 DOUBLE PRECISION, INTENT(IN)             :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: y(1)
 DOUBLE PRECISION, INTENT(IN)             :: a
 DOUBLE PRECISION, INTENT(IN)             :: b
 
 COMMON   /invpwx/  xx        ,ncol
 
 DO  i = 1,ncol
   y(i) = x(i)*a - y(i)*b
 END DO
 RETURN
END SUBROUTINE sub
