SUBROUTINE sub1(x,y,a,b)
!     SUBROUTINE SUB(X,Y,A,B)
!*******
!     SUB WILL FORM Y = A*X - B*Y  WHERE A AND B ARE SCALAR MULTIPLIERS
!     FOR THE VECTORS X AND Y
!*******
!     DOUBLE PRECISION   X(1)      ,Y(1)     ,A        ,B
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(OUT)                        :: y(1)
 DOUBLE PRECISION, INTENT(IN)             :: a
 DOUBLE PRECISION, INTENT(IN)             :: b
 
 
 COMMON   /invpwx/  xx        ,ncol
 
 a1 = a
 b1 = b
 DO  i = 1,ncol
   y(i) = x(i)*a1- y(i)*b1
 END DO
 RETURN
END SUBROUTINE sub1
