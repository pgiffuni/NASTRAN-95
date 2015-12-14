SUBROUTINE csub (x,y,z,a,b)
!*******
!     CSUB WILL FORM Z = A*X - B*Y WHERE A AND B ARE SCALAR
!     MULTIPLIERS FOR THE COMPLEX VECTORS X AND Y
!*******
 
 DOUBLE PRECISION, INTENT(IN)             :: x(2)
 DOUBLE PRECISION, INTENT(IN)             :: y(2)
 DOUBLE PRECISION, INTENT(OUT)            :: z(1)
 DOUBLE PRECISION, INTENT(IN)             :: a(2)
 DOUBLE PRECISION, INTENT(IN)             :: b(2)
 DOUBLE PRECISION :: , dum
 COMMON   /cinvpx/  aaa       ,ncol
 
 ncol2 = ncol+ncol
 DO  i = 1,ncol2,2
   dum  = x(i)*a(1) - x(i+1)*a(2) - y(i)*b(1) + y(i+1)*b(2)
   z(i+1) = x(i)*a(2) + x(i+1)*a(1) - y(i+1)*b(1) - y(i)*b(2)
   z(i) = dum
 END DO
 RETURN
END SUBROUTINE csub
