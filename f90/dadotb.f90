DOUBLE PRECISION FUNCTION dadotb( a, b )
     
 
 DOUBLE PRECISION, INTENT(IN)             :: a(3)
 DOUBLE PRECISION, INTENT(IN)             :: b(3)
 
!*****
!  DOUBLE PRECISION VERSION
 
!  DOT PRODUCT  A . B
!*****
 dadotb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
 RETURN
END FUNCTION dadotb
