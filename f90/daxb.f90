SUBROUTINE daxb(a,b,c)
     
 
 DOUBLE PRECISION, INTENT(IN)             :: a(3)
 DOUBLE PRECISION, INTENT(IN)             :: b(3)
 DOUBLE PRECISION, INTENT(OUT)            :: c(3)
 DOUBLE PRECISION :: d(3)
!*****
!  DOUBLE PRECISION VERSION
 
!  THIS ROUTINE PERFORMS A X B INTO C  (C MAY OVERLAP A OR B IN CORE)
!*****
 d(1) = a(2)*b(3) - a(3)*b(2)
 d(2) = a(3)*b(1) - a(1)*b(3)
 d(3) = a(1)*b(2) - a(2)*b(1)
 c(1) = d(1)
 c(2) = d(2)
 c(3) = d(3)
 RETURN
END SUBROUTINE daxb
