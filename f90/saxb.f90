SUBROUTINE saxb( a, b, c )
     
 
 REAL, INTENT(IN)                         :: a(3)
 REAL, INTENT(IN)                         :: b(3)
 REAL, INTENT(OUT)                        :: c(3)
 REAL :: d(3)
!*****
!  SINGLE-PRECISION VERSION
 
!  THIS ROUTINE PERFORMS A X B INTO C.  (C MAY OVERLAP A OR B IN CORE.)
!*****
 d(1) = a(2)*b(3) - a(3)*b(2)
 d(2) = a(3)*b(1) - a(1)*b(3)
 d(3) = a(1)*b(2) - a(2)*b(1)
 c(1) = d(1)
 c(2) = d(2)
 c(3) = d(3)
 RETURN
 
!*****
 ENTRY sapb( a, b, c )
 
!  THIS ROUTINE PERFORMS A + B INTO C.
!*****
 c(1) = a(1) + b(1)
 c(2) = a(2) + b(2)
 c(3) = a(3) + b(3)
 RETURN
 
!*****
 ENTRY samb( a, b, c )
 
!  THIS ROUTINE PERFORMS A - B INTO C.
!*****
 c(1) = a(1) - b(1)
 c(2) = a(2) - b(2)
 c(3) = a(3) - b(3)
 RETURN
END SUBROUTINE saxb
