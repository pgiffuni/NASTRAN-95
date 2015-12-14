FUNCTION sadotb( a, b )
     
 
 REAL, INTENT(IN)                         :: a(3)
 REAL, INTENT(IN)                         :: b(3)
 
!*****
!  SINGLE-PRECISION VERSION
 
!  DOT PRODUCT A . B
!*****
 sadotb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
 RETURN
END FUNCTION sadotb
