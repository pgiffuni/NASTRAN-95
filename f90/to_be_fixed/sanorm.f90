SUBROUTINE sanorm (*,a)
     
 
 , INTENT(IN OUT)                         :: *
 REAL, INTENT(IN OUT)                     :: a(3)
 
 
!     VECTOR NORMALIZATION AND VECTOR LENGTH
 
 xl=a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
 IF (xl <= 0.0) RETURN 1
 xl   = SQRT(xl)
 a(1) = a(1)/xl
 a(2) = a(2)/xl
 a(3) = a(3)/xl
 RETURN
END SUBROUTINE sanorm
