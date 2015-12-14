SUBROUTINE cdivid(a,b,d,ncol)
     
 
 DOUBLE PRECISION, INTENT(IN)             :: a(1)
 DOUBLE PRECISION, INTENT(OUT)            :: b(1)
 DOUBLE PRECISION, INTENT(IN)             :: d(2)
 INTEGER, INTENT(IN)                      :: ncol
 DOUBLE PRECISION :: dtemp,denm
 
!     THIS ROUTINE DIVIDES THE VECTOR A BY D AND STORE RESULT IN B
 
 denm = d(1)**2 + d(2)**2
 DO  i = 1,ncol,2
   dtemp = (a(i)*d(1) +a(i+1)*d(2))/denm
   b(i+1) = (a(i+1)*d(1) -a(i) * d(2))/denm
   b(i) = dtemp
 END DO
 RETURN
END SUBROUTINE cdivid
