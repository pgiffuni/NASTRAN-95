SUBROUTINE norm1 (x,div)
!*******
!     NORM WILL NORMALIZE X TO MAXIMUM ELEMENT EQUAL TO ONE AND STORE TH
!     DIVISOR IN MAX
!*******
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: div
 DOUBLE PRECISION :: MAX
 COMMON   /invpwx/  filek(7)
 EQUIVALENCE        (ncol,filek(2))
 DATA ind1 /1/
 
 MAX = 0.d0
 DO  i=1,ncol
   div = DABS( x(i) )
   IF( div <= MAX ) CYCLE
   MAX = div
   ind = i
 END DO
 IF( x(ind) < 0.d0 ) ind = -ind
 i = IABS(ind1)
 xx = x(i)
 div = SIGN(1.,xx)*FLOAT(ISIGN(1,ind1))*MAX
 xx = div
 ind1 = ind*IFIX(SIGN(1.,xx))
 MAX = 1.d0/div
 DO  i=1,ncol
   xi = x(i)*MAX
   IF (ABS(xi) < 1.e-36) xi=0.
   x(i)= xi
 END DO
 RETURN
END SUBROUTINE norm1
