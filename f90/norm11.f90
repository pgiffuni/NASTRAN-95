SUBROUTINE norm11(x,div)
     
 
 REAL, INTENT(IN OUT)                     :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: div
 
 REAL :: MAX
 COMMON   /invpwx/  filek(7)
 EQUIVALENCE        (ncol,filek(2))
 DATA ind1 /1/
 
 MAX = 0.0
 DO  i=1,ncol
   xx = ABS( x(i) )
   IF( xx <= MAX ) CYCLE
   MAX = xx
   ind = i
 END DO
 IF( x(ind) < 0.0 ) ind = -ind
 i = IABS(ind1)
 xx = x(i)
 div = SIGN(1.,xx)*FLOAT(ISIGN(1,ind1))*MAX
 xx = div
 ind1 = ind*IFIX(SIGN(1.,xx))
 MAX = 1.0 /div
 DO  i=1,ncol
   xi = x(i)*MAX
   IF (ABS(xi) < 1.e-36) xi = 0.0
   x(i)= xi
 END DO
 RETURN
END SUBROUTINE norm11
