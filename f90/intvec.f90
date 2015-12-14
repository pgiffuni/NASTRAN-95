SUBROUTINE intvec (vector)
     
 
 INTEGER, INTENT(IN OUT)                  :: vector
 INTEGER :: xyzr(4),CHAR,vec(4),vecwrd
 COMMON /system/ skip(40), ncpw
 DATA    xyzr  / 1HX,1HY,1HZ,1HR /
 DATA    n     / 1HN/
 
 nshape = 0
 vecwrd = vector
 IF (vecwrd == 0) GO TO 125
 DO  i = 1,4
   vec(i) = 0
 END DO
 
!     SEPARATE THE FOUR CHARACTERS IN -VECWRD- (ANY COMBINATION OF THE
!     CHARACTERS X, Y, Z, AND R.
 
 DO  k = 1,4
   CHAR = klshft(vecwrd,(k-1))
   CHAR = krshft(CHAR,(ncpw-1))
   DO  i = 1,4
     IF (CHAR == krshft(xyzr(i),(ncpw-1))) GO TO 115
   END DO
   IF(CHAR == krshft(n,(ncpw-1))) nshape = 1
   CYCLE
   115 vec(i) = 1
 END DO
 
 vector = vec(1) + 2*vec(2) + 4*vec(3) + 8*vec(4)
 IF (vector == 8) vector = 15
 IF (nshape == 1) vector =-vector
 125 RETURN
END SUBROUTINE intvec
