SUBROUTINE rotate (da,row,row1,row2, o,SIN,COS)
     
!     THIS ROUTINE IS CALLED ONLY BY TRIDI SUBROUTINE, WHICH IS CALLED
!     ONLY BY VALVEC
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: da(1)
 INTEGER, INTENT(IN OUT)                  :: row
 INTEGER, INTENT(IN)                      :: row1
 INTEGER, INTENT(IN)                      :: row2
 DOUBLE PRECISION, INTENT(IN OUT)         :: o(1)
 DOUBLE PRECISION, INTENT(IN)             :: SIN(1)
 DOUBLE PRECISION, INTENT(IN)             :: COS(1)
 
!    1,                CHECK
 DOUBLE PRECISION :: sine,cosine,x,y,z
 COMMON /givn  /  title(100),n
 
!     O     = 2ND ROW OF THE COMPLETE MATRIX.
!     SIN   = SINES.
!     COS   = COSINES.
!     DA = MATRIX PARTITION (TRIANGULAR) - DOUBLE PRECISION
 
 m    = 0
 200 DO  j = row1,row2
   sine = SIN(j)
   cosine = COS(j)
   m    = m + 1
   IF (sine == 0.0D0) GO TO 210
   x    = o(row+1)*cosine + o(j)*sine
   y    = da(m)   *sine   + o(j)*cosine
   z    = x       *cosine + y   *sine
   o(j) = y       *cosine - x   *sine
   da(m)= o(row+1)+ da(m) - z
   o(row+1) = z
   210 IF (j == n) CYCLE
   jp1  = j + 1
   DO  i = jp1,n
     m    = m + 1
     x    = da(m)*cosine - o(i)*sine
     o(i) = o(i)*cosine  + da(m)*sine
     y    = COS(i)*o(j)  + SIN(i)*x
     da(m)= COS(i)*x     - SIN(i)*o(j)
     o(j) = y
   END DO
 END DO
 RETURN
END SUBROUTINE rotate
