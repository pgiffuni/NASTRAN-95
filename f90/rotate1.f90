SUBROUTINE rotate1 (a,row,row1,row2, o,SIN,COS)
     
!     ROTATION OF A MATRIX PARTITION.
!     THIS ROUTINE IS CALLED ONLY BY TRIDI SUBROUTINE, WHICH IS CALLED
!     ONLY BY VALVEC
 
 
 REAL, INTENT(IN OUT)                     :: a(1)
 INTEGER, INTENT(IN OUT)                  :: row
 INTEGER, INTENT(IN)                      :: row1
 INTEGER, INTENT(IN)                      :: row2
 REAL, INTENT(IN OUT)                     :: o(1)
 REAL, INTENT(IN)                         :: SIN(1)
 REAL, INTENT(IN)                         :: COS(1)
 
 
 REAL :: sine,cosine,x,y,z
 COMMON /givn  /  title(100),n
 
!     O     = 2ND ROW OF THE COMPLETE MATRIX.
!     SIN   = SINES.
!     COS   = COSINES.
!     A  = MATRIX PARTITION (TRIANGULAR) - SINGLE PRECISION
 
 m    = 0
 DO  j = row1,row2
   sine = SIN(j)
   cosine = COS(j)
   m    = m + 1
   IF (sine == 0.) GO TO 101
   x    = o(row+1)*cosine + o(j)*sine
   y    = a(m)    * sine  + o(j)*cosine
   z    = x       *cosine + y   *sine
   o(j) = y       *cosine - x   *sine
   a(m) = o(row+1) + a(m) - z
   o(row+1) = z
   101 IF (j == n) CYCLE
   jp1  = j + 1
   DO   i = jp1,n
     m    = m + 1
     x    = a(m)*cosine - o(i)*sine
     o(i) = o(i)*cosine + a(m)*sine
     y    = COS(i)*o(j) + SIN(i)*x
     a(m) = COS(i)*x    - SIN(i)*o(j)
     o(j) = y
   END DO
 END DO
 RETURN
END SUBROUTINE rotate1
