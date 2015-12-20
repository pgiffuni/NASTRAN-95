SUBROUTINE invers (ndim,a,n,b,m,determ,ising,INDEX)
     
!     INVERSE, OR LINEAR EQUATIONS SOLVER
 
!     NDIM IS THE ACTUAL SIZE OF A IN CALLING PROGRAM. E.G. A(NDIM,NDIM)
!     A IS SQUARE MATRIX TO BE INVERTED.
!     N IS SIZE OF UPPER LEFT PORTION BEING INVERTED.
!     B IS COLUMN OF CONSTANTS (OPTIONAL INPUT). SUPPLY SPACE B(NDIM,1)
!     M IS THE NUMBER OF COLUMNS OF CONSTANTS
!     DETERM RETURNS THE VALUE OF DETERMINANT IF NON-SINGULAR
!     ISING RETURNS 2, IF MATRIX A(N,N) IS SINGULAR
!                   1, IF MATRIX A(N,N) IS NON-SINGULAR
!     (IF ISING IS SET TO .LT. 0 UPON INPUT, DETERM IS NOT CALCULATED)
!     INVERSE RETURNS  IN A
!     SOLUTION VECTORS RETURN IN B
!     INDEX IS WORKING STORAGE (N,3)
 
 
 INTEGER, INTENT(IN OUT)                  :: ndim
 REAL, INTENT(IN OUT)                     :: a(ndim,1)
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN OUT)                     :: b(ndim,1)
 INTEGER, INTENT(IN)                      :: m
 REAL, INTENT(OUT)                        :: determ
 INTEGER, INTENT(OUT)                     :: ising
 INTEGER, INTENT(OUT)                     :: INDEX(n,3)
 
 COMMON /machin/ mach
 EQUIVALENCE     (irow,jrow), (icolum,jcolum), (amax, t, swap)
 DATA    epsi  / 1.0E-30 /
 
!     INITIALIZE
 
 IF (mach  == 5) epsi = 1.e-18
 determ = 1.0
 IF (ising < 0) determ = 0.0
 DO  j = 1,n
   INDEX(j,3) = 0
 END DO
 DO  i = 1,n
   
!     SEARCH FOR PIVOT
   
   amax = 0.0
   DO  j = 1,n
     IF (INDEX(j,3) == 1) CYCLE
     DO  k = 1,n
       IF (INDEX(k,3) - 1 < 0) THEN
         GO TO    20
       ELSE IF (INDEX(k,3) - 1 == 0) THEN
         GO TO    30
       ELSE
         GO TO   190
       END IF
       20 IF (ABS(a(j,k)) <= amax) CYCLE
       irow   = j
       icolum = k
       amax   = ABS(a(j,k))
       30 CONTINUE
     END DO
   END DO
   INDEX(icolum,3) = INDEX(icolum,3) + 1
   INDEX(i,1) = irow
   INDEX(i,2) = icolum
   
!     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
   
   IF (irow == icolum) GO TO 70
   determ = -determ
   DO  l = 1,n
     swap = a(irow,l)
     a(irow,l  ) = a(icolum,l)
     a(icolum,l) = swap
   END DO
   IF (m <= 0) GO TO 70
   DO  l = 1,m
     swap = b(irow,l)
     b(irow,l  ) = b(icolum,l)
     b(icolum,l) = swap
   END DO
   
!     DIVIDE PIVOT ROW BY PIVOT ELEMENT
   
   70 pivot  = a(icolum,icolum)
   determ = determ*pivot
   IF (ABS(pivot) < epsi) GO TO 190
   a(icolum,icolum) = 1.0
   DO  l = 1,n
     a(icolum,l) = a(icolum,l)/pivot
   END DO
   IF (m <= 0) GO TO 100
   DO  l=1,m
     b(icolum,l) = b(icolum,l)/pivot
   END DO
   
!     REDUCE NON PIVOT ROWS
   
   100 DO  l1 = 1,n
     IF (l1 == icolum) CYCLE
     t = a(l1,icolum)
     a(l1,icolum) = 0.0
     IF (ABS(t) < epsi) CYCLE
     DO  l = 1,n
       a(l1,l) = a(l1,l) - a(icolum,l)*t
     END DO
     IF (m <= 0) CYCLE
     DO  l = 1,m
       b(l1,l) = b(l1,l) - b(icolum,l)*t
     END DO
   END DO
 END DO
 
!     INTERCHANGE COLUMNS
 
 DO  i = 1,n
   l = n + 1 - i
   IF (INDEX(l,1) == INDEX(l,2)) CYCLE
   jrow   = INDEX(l,1)
   jcolum = INDEX(l,2)
   DO  k = 1,n
     swap = a(k,jrow)
     a(k,jrow  ) = a(k,jcolum)
     a(k,jcolum) = swap
   END DO
 END DO
 DO  k = 1,n
   IF (INDEX(k,3) == 1) GO TO 160
   ising = 2
   GO TO 180
   160 CONTINUE
 END DO
 ising = 1
 180 RETURN
 190 ising = 2
 RETURN
END SUBROUTINE invers
