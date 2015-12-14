SUBROUTINE egnvct (c1,c2,eigen,c3,n1,n2,n)
     
!     SUBROUTINE TO OBTAIN EIGENVECTOR FROM REAL NON-SYMMETRIC
!     MATRICES FOR WHICH THE EIGENVALUE IS KNOWN.  THE METHOD
!     USED IS THE DIRECT METHOD OUTLINED IN ERR-FW-   BY DR.
!     A. M. CUNNINGHAM.
 
 
 COMPLEX, INTENT(OUT)                     :: c1(n,n)
 COMPLEX, INTENT(OUT)                     :: c2(n)
 COMPLEX, INTENT(IN)                      :: eigen
 COMPLEX, INTENT(OUT)                     :: c3(n)
 INTEGER, INTENT(OUT)                     :: n1(n)
 INTEGER, INTENT(OUT)                     :: n2(n)
 INTEGER, INTENT(IN)                      :: n
 
 COMPLEX :: d1,d2,d3,d4,d5,d6,d8
 
 ii3 = n
 ii2 = n - 1
 x1  = 0.0
 DO  j = 1,n
   n1(j) = j
   n2(j) = j
   c1(j,j) = c1(j,j) - eigen
   DO  i = 1,n
     x2 = cabs(c1(i,j))
     IF (x1-x2 < 0.0) THEN
       GO TO     5
     ELSE
       GO TO    10
     END IF
     5 x1 = x2
     i1 = i
     j1 = j
   END DO
 END DO
 DO  k6 = 2,n
   IF (cabs(c1(i1,j1)) == 0.0) THEN
     GO TO    30
   ELSE
     GO TO    50
   END IF
   30 k5 = k6 - 1
   
!     SINGULAR MATRIX RETURN ZERO
   
   DO  i = 1,n
     c3(i) = 0.0
   END DO
   GO TO 250
   
   50 d1 = (1.0,0.0)/c1(i1,j1)
   d2 = c1(i1,ii3)
   d3 = c1(ii3,j1)
   d4 = c1(ii3,ii3)
   DO  i = 1,ii2
     c3(i    ) = c1(i,j1)
     c1(i,j1 ) = c1(i,ii3)
     c1(i,ii3) =-c3(i)*d1
     d5 = -c1(i1,i)*d1
     c1(i1 ,i) = c1(ii3,i)
     c1(ii3,i) = d5
   END DO
   c3(i1) = d3
   c1(i1  ,j1) = d4
   c1(ii3 ,j1) =-d2*d1
   c1(i1 ,ii3) =-d3*d1
   c1(ii3,ii3) = d1
   IF (ii3 == n) GO TO 80
   ii4 = ii3 + 1
   DO  i = ii4,n
     d6 = c1(i1,i)
     c1(i1 ,i) = c1(ii3,i)
     c1(ii3,i) = d6
     c3(i    ) = c1(i,j1)
     c1(i,j1 ) = c1(i,ii3)
     c1(i,ii3) = c3(i)
   END DO
   80 i = n1(j1)
   n1(j1 ) = n1(ii3)
   n1(ii3) = i
   i = n2(i1)
   n2(i1 ) = n2(ii3)
   n2(ii3) = i
   x1 = 0.0
   DO  j = 1,ii2
     d8 = c1(ii3,j)
     DO  i = 1,ii2
       c1(i,j) = c1(i,j) + c3(i)*d8
       x2 = cabs(c1(i,j))
       IF(x1-x2 < 0.0) THEN
         GO TO   120
       ELSE
         GO TO   130
       END IF
       120 x1 = x2
       i1 = i
       j1 = j
     END DO
   END DO
   ii3 = ii3 - 1
   ii2 = ii2 - 1
 END DO
 
 c3(2) = c1(2,1)
 c3(1) = (1.0,0.0)
 DO  j = 3,n
   c3(j) = (0.0,0.0)
   j1 = j - 1
   DO  i = 1,j1
     c3(j) = c3(j) + c3(i)*c1(j,i)
   END DO
 END DO
 IF (cabs(c1(1,1)) < 1.0E-20) GO TO 202
 DO  k6 = 1,2
   
   loop184:  DO  j = 1,n
     i1 = n2(j)
     DO  i = 1,n
       IF (i1 == n1(i)) GO TO 184
     END DO
     c2(j) = c3(i)
   END DO loop184
   
   DO  j = 2,n
     i1 = n - j + 1
     j1 = i1 + 1
     DO  i = 1,i1
       c2(i) = c2(i) + c1(i,j1)*c2(j1)
     END DO
   END DO
   d1 = c1(1,1)/c2(1)
   c3(1) = (1.0,0.0)
   DO  j = 2,n
     i1 = j - 1
     c3(j) = c2(j)*c1(j,j)*d1
     DO  i = 1,i1
       c3(j) = c3(j) + c1(j,i)*c3(i)
     END DO
   END DO
 END DO
 
!     C3(I) NOW CONTAINS THE EIGENVECTOR WHICH MUST BE RE-ARRANGED
!     ACCORDING TO THE ORDER DICTATED BY N1(I) BACK TO THE ORIGINAL
!     ORDER.
 
 202 DO  i = 1,n
   i1 = n1(i)
   n1(i) = i
   205 IF (i1-i == 0) THEN
     GO TO   230
   END IF
   210 d1 = c3(i1)
   c3(i1) = c3(i)
   c3(i ) = d1
   k = n1(i1)
   n1(i1) = i1
   i1 = k
   GO TO 205
 END DO
 n1(1) = 2
 
 250 RETURN
END SUBROUTINE egnvct
