SUBROUTINE incore(a,n,b,cx,ix)
     
!     IN-CORE DECOMPOSITION OF SQUARE, COMPLEX, NXN MATRIX,A.
!     AX = B.
!     CX = X
!     IX = NUMBER OF B VECTORS SPECIFIED.
 
 
 COMPLEX, INTENT(IN OUT)                  :: a(n,n)
 INTEGER, INTENT(IN)                      :: n
 COMPLEX, INTENT(IN OUT)                  :: b(ix,n)
 COMPLEX, INTENT(OUT)                     :: cx(ix,n)
 INTEGER, INTENT(IN)                      :: ix
 COMPLEX :: cmax, scrch
 COMPLEX :: t1,t2,t3
 COMPLEX :: csum
 
 
 IF(n == 2) GO TO 500
 IF(n == 1) GO TO 600
 nm1 = n-1
 
!     PIVOT MAYBE.
 
 DO  j=1,nm1
   cmax = a(j,j)
   jp1 = j + 1
   jmax = j
   DO  jj= jp1,n
     IF(cabs(a(j,jj)) <= cabs(cmax)) CYCLE
     cmax = a(j,jj)
     irow = jj
     jmax = jj
   END DO
   
!     IROW = ROW WITH LARGEST ELEMENT IN COLUMN J.
!     MOVE PIVOT ROW TO TOP OF ELIMINATION
   
   amax = cabs(cmax)
   IF(amax == 0.) CYCLE
   IF(jmax == j) GO TO 120
   DO  jj= j,n
     scrch = a(jj,j)
     a(jj,j) = a(jj,irow)
     a(jj,irow) = scrch
   END DO
   
!     INTERCHANGE B VECTOR
   
   DO  jj = 1,ix
     scrch = b(jj,j)
     b(jj,j) = b(jj,irow)
     b(jj,irow) = scrch
   END DO
   
!     ELIMINATE COLUMN
   
   120 CONTINUE
   a(j,j) = (1.0,0.0) / a(j,j)
   t1 = a(j,j)
   DO  i=jp1,n
     t2 = a(i,j)
     IF(cabs(t2) < (1.0E-19)) CYCLE
     t2 = -t2*t1
     a(i,j) = t2
     DO  l = jp1,n
       t3 = a(j,l)
       IF(cabs(t3) < (1.0E-19)) CYCLE
       a(i,l) = a(i,l) + t3*t2
     END DO
   END DO
   
!     HANDLE B ELIMINATION.
   
   DO  jj = 1,ix
     b(jj,j) = b(jj,j) * t1
   END DO
   DO  jj = 1,ix
     DO  k=jp1,n
       b(jj,k) = b(jj,k) - b(jj,j)*a(j,k)
     END DO
   END DO
 END DO
 
!     BACKWARD PASS.
 
 DO  jj = 1,ix
   cx(jj,n) = b(jj,n)/a(n,n)
 END DO
 DO  jj = 1,ix
   i = n
   190 CONTINUE
   csum = (0.,0.)
   k = i-1
   DO  j=i,n
     csum = csum + cx(jj,j)*a(j,k)
   END DO
   cx(jj, k)    = b(jj, k)    + csum
   IF(i <= 2) CYCLE
   i = i-1
   GO TO 190
 END DO
 RETURN
 500 CONTINUE
 DO  i = 1,ix
   cx(i ,2) = (b(i ,2)-(b(i ,1)*a(1,2)/a(1,1)))/(a(2,2)-(a(2,1)  &
       *a(1,2)/a(1,1)))
   cx(i ,1) = b(i ,1)/a(1,1)-a(2,1)*cx(i ,2)/a(1,1)
 END DO
 RETURN
 600 CONTINUE
 DO  i=1,ix
   cx(i ,1) = b(i ,1)/a(1,1)
 END DO
 RETURN
END SUBROUTINE incore
