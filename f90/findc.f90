INTEGER FUNCTION findc (b,bbar,n,ix,jx)
     
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN OUT)                  :: bbar
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: ix(1)
 INTEGER, INTENT(OUT)                     :: jx(1)
 
 
!*******
!     PICK OUT PAIRS OF NUMBERS FOR ACTIVE ROWS
!*******
 icc = 0
 j = 1
 DO  i=1,n
   IF (i-ix(i) <= bbar) CYCLE
   jx(j) = i+b-1
   jx(j+1) = ix(i)
   j = j+2
 END DO
 j = j-1
 IF(j == 0) GO TO 31
 DO  k = 1,j,2
   IF((j-k-1)/2 < icc) EXIT
   ic = 0
   DO  l=k,j,2
     IF(jx(k) < jx(l+1)) CYCLE
     ic = ic+1
   END DO
   icc = MAX0(icc,ic)
 END DO
 31 findc = icc
 RETURN
END FUNCTION findc
