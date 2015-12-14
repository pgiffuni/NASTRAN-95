FUNCTION iunion(i1,i2)
     
!     I1 AND I2 ARE .GT. 0 BUT .LE. 654321 AND CONSIST OF ANY UNIQUE
!     COMBINATION OF THE DIGITS 1 THRU 6
 
 
 
 INTEGER, INTENT(IN)                      :: i1
 INTEGER, INTENT(IN)                      :: i2
 INTEGER :: k(6,2) , kk(6) , r
 
!              DECODE I1 INTO K(*,1)
 i=1
 ii=i1
 ASSIGN 10 TO r
 GO TO 100
 
!              DECODE I2 INTO K(*,2)
 10 i=2
 ii=i2
 ASSIGN 20 TO r
 GO TO 100
 
!              FORM UNION OF K(*,1) AND K(*,2) IN KK(*)
 20 DO  i=1,6
   kk(i)=0
   IF(k(i,1) == i .OR. k(i,2) == i) kk(i)=i
 END DO
 
!              PACK KK(*) INTO IUNION
 j=1
 l=0
 DO  i=1,6
   IF(kk(i) == 0) CYCLE
   IF(l > 0) j=10*j
   l=l+j*i
 END DO
 
 iunion=l
 
 RETURN
 
 
 100 DO  j=1,6
   k(j,i)=0
 END DO
 DO  j=1,6
   l=ii-10*(ii/10)
   ii=(ii-l)/10
   IF(l == 0) GO TO 130
   k(l,i)=l
 END DO
 130 GO TO r,(10,20)
 
END FUNCTION iunion
