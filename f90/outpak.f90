SUBROUTINE outpak (ii,iout,isn)
     
 
 INTEGER, INTENT(OUT)                     :: ii
 INTEGER, INTENT(OUT)                     :: iout(1)
 INTEGER, INTENT(OUT)                     :: isn
 EXTERNAL        orf
 INTEGER :: NUMBER(10),idig(4),orf,op
 COMMON /system/ ksys(65)
 EQUIVALENCE     (ksys(2),op),(ksys(9),nlpp),(ksys(12),nline), (ksys(41),ncpw)
 DATA    NUMBER/ 1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H0/
 DATA    nblank/ 4H    /
 
 kcode = 0
 IF (isn < 0) kcode = 1
 isn = IABS(isn)
 icode = 0
 IF (ii > 32) icode = 1
 IF (icode == 1) GO TO 50
 
!     TRANSLATE ISN TO DIGITS
 
 10 idig(1) = isn/1000
 idig(2) = (isn-idig(1)*1000)/100
 idig(3) = (isn-idig(1)*1000 - idig(2)*100)/10
 idig(4) = isn - idig(1)*1000 - idig(2)*100 - idig(3)*10
 DO  i = 1,4
   IF (idig(i) == 0) idig(i) = 10
 END DO
 
!     FORM WORD AND STORE IN IOUT ARRAY
 
 k = 0
 DO  i = 1,4
   j = idig(i)
   k = orf(klshft(krshft(NUMBER(j),ncpw-1),ncpw-i),k)
 END DO
 iout(ii) = k
 GO TO 80
 50 nline = nline + 1
 IF (nline <= nlpp) GO TO 60
 CALL page
 nline = nline + 1
 WRITE (op,100)
 nline = nline + 1
 60 WRITE (op,90) (iout(i),i=2,32)
 ii = 5
 IF (kcode == 1) ii = 7
 DO  ll = 2,32
   iout(ll) = nblank
 END DO
 GO TO 10
 80 RETURN
 
 90 FORMAT (5X,31A4)
 100 FORMAT (/,1H )
END SUBROUTINE outpak
