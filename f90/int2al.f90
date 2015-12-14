SUBROUTINE int2al (INT,alf,ch)
!     ----------
!     THIS ROUTINE CONVERTS AN INTEGER TO ALPHA-NUMERIC WORD. THE
!     NUMBER IS LEFT JUSTIFIED WITH NO BLANKS.
 
!     INPUT/OUTPUT
 
!     INT - INTEGER - INPUT - NOT CHANGED
!     ALF - BCD 2 WORDS - OUTPUT - 2A4 MAY BE USED FOR PRINTING
!     CH  - BCD 9 WORDS - OUTPUT - CH(1) .EQ. NUMBER OF CHARACTERS
!           NEEDED TO CREATE INT. MAY BE PRINTED BY CH(I), I=2,CH(1)
!           IN A1 FORMAT.
 
!     NOTE - ANY INPUT NUMBER OUTSIDE THE RANGE OF -9999999 AND +9999999
!            (I.E. MORE THAN 8 DIGITS) IS SET TO ZERO IN OUTPUT.
!     ----------
 
 
 INTEGER, INTENT(IN)                      :: INT
 INTEGER, INTENT(OUT)                     :: alf(2)
 INTEGER, INTENT(OUT)                     :: ch(9)
 INTEGER :: zero,   BLANK
 CHARACTER (LEN=8) :: k8
 DATA        BLANK,  zero /   1H ,    1H0 /
 
 IF (INT < -9999999 .OR. INT > +99999999) GO TO 50
 CALL int2k8 (*50,INT,k8)
 READ (k8,10) alf
 READ (k8,20) (ch(j),j=2,9)
 10   FORMAT (2A4)
 20   FORMAT (8A1)
 DO  j=2,9
   IF (ch(j) == BLANK) GO TO 40
 END DO
 j=10
 40   ch(1)=j-2
 RETURN
 
 50   ch(1) =1
 ch(2) =zero
 alf(1)=zero
 alf(2)=BLANK
 RETURN
END SUBROUTINE int2al
