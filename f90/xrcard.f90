SUBROUTINE xrcard (out,nflag,in)
     
!     MODIFIED BY G.CHAN/UNISYS FOR EFFICIENCY,        2/1988
!     LAST REVISED, 8/1989, IMPROVED EFFICIENCY BY REDUCING CHARACTER
!     OPERATIONS (VERY IMPORTANT FOR CDC MACHINE)
 
!     REVERT TO NASTRAN ORIGIANL XRCARD ROUTINE (NOW CALLED YRCARD)
!     IF DIAG 42 IS TURNED ON
!     (THIS NEW XRCARD IS SEVERAL TIMES FASTER THAN YRCARD)
 
 
 
 INTEGER, INTENT(OUT)                     :: out(2)
 INTEGER, INTENT(IN OUT)                  :: nflag
 INTEGER, INTENT(IN OUT)                  :: in(18)
 EXTERNAL        lshift   ,rshift   ,complf
 LOGICAL :: alpha    ,delim    ,expont   ,power    ,minus   , nogo     ,debug
 INTEGER :: idoubl(2),TYPE(72) ,nt(15)
 INTEGER :: plus1    ,minus1   ,blank1   ,dot1     ,e1      ,  &
     d1       ,num1(10) ,char1(72),dollr1   ,comma1  ,  &
     equal1   ,slash1   ,oparn1   ,cparn1   ,astk1   ,  &
     char1n   ,char1i   ,charn2   ,chars(23),zero1   , precis   ,psign
 REAL :: flpt
 DOUBLE PRECISION :: ddoubl
 CHARACTER (LEN=1) :: save1(8) ,khar1(72)
 CHARACTER (LEN=8) :: save8    ,blank8   ,numric
 CHARACTER (LEN=23) :: char23   ,ufm
 CHARACTER (LEN=72) :: char72
 COMMON /xmssg / ufm
 COMMON /lhpwx / lowpw    ,highpw
 COMMON /system/ bufsz    ,nout     ,nogo
 EQUIVALENCE     (khar1(1),char72), (save8 ,save1( 1)),  &
     (flpt    ,intg  ), (ddoubl,idoubl(1)),  &
     (chars(1),dollr1), (chars( 8),cparn1),  &
     (chars(2),plus1 ), (chars( 9),e1    ),  &
     (chars(3),equal1), (chars(10),d1    ),  &
     (chars(4),minus1), (chars(11),dot1  ),  &
     (chars(5),comma1), (chars(12),blank1),  &
     (chars(6),slash1), (chars(13),astk1 ),  &
     (chars(7),oparn1), (chars(14),num1(1)), (num1(10),zero1 )
 
 DATA    char23/ '$+=-,/()ED. *1234567890'/,   blank8/ '       ' /
 DATA    dollr1,  blank4,  diag  ,   debug ,   numric            /  &
     0,       4H    ,  4HDIAG,  .false.,  'NUMERIC'          /
 DATA    equal4,  slash4,  oparn4,   astk4                       /  &
     4H=   ,  4H/   ,  4H(   ,   4H*                         /
 
 IF (dollr1 /= 0) GO TO 20
 CALL k2b (char23,chars,23)
 a77777 = complf(0)
 a67777 = rshift(lshift(a77777,1),1)
 prev   = blank4
 CALL sswtch (42,l42)
 IF (debug) WRITE (nout,10)
 10 FORMAT (//5X,'INPUT DEBUG IN XRCARD ROUTINE')
 
!         KXX=0
 
!     USE ORIGINAL XRCARD ROUTINE IF DIAG 42 IS TURNED ON
 
 20 IF (prev == diag) CALL sswtch (42,l42)
 IF (l42 == 0) GO TO 30
 CALL yrcard (out,nflag,in)
 RETURN
 
!     CONVERT 18 BCD WORDS IN 'IN' TO 72 CHARACTER STRING, AND SET TYPE
 
 30 CALL bcdkh7 (in,char72)
 CALL k2b (char72,char1,72)
 IF (debug) WRITE (nout,45) char72
 45 FORMAT (/,' INPUT- ',a72)
 
 
 DO  n = 1,72
   char1n = char1(n)
   IF (char1n == blank1 ) GO TO 60
   DO  k = 1,10
     IF (char1n == num1(k)) GO TO 70
   END DO
   TYPE(n) =-1
   CYCLE
   60 TYPE(n) = 0
   CYCLE
   70 TYPE(n) = 1
 END DO
 n = 73
 90 n = n - 1
 IF (TYPE(n) == 0) GO TO 90
 last  = n + 1
 IF (last > 72) last = 72
 alpha = .false.
 delim = .true.
 iout  = 0
 n     = 0
 asave = 1
 out(asave) = 0
 save8 = blank8
 100 IF (n == last) GO TO 510
 IF (nflag-iout < 5) GO TO 660
 minus = .false.
 n = n + 1
 char1n = char1(n)
 IF (TYPE(n) < 0) THEN
   GO TO   110
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   100
 ELSE
   GO TO   210
 END IF
!                  BCD  BLANK NUMERIC
 
 110 IF (char1n == plus1 .OR. char1n == minus1 .OR. char1n == dot1) GO TO 200
 IF (char1n == dollr1) GO TO 180
 
!     GOOD ALPHA FIELD OR DELIMITER
 
 IF (alpha) GO TO 120
 IF ((char1n == comma1 .OR. char1n == dollr1) .AND. (.NOT.delim)) GO TO 180
 IF (char1n == cparn1 .AND. .NOT.delim) GO TO 180
 iout  = iout + 1
 asave = iout
 out(asave) = 0
 alpha = .true.
 120 IF  (char1n == oparn1 .OR. char1n == slash1 .OR. char1n == equal1  &
     .OR. char1n == comma1 .OR. char1n == astk1  .OR. char1n == dollr1) GO TO 180
 IF (char1n == cparn1) GO TO 180
 ASSIGN 125 TO irtn
 imhere = 125
 GO TO 170
 125 out(asave) = out(asave) + 1
 iout  = iout + 2
 delim = .false.
 out(iout-1) = blank4
 out(iout  ) = blank4
 ICHAR = 0
 GO TO 150
 130 IF (n == last) GO TO 510
 n = n + 1
 char1n = char1(n)
 IF (TYPE(n) < 0) THEN
   GO TO   140
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   160
 ELSE
   GO TO   150
 END IF
!                  BCD BLANK NUMERIC
 
 140 IF  (char1n == oparn1 .OR. char1n == slash1 .OR. char1n == equal1  &
     .OR. char1n == comma1 .OR. char1n == astk1  .OR. char1n == dollr1) GO TO 180
 IF (char1n == cparn1) GO TO 180
 
!     RECONSTRUCT CHARACTERS INTO SAVE1 SPACE, UP TO 8 CHARACTERS ONLY
 
 150 IF (ICHAR == 8) GO TO 130
 ICHAR = ICHAR + 1
 save1(ICHAR) = khar1(n)
 imhere = 150
 IF (debug) WRITE (nout,171) save8,imhere,ICHAR,iout
 
!     GO FOR NEXT CHARACTER
 
 IF (ICHAR < 8) GO TO 130
 ASSIGN 130 TO irtn
 imhere = 155
 GO TO 170
 
!     A BLANK CHARACTER IS ENCOUNTERED WHILE PROCESSING ALPHA STRING
!     IF THIS IS AT THE BEGINNING OF A NEW BCD WORD, GO TO 100
!     IF THIS IS AT THE END OF A BCD WORD, GO TO 170 TO WRAP IT UP
 
 160 IF (ICHAR == 0) GO TO 100
 ASSIGN 100 TO irtn
 imhere = 160
 
!     MOVE CHARACTER DATA IN SAVE8 TO OUT(IOUT-1) AND OUT(IOUT) IN BCD
!     WORDS
 
 170 IF (save8 /= blank8) CALL khrbc2 (save8,out(iout-1))
 IF (.NOT.debug) GO TO 175
 WRITE (nout,171) save8,imhere,ICHAR,iout
 WRITE (nout,172) iout,out(iout-1),out(iout),delim, iout,out(iout-1),out(iout)
 171 FORMAT ('   SAVE8= /',a8,'/   @',i3,',  ICHAR,IOUT=',2I3)
 172 FORMAT ('   IOUT,OUT  =',i4,2H /,2A4,'/  DELIM=',l1,  &
     /14X,    '=',i4,2H /,2I25,'/')
 175 save8 = blank8
 GO TO irtn, (100,125,130,185)
 
!     DELIMITER HIT
 
 180 ASSIGN 185 TO irtn
 imhere = 180
 GO TO 170
 185 IF (.NOT. delim) GO TO 190
 IF (iout == 0) iout = 1
 iout = iout + 2
 out(asave)  = out(asave) + 1
 out(iout-1) = blank4
 out(iout  ) = blank4
 190 IF (char1n == dollr1) GO TO 520
 delim = .true.
 IF (char1n == cparn1) delim = .false.
 IF (char1n == comma1) GO TO 100
 IF (char1n == cparn1) GO TO 100
 
!     OUTPUT DELIMITER
 
 iout = iout + 2
 out(asave) = out(asave) + 1
 out(iout) = blank4
 IF (char1n == oparn1) out(iout) = oparn4
 IF (char1n == slash1) out(iout) = slash4
 IF (char1n == equal1) out(iout) = equal4
 IF (char1n ==  astk1) out(iout) =  astk4
 IF (out(iout) == blank4) GO TO 590
 out(iout-1) = a77777
 IF (debug) WRITE (nout,195) iout,out(iout),delim,char1n
 195 FORMAT (5X,'IOUT,OUT/@195 =',i4,2H ',A4,8H' delim=,l1,2H ',A1, 1H')
 save8 = blank8
 GO TO 100
 
!     PLUS, MINUS, OR DOT ENCOUNTERED
 
 200 IF (char1n == minus1) minus = .true.
 IF (char1n /=   dot1) n = n + 1
 IF (n > last) GO TO 530
 
!     NUMERIC
 
 210 alpha = .false.
 delim = .false.
 it    = 0
 nt(1) = 0
 DO  i = n,last
   IF (TYPE(i) < 0) THEN
     GO TO   290
   ELSE IF (TYPE(i) == 0) THEN
     GO TO   270
   END IF
   
!     INTEGER CHARACTER
   
   220 char1i = char1(i)
   DO  k = 1,9
     IF (char1i == num1(k)) GO TO 250
   END DO
   k  = 0
   250 it = it + 1
   IF (it < 16) nt(it) = k
 END DO
 
!     FALL HERE IMPLIES WE HAVE A SIMPLE INTEGER
 
 270 NUMBER = 0
 DO  i = 1,it
   IF (((a67777-nt(i))/10) < NUMBER) GO TO 550
   NUMBER = NUMBER*10  +  nt(i)
 END DO
 IF (minus) NUMBER = - NUMBER
 iout = iout + 2
 out(iout-1) =-1
 out(iout  ) = NUMBER
 IF (.NOT.debug) GO TO 285
 imhere = 280
 WRITE (nout,171) numric,imhere
 WRITE (nout,282) iout,out(iout-1),out(iout),delim
 282 FORMAT (10X,i4,1H),2I8,'    DELIM=',l1)
 285 n = n + it - 1
 GO TO 100
 
!     FLOATING PT. NUMBER, DELIMITER, OR ERROR IF FALL HERE
 
!     COUNT THE NUMBER OF DIGITS LEFT BEFORE CARD END OR DELIMITER HIT
 
 290 n1 = i
 DO  n2 = n1,last
   charn2 = char1(n2)
   IF (charn2 == oparn1 .OR. charn2 == slash1 .OR.  &
       charn2 == equal1 .OR. charn2 == comma1 .OR.  &
       charn2 == dollr1 .OR. TYPE(n2) == 0) GO TO 310
   IF (charn2 == cparn1) GO TO 310
 END DO
 n2 = last + 1
 310 IF (n1 == n2) GO TO 270
 
!     CHARACTER N1 NOW MUST BE A DECIMAL FOR NO ERROR
 
 IF (char1(n1) /= dot1) GO TO 570
 power = .false.
 n1 = n1 + 1
 n2 = n2 - 1
 places = 0
 expont = .false.
 ipower = 0
 psign  = zero1
 precis = zero1
 IF (n2 < n1) GO TO 410
 DO  i = n1,n2
   char1i = char1(i)
   IF (TYPE(i) < 0) THEN
     GO TO   360
   ELSE IF (TYPE(i) == 0) THEN
     GO TO   570
   END IF
   
!     FLOATING PT. NUMBER
   
   320 DO  k = 1,9
     IF (char1i == num1(k)) GO TO 340
   END DO
   k = 0
   340 IF (expont) GO TO 350
   it = it + 1
   IF (it < 16) nt(it) = k
   places = places + 1
   CYCLE
   
!     BUILD POWER HERE
   
   350 power  = .true.
   ipower = ipower*10 + k
   IF (ipower > 1000) GO TO 630
   CYCLE
   
!     START EXPONENTS HERE
   
   360 IF (expont) GO TO 380
   expont = .true.
   IF (char1i /= plus1 .AND. char1i /= minus1) GO TO 370
   precis = e1
   psign  = char1i
   GO TO 390
   370 IF (char1i /= e1 .AND. char1i /= d1) GO TO 610
   precis = char1i
   GO TO 390
   
!     SIGN OF POWER
   
   380 IF (power) GO TO 610
   IF (psign /= zero1 .OR.  &
       (char1i /= plus1 .AND. char1i /= minus1)) GO TO 610
   psign = char1i
   power = .true.
   390 IF (i == last) GO TO 530
 END DO
 410 n = n2
 
!     ALL DATA COMPLETE FOR FLOATING POINT NUMBER
!     ONLY 15 FIGURES WILL BE ACCEPTED
 
 IF (it <= 15) GO TO 420
 ipower = ipower + it - 15
 it = 15
 420 IF (psign == minus1) ipower = -ipower
 ipower = ipower - places
 NUMBER = 0
 IF (it < 7) GO TO 430
 n2 = 7
 GO TO 440
 430 n2 = it
 440 DO  i = 1,n2
   NUMBER = NUMBER*10 + nt(i)
 END DO
 ddoubl = DBLE(FLOAT(NUMBER))
 IF (it <= 7) GO TO 470
 NUMBER = 0
 n2 = it - 7
 DO  i = 1,n2
   it = i + 7
   NUMBER = NUMBER*10 + nt(it)
 END DO
 ddoubl = ddoubl*10.0D0**n2 + DBLE(FLOAT(NUMBER))
 470 IF (minus) ddoubl = -ddoubl
 
!     POWER HAS TO BE WITHIN RANGE OF MACHINE
 
 ichek = ipower + it
 IF (ddoubl == 0.0D0) GO TO 490
 IF (ichek < lowpw+1 .OR. ichek > highpw-1 .OR.  &
     ipower < lowpw+1 .OR. ipower > highpw-1) GO TO 640
 ddoubl = ddoubl*10.0D0**ipower
 490 IF (precis == d1) GO TO 500
 flpt = ddoubl
 iout = iout + 2
 out(iout-1) =-2
 out(iout  ) = intg
 GO TO 100
 500 iout = iout + 3
 out(iout-2) =-4
 out(iout-1) = idoubl(1)
 out(iout  ) = idoubl(2)
 GO TO 100
 
!     PREPARE TO RETURN
 
 510 IF (.NOT. delim) GO TO 520
 IF (save8 /= blank8) CALL khrbc2 (save8,out(iout-1))
 out(iout+1) = 0
 GO TO 525
 520 out(iout+1) = a67777
 525 prev = out(2)
 RETURN
 
!     ERRORS
 
 530 WRITE  (nout,540) ufm
 540 FORMAT (a23,'300, INVALID DATA COLUMN 72')
 GO TO  680
 550 WRITE  (nout,560) ufm
 560 FORMAT (a23,'300, INTEGER DATA OUT OF MACHINE RANGE')
 GO TO  680
 570 WRITE  (nout,580) ufm,n1
 580 FORMAT (a23,'300, INVALID CHARACTER FOLLOWING INTEGER IN COLUMN', i4)
 GO TO  680
 590 WRITE  (nout,600) ufm,char1n
 600 FORMAT (a23,'300, FORGOTTEN DELIMITER - ',a1,',  PROGRAM ERROR')
 GO TO  680
 610 WRITE  (nout,620) ufm,i
 620 FORMAT (a23,'300, DATA ERROR-UNANTICIPATED CHARACTER IN COLUMN', i4)
 GO TO  680
 630 CONTINUE
 640 WRITE  (nout,650) ufm
 650 FORMAT (a23,'300, DATA ERROR - MISSING DELIMITER OR REAL POWER ',  &
     'OUT OF MACHINE RANGE')
 GO TO  680
 660 WRITE  (nout,670) ufm
 670 FORMAT (a23,'300, ROUTINE XRCARD FINDS OUTPUT BUFFER TOO SMALL ',  &
     'TO PROCESS CARD COMPLETELY')
 680 nogo = .true.
 WRITE  (nout,690) char72
 690 FORMAT (/5X,1H',A72,1H','  ERROR IN XRCARD ROUTINE')
 out(1) = 0
 RETURN
 
END SUBROUTINE xrcard
