SUBROUTINE yrcard (out,nflag,in)
     
!     THIS WAS NASTRAN ORIGINAL XRCARD ROUTINE, AND IS NOW RENAMED
!     YRCARD
!     THIS ROUTINE IS CALLED ONLY BY XRCARD
!     THIS ROUTINE CAN BE DELETED IF THE NEW XRCARD ROUTINE PASSES
!     ALL RELIABILITY TESTS                 G.CHAN/UNISYS,  2/1988
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(OUT)                     :: out(1)
 INTEGER, INTENT(IN OUT)                  :: nflag
 INTEGER, INTENT(IN)                      :: in(18)
 EXTERNAL         lshift,rshift,complf
 LOGICAL :: alpha,delim,expont,power,lminus,pass,nogo
 INTEGER :: f6
 REAL :: fl1
 DOUBLE PRECISION :: xdoubl
 DIMENSION        ndoubl(2),num(10),TYPE(72),CHAR(72),  nt(15),chars(13)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /lhpwx /  lowpw,highpw
 COMMON /system/  ibufsz,f6,nogo,dum1(7),npages,nlines
 EQUIVALENCE      (fl1      ,int1  ), (xdoubl,ndoubl(1)),  &
     (num(10)  ,zero  ), (chars( 1),dollar),  &
     (chars( 2),plus  ), (chars( 3),equal ),  &
     (chars( 4),minus ), (chars( 5),comma ),  &
     (chars( 6),slash ), (chars( 7),oparen),  &
     (chars( 8),cparen), (chars( 9),e     ),  &
     (chars(10),d     ), (chars(11),period),  &
     (chars(12),BLANK ), (chars(13),astk  )
 DATA    blanks/  4H     /, BLANK / 4H     /, dollar/ 4H$    /,  &
     equal /  1H=    /, astk  / 1H*    /, comma / 1H,    /,  &
     slash /  1H/    /, cparen/ 1H)    /, oparen/ 1H(    /,  &
     plus  /  1H+    /, minus / 1H-    /, period/ 1H.    /,  &
     e     /  1HE    /, d     / 1HD    /, pass  / .false./,  &
     num   /  1H1, 1H2, 1H3,1H4,1H5, 1H6, 1H7,1H8,1H9,1H0/
 
 IF (pass) GO TO 50
 pass   = .true.
 a77777 = complf(0)
 a67777 = rshift(lshift(a77777,1),1)
 
!     READ AND TYPE 72 CHARACTERS
 
 50 n = 0
 DO  i = 1,18
   DO  j = 1,4
     n = n + 1
     charac = khrfn1(blanks,1,in(i),j)
     IF (charac == BLANK) GO TO 70
     DO  k = 1,10
       IF (charac == num(k)) GO TO 80
     END DO
     TYPE(n) = -1
     GO TO 90
     70 TYPE(n) = 0
     GO TO 90
     80 TYPE(n) = 1
     90 CHAR(n) = charac
   END DO
 END DO
 alpha = .false.
 delim = .true.
 iout  = 0
 n     = 0
 asave = 1
 out(asave) = 0
 100 IF (n == 72) GO TO 510
 IF (nflag-iout < 5) GO TO 660
 lminus = .false.
 n = n + 1
 nchar = CHAR(n)
 IF (TYPE(n) < 0) THEN
   GO TO   110
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   100
 ELSE
   GO TO   210
 END IF
 110 IF (nchar == plus .OR. nchar == minus .OR. nchar == period) GO TO 200
 IF (nchar == dollar) GO TO 180
 
!     GOOD ALPHA FIELD OR DELIMETER
 
 IF (alpha) GO TO 120
 IF ((nchar == comma .OR. nchar == dollar) .AND. (.NOT.delim)) GO TO 180
 IF (nchar == cparen .AND. .NOT.delim) GO TO 180
 iout  = iout + 1
 asave = iout
 out(asave) = 0
 alpha = .true.
 120 IF (nchar == oparen .OR. nchar == slash .OR. nchar == equal .OR.  &
     nchar == comma  .OR. nchar == astk  .OR. nchar == dollar) GO TO 180
 IF (nchar == cparen) GO TO 180
 out(asave) = out(asave) + 1
 iout  = iout + 2
 delim = .false.
 out(iout-1) = blanks
 out(iout  ) = blanks
 ICHAR = 0
 GO TO 150
 130 IF (n == 72) GO TO 510
 n = n + 1
 nchar = CHAR(n)
 IF (TYPE(n) < 0) THEN
   GO TO   140
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   100
 ELSE
   GO TO   150
 END IF
 140 IF (nchar == oparen .OR. nchar == slash .OR. nchar == equal .OR.  &
     nchar == comma  .OR. nchar == astk  .OR. nchar == dollar) GO TO 180
 IF (nchar == cparen) GO TO 180
 150 IF (ICHAR == 8) GO TO 130
 ICHAR = ICHAR + 1
 IF (ICHAR <= 4) GO TO 160
 ipos = ICHAR - 4
 word = iout
 GO TO 170
 160 ipos = ICHAR
 word = iout - 1
 
!     CLEAR SPOT IN WORD FOR CHAR(N) AND PUT CHAR(N) IN IT
 
 170 out(word) = khrfn1(out(word),ipos,nchar,1)
 
!     GO FOR NEXT CHARACTER
 
 GO TO 130
 
 
!     DELIMETER HIT
 
 180 IF (.NOT. delim) GO TO 190
 IF (iout == 0) iout = 1
 iout = iout + 2
 out(asave)  = out(asave) + 1
 out(iout-1) = blanks
 out(iout  ) = blanks
 190 IF (nchar == dollar) GO TO 520
 delim = .true.
 IF (nchar == cparen) delim = .false.
 IF (nchar ==  comma) GO TO 100
 IF (nchar == cparen) GO TO 100
 
!     OUTPUT DELIMETER
 
 iout = iout + 2
 out(asave ) = out(asave) + 1
 out(iout-1) = a77777
 out(iout  ) = khrfn1(blanks,1,nchar,1)
 GO TO 100
 
 
 200 IF (nchar ==  minus) lminus = .true.
 IF (nchar /= period) n = n + 1
 IF (n > 72) GO TO 530
 
 210 alpha = .false.
 delim = .false.
 it = 0
 nt(1) = 0
 DO  i = n,72
   IF (TYPE(i) < 0) THEN
     GO TO   290
   ELSE IF (TYPE(i) == 0) THEN
     GO TO   270
   END IF
   
!     INTEGER CHARACTER
   
   220 DO  k = 1,9
     IF (CHAR(i) == num(k)) GO TO 250
   END DO
   k  = 0
   250 it = it + 1
   IF (it < 16) nt(it) = k
 END DO
 
!     FALL HERE IMPLIES WE HAVE A SIMPLE INTEGER
 
 270 NUMBER = 0
 DO  i = 1,it
   IF (((a67777-nt(i))/10) < NUMBER) GO TO 550
   NUMBER = NUMBER*10 + nt(i)
 END DO
 IF (lminus) NUMBER = - NUMBER
 iout = iout + 2
 out(iout-1) = -1
 out(iout  ) = NUMBER
 n = n + it - 1
 GO TO 100
 
!     REAL NUMBER, DELIMETER, OR ERROR IF FALL HERE
 
!     COUNT THE NUMBER OF DIGITS LEFT BEFORE CARD END OR DELIMETER HIT
 
 290 n1 = i
 DO  n2 = n1,72
   IF (CHAR(n2) == oparen .OR. CHAR(n2) == slash .OR.  &
       CHAR(n2) == equal  .OR. CHAR(n2) == comma .OR.  &
       CHAR(n2) == dollar .OR. TYPE(n2) == 0) GO TO 310
   IF (CHAR(n2) == cparen) GO TO 310
 END DO
 n2 = 73
 310 IF (n1 == n2) GO TO 270
 
!     CHARACTER N1 NOW MUST BE A DECIMAL FOR NO ERROR
 
 IF (CHAR(n1) /= period) GO TO 570
 power = .false.
 n1 = n1 + 1
 n2 = n2 - 1
 places = 0
 psign  = 0
 expont = .false.
 ipower = 0
 precis = 0
 IF (n2 < n1) GO TO 410
 DO  i = n1,n2
   IF (TYPE(i) < 0) THEN
     GO TO   360
   ELSE IF (TYPE(i) == 0) THEN
     GO TO   570
   END IF
   
!     NUMERIC
   
   320 DO  k = 1,9
     IF (CHAR(i) == num(k)) GO TO 340
   END DO
   k  = 0
   340 IF (expont) GO TO 350
   it = it + 1
   IF (it < 16) nt(it) = k
   places = places + 1
   CYCLE
   
!     BUILD IPOWER HERE
   
   350 power  = .true.
   ipower = ipower*10 + k
   IF (ipower > 1000) GO TO 630
   CYCLE
   
!     START EXPONENTS HERE
   
   360 IF (expont) GO TO 380
   expont = .true.
   IF (CHAR(i) /= plus .AND. CHAR(i) /= minus) GO TO 370
   precis = e
   psign = CHAR(i)
   GO TO 390
   370 IF (CHAR(i) /= e .AND. CHAR(i) /= d) GO TO 600
   precis = CHAR(i)
   GO TO 390
   
!     SIGN OF POWER
   
   380 IF (power) GO TO 590
   IF (psign /= 0 .OR.(CHAR(i) /= plus .AND. CHAR(i) /= minus)) GO TO 610
   psign = CHAR(i)
   power = .true.
   390 IF (i == 72) GO TO 530
 END DO
 410 n = n2
 
!     ALL DATA COMPLETE FOR FLOATING POINT NUMBER
!     15 FIGURES WILL BE ACCEPTED ONLY
 
 IF (it <= 15) GO TO 420
 ipower = ipower + it - 15
 it = 15
 420 IF (psign == minus) ipower = -ipower
 ipower = ipower - places
 NUMBER = 0
 IF (it < 7) GO TO 430
 n2 = 7
 GO TO 440
 430 n2 = it
 440 DO  i = 1,n2
   NUMBER = NUMBER*10 + nt(i)
 END DO
 xdoubl = DBLE(FLOAT(NUMBER))
 IF (it <= 7) GO TO 470
 NUMBER = 0
 n2 = it - 7
 DO  i = 1,n2
   it = i + 7
   NUMBER = NUMBER*10 + nt(it)
 END DO
 xdoubl = xdoubl*10.0D0**n2 + DBLE(FLOAT(NUMBER))
 470 IF (lminus) xdoubl = -xdoubl
 
!     POWER HAS TO BE WITHIN RANGE OF MACHINE
 
 ichek = ipower + it
 IF (xdoubl == 0.0D0) GO TO 490
 IF (ichek < lowpw+1 .OR. ichek > highpw-1 .OR.  &
     ipower < lowpw+1 .OR. ipower > highpw-1) GO TO 640
 xdoubl = xdoubl*10.0D0**ipower
 490 IF (precis == d) GO TO 500
 fl1  = xdoubl
 iout = iout + 2
 out(iout-1) =-2
 out(iout  ) = int1
 GO TO 100
 500 iout = iout + 3
 out(iout-2) =-4
 out(iout-1) = ndoubl(1)
 out(iout  ) = ndoubl(2)
 GO TO 100
 
 
!     PREPARE TO RETURN
 
 510 IF (.NOT. delim) GO TO 520
 out(iout+1) = 0
 RETURN
 520 out(iout+1) = a67777
 RETURN
 
!     ERRORS
 
 530 WRITE  (f6,540) ufm
 540 FORMAT (a23,' 300 *** INVALID DATA COLUMN 72')
 GO TO  680
 550 WRITE  (f6,560) ufm
 560 FORMAT (a23,' 300 *** INTEGER DATA OUT OF MACHINE RANGE')
 GO TO  680
 570 WRITE  (f6,580) ufm,n1
 580 FORMAT (a23,' 300 *** INVALID CHARACTER FOLLOWING INTEGER IN ',  &
     'COLUMN',i3)
 GO TO  680
 590 CONTINUE
 600 CONTINUE
 610 WRITE  (f6,620) ufm,i
 620 FORMAT (a23,' 300 *** DATA ERROR-UNANTICIPATED CHARACTER IN ',  &
     'COLUMN',i3)
 GO TO  680
 630 CONTINUE
 640 WRITE  (f6,650) ufm
 650 FORMAT (a23,' 300 *** DATA ERROR - MISSING DELIMITER OR REAL ',  &
     'POWER OUT OF MACHINE RANGE')
 GO TO  680
 660 WRITE  (f6,670) ufm
 670 FORMAT (a23,' 300 *** ROUTINE XRCARD FINDS OUTPUT BUFFER TOO ',  &
     'SMALL TO PROCESS CARD COMPLETELY')
 680 nogo = .true.
 WRITE  (f6,690) CHAR
 690 FORMAT (/5X,1H',72A1,1H')
 out(1) = 0
 
 RETURN
END SUBROUTINE yrcard
