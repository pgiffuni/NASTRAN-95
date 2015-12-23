SUBROUTINE rcard (out,FMT,nflag,in)
!DIR$ INTEGER=64
 
!     CDIR$ IS CRAY COMPILER DIRECTIVE. 64 BIT INTEGER IS USED LOCALLY
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(OUT)                     :: out(1)
 INTEGER, INTENT(OUT)                     :: FMT(1)
 INTEGER, INTENT(OUT)                     :: nflag
 INTEGER, INTENT(IN)                      :: in(1)
 EXTERNAL         lshift,rshift,complf
 LOGICAL :: pass,seqgp,deciml,lminus,expont,DOUBLE,blkon,nogo
 REAL    :: fl1
 INTEGER :: field, f6
 DOUBLE PRECISION :: xdoubl
 DIMENSION        bcd(16),val(16),num(10), TYPE(16),  &
                  nt(16),ndoubl(2),line(20),chars(7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /lhpwx /  lowpw,highpw
 COMMON /system/  ibufsz,f6,nogo,dum1(7),npages,nlines,dum2(10), lsystm
 EQUIVALENCE      (fl1     ,int1 ), (xdoubl,ndoubl(1)),  &
     (num(10) ,zero ), (chars(1),BLANK  ),  &
     (chars(2),star ), (chars(3),plus   ),  &
     (chars(4),minus), (chars(5),period ), (chars(6),e    ), (chars(7),d      )
 DATA    pass  /.false. /, blanks/ 4H     /, stars / 4H**** /,  &
     line  / 20*4H           /         , BLANK / 1H     /,  &
     star  / 1H*    /, plus  / 1H+    /, minus / 1H-    /,  &
     period/ 1H.    /, e     / 1HE    /, d     / 1HD    /,  &
     seq   / 3HSEQ  /, p     / 4HP    /, izero / 0      /
 DATA    num   / 1H1,1H2 ,1H3,1H4, 1H5,1H6,1H7,1H8,1H9,1H0  /
 
 IF (pass) GO TO 40
 pass   = .true.
 a67777 = complf(0)
 a67777 = rshift(lshift(a67777,1),1)
 
 DO  i = 1,10
   num(i) = khrfn1(izero,4,num(i),1)
 END DO
 
 DO  i = 1,7
   chars(i) = khrfn1(izero,4,chars(i),1)
 END DO
 seq = khrfn3(izero,seq,1,0)
 p   = khrfn3(izero,p  ,0,0)
 
 40 field = 0
 nwords= 2
 n_8_or_16 = 8
 word  = 0
 iout  = 0
 ifmt  = 0
 seqgp = .false.
 50 IF (word == 18) GO TO 680
 
!     OPERATE ON 1 FIELD  (2 OR 4 WORDS), GET FIRST NON-BLANK CHARACTER.
 
 field  = field + 1
 deciml = .false.
 lminus = .false.
 expont = .false.
 DOUBLE = .false.
 blkon  = .false.
 places = 0
 it     = 0
 SIGN   = 0
 power  = 0
 
!     READ 8 OR 16 CHARACTERS OF ONE FIELD
 
!     TYPE AS 0 = BLANK, -1 = BCD, +1 = INTEGER
 
 n = 0
 word1 = word + 1
 word  = word + nwords
 DO  i = word1,word
   DO  j = 1,4
     n = n + 1
     charac = khrfn1(izero,4,in(i),j)
     IF (charac == BLANK) GO TO 70
     IF (charac == zero ) GO TO 80
     DO  k = 1,9
       IF (charac == num(k)) GO TO 90
     END DO
     TYPE(n)= -1
     val(n) = charac
     GO TO 100
     70 TYPE(n)= 0
     val(n) = BLANK
     GO TO 100
     80 k = 0
     90 TYPE(n)= 1
     val(n) = k
     100 bcd(n) = charac
   END DO
 END DO
 
!     BCD, INTEGER TRANSFER ON FIRST NON-BLANK CHARACTER
 
 IF (.NOT.seqgp) GO TO 120
 SELECT CASE ( field )
   CASE (    1)
     GO TO 120
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 690
   CASE (    4)
     GO TO 120
   CASE (    5)
     GO TO 690
   CASE (    6)
     GO TO 120
   CASE (    7)
     GO TO 690
   CASE (    8)
     GO TO 120
   CASE (    9)
     GO TO 690
 END SELECT
 
 120 DO  n = 1,n8or16
   IF (TYPE(n) < 0) THEN
     GO TO   150
   ELSE IF (TYPE(n) == 0) THEN
     GO TO   130
   ELSE
     GO TO   320
   END IF
   130 CONTINUE
 END DO
 
!     ALL BLANK FIELD IF FALL HERE
 
 IF (field == 1) GO TO 160
 140 iout = iout + 1
 out(iout) = 0
 ifmt = ifmt + 1
 FMT(ifmt) = 0
 GO TO 50
 
!     **********************************************
 
!     ALPHA HIT FOR FIRST CHARACTER
 
 150 IF (field == 1 .AND. val(n) == star) GO TO 270
 IF (val(n) == plus  ) GO TO 290
 IF (val(n) == period) GO TO 300
 IF (val(n) == minus ) GO TO 310
 
!     PLAIN ALPHA BCD FIELD NOW ASSUMED.
!     CHECKING FOR DOULBE-FIELD * IF WE ARE IN FIELD 1
 
 IF (field /= 1) GO TO 160
 IF (bcd(8) /= star .OR. TYPE(8) /= -1) GO TO 160
 nwords = 4
 n8or16 = 16
 
!     REMOVE STAR BEFORE PUTTING 2 BCD WORDS INTO OUT
 
 bcd(8) = BLANK
 160 iout   = iout + 2
 IF (TYPE(1)) 170,180,170
 170 IF (nwords == 4 .AND. field == 1) GO TO 180
 n = word - nwords
 out(iout-1) = in(n+1)
 out(iout  ) = in(n+2)
 GO TO 260
 
!     CHARACTER N WAS FIRST NON-BLANK CHARACTER
 
 180 MAX = n_8_or_16  -  n  +  1
 DO  i = 1,MAX
   bcd(i) = bcd(n)
   n = n + 1
 END DO
 200 IF (MAX >= 8) GO TO 210
 MAX = MAX + 1
 bcd(MAX) = BLANK
 GO TO 200
 210 word1 = 0
 word2 = 0
 DO  i = 1,4
   word1 = khrfn3(bcd(i  ),word1,1,1)
   word2 = khrfn3(bcd(i+4),word2,1,1)
 END DO
 out(iout-1) = word1
 out(iout  ) = word2
 260 ifmt = ifmt + 1
 FMT(ifmt) = 3
 IF (field /= 1) GO TO 50
 IF (khrfn3(izero,out(iout-1),1,0) == seq .AND.  &
     khrfn3(izero,out(iout  ),0,0) == p) seqgp = .true.
 GO TO 50
 
!     **********************************************
 
!     FIRST CHARACTER  ON CARD IS A STAR
 
 270 nwords = 4
 n_8_or_16 = 16
 280 iout = iout + 2
 out(iout-1) = 0
 out(iout  ) = 0
 ifmt = ifmt + 1
 FMT(ifmt) = 3
 GO TO 50
 
!     **********************************************
 
!     FIRST CHARACTER IN FIELD IS A PLUS
 
 290 IF (field == 1) GO TO 280
 
!     IGNORING PLUS SIGN AND NOW ASSUMING FIELD IS NUMBERIC
 
 GO TO 340
 
!     **********************************************
 
!     FIRST CHARACTER IN FIELD IS A PERIOD
 
 300 deciml = .true.
 places = 0
 GO TO 340
 
!     **********************************************
 
!     FIRST CHARACTER IN FIELD IS A MINUS
 
 310 lminus = .true.
 GO TO 340
 
!     **********************************************
 
!     0 TO 9 NUMERIC HIT
 
 320 IF (val(n) == 0.0) THEN
   GO TO   340
 END IF
 
!     NON-ZERO NUMBER.  SAVING IT NOW IN TABLE NT
 
 330 nt(1) = val(n)
 it = 1
 340 IF (n == n_8_or_16) GO TO 380
 
!     PROCESS REMAINING DIGITS
 
 nnn = n + 1
 DO  n = nnn,n8or16
   IF ((TYPE(n) == 0 .OR. val(n) == zero) .AND. it == 0 .AND.  &
       .NOT.deciml) CYCLE
   IF (TYPE(n) > 0) THEN
     GO TO   360
   END IF
   
!     FALL THRU IMPLIES NON 0 TO 9 CHARACTER
   
   350 IF (val(n) /= period) GO TO 430
   IF (deciml) GO TO 910
   places = 0
   deciml = .true.
   CYCLE
   
!     0 TO 9 CHARACTER HIT. SAVE IT.
   
   360 it = it + 1
   nt(it) = val(n)
   IF (deciml) places = places + 1
 END DO
 
!     NUMERIC WORD COMPLETED
!     IF DECIML IS .FALSE. NUMERIC IS A SIMPLE INTEGER
 
 380 IF (deciml) GO TO 570
 
!     **********************************************
 
!     SIMPLE INTEGER
 
 390 NUMBER = 0
 IF (it == 0) GO TO 410
 DO  i = 1,it
   IF (((a67777-nt(i))/10) < NUMBER) GO TO 890
   NUMBER = NUMBER*10 + nt(i)
 END DO
 410 IF (lminus) NUMBER = - NUMBER
 420 iout = iout + 1
 out(iout) = NUMBER
 ifmt = ifmt + 1
 FMT(ifmt) = 1
 GO TO 50
 
!     **********************************************
 
!     PROBABLE (E, D, +, -) EXPONENT HIT OR BLANK
 
 430 IF (TYPE(n) == 0) THEN
   GO TO   440
 ELSE
   GO TO   460
 END IF
 
!     BLANK HIT THUS ONLY AN EXPONENT OR BLANKS PERMITTED FOR BALANCE
!     OF FIELD
 
 440 IF (n == n_8_or_16) GO TO 450
 n = n + 1
 IF (TYPE(n) < 0) THEN
   GO TO   460
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   440
 ELSE
   GO TO   960
 END IF
 
!     FALL THRU ABOVE LOOP IMPLIES BALANCE OF FIELD WAS BLANKS
 
 450 IF (deciml) GO TO 570
 GO TO 390
 
!     **********************************************
 
!     COMING HERE IMPLIES A NON-BLANK CHARACTER HAS BEEN HIT BEGINNING
!     AN EXPONENT. IT HAS TO BE A (+, -, D, OR E ) FOR NO ERROR
 
 460 IF (val(n) /= plus) GO TO 470
 expont= .true.
 SIGN  = plus
 GO TO 500
 470 IF (val(n) /= e) GO TO 480
 expont= .true.
 GO TO 500
 480 IF (val(n) /= minus) GO TO 490
 expont= .true.
 SIGN  = minus
 GO TO 500
 490 IF (val(n) /= d) GO TO 960
 expont= .true.
 DOUBLE= .true.
 
!     READ INTEGER POWER, WITH OR WITHOUT SIGN
 
 500 IF (n == n_8_or_16) GO TO 950
 n = n + 1
 
 IF (TYPE(n) < 0) THEN
   GO TO   510
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   500
 ELSE
   GO TO   520
 END IF
 510 IF (val(n) /= plus .AND. val(n) /= minus) GO TO 520
 IF (SIGN /= 0) GO TO 940
 SIGN = val(n)
 GO TO 500
 
!     FIRST DIGIT OF INTEGER POWER AT HAND NOW
 
 520 power = 0
 blkon = .false.
 
 530 IF (TYPE(n) > 0) THEN
   GO TO   540
 ELSE
   GO TO   930
 END IF
 540 power = power*10 + val(n)
 
!     GET ANY MORE DIGITS IF PRESENT
 
 550 IF (n == n_8_or_16) GO TO 570
 n = n + 1
 IF (blkon) IF (TYPE(n)) 980,550,980
 IF (TYPE(n) == 0) THEN
   GO TO   560
 ELSE
   GO TO   530
 END IF
 
!     BLANK HIT, BALANCE OF FIELD MUST BE BLANKS
 
 560 blkon = .true.
 GO TO 550
 
!     **********************************************
 
!     SINGLE OR DOUBLE PRECISION FLOATING POINT NUMBER
!     COMPLETE AND OUTPUT IT.
 
!     15 SIGNIFICANT FIGURES POSSIBLE ON INPUT
!     CONSIDERED SINGLE PRECISION UNLESS D EXPONENT IS PRESENT
 
 570 IF (SIGN == minus) power = -power
 power = power - places
 
 NUMBER = 0
 IF (it == 0) THEN
   GO TO   620
 END IF
 580 IF (it < 7) GO TO 590
 n = 7
 GO TO 600
 590 n = it
 600 DO  i = 1,n
   NUMBER = NUMBER*10 + nt(i)
 END DO
 620 xdoubl = DBLE(FLOAT(NUMBER))
 IF (it <= 7) GO TO 640
 NUMBER = 0
 DO  i = 8,it
   NUMBER = NUMBER*10 + nt(i)
 END DO
 xdoubl = xdoubl*10.0D0**(it-7) + DBLE(FLOAT(NUMBER))
 640 IF (lminus) xdoubl = -xdoubl
 
!     CHECK FOR POWER IN RANGE OF MACHINE
 
 ichek = power + it
 IF (xdoubl == 0.0D0) GO TO 660
 IF (ichek < lowpw+1 .OR. ichek > highpw-1 .OR.  &
     power < lowpw+1 .OR. power > highpw-1) GO TO 860
 
 xdoubl = xdoubl*10.0D0**power
 660 ifmt = ifmt + 1
 IF (DOUBLE) GO TO 670
 fl1  = xdoubl
 iout = iout + 1
 out(iout) = int1
 FMT(ifmt) = 2
 GO TO 50
 670 iout = iout + 2
 out(iout-1) = ndoubl(1)
 out(iout  ) = ndoubl(2)
 FMT(ifmt) = 4
 GO TO 50
 680 nflag = iout
 FMT(ifmt+1) = -1
 RETURN
 
!     **********************************************
 
!     FIRST CHARACTER OF FIELD 3, 5, 7,  OR 9 ON SEQGP CARD ENCOUNTERED.
 
!     IT HAS TO BE A 1 TO 9 FOR NO ERROR
 
 690 DO  n = 1,n8or16
   IF (TYPE(n) < 0) THEN
     GO TO  1000
   ELSE IF (TYPE(n) == 0) THEN
     GO TO   700
   ELSE
     GO TO   710
   END IF
   700 CONTINUE
 END DO
 GO TO 140
 
!     STORE NUMBER IN NT
 
 710 npoint = 0
 720 it = it + 1
 nt(it) = val(n)
 730 IF (n == n_8_or_16) GO TO 800
 n = n + 1
 
!     GET NEXT CHARACTER
 
 IF (npoint > 0 .AND. .NOT.deciml .AND. .NOT.blkon) GO TO 790
 IF (deciml) GO TO 770
 IF (blkon ) GO TO 750
 IF (TYPE(n) > 0) THEN
   GO TO   720
 END IF
 740 IF (val(n)  == period) GO TO 760
 750 IF (TYPE(n) /=      0) GO TO 1020
 blkon = .true.
 GO TO 730
 
 760 deciml = .true.
 npoint = npoint + 1
 GO TO 730
 
 770 IF (TYPE(n) > 0) THEN
   GO TO   780
 ELSE
   GO TO  1020
 END IF
 
 780 deciml = .false.
 GO TO 720
 
 790 IF (val(n) == period .AND. TYPE(n) < 0) GO TO 760
 GO TO 750
 
!     READY TO COMPUTE INTEGER VALUE OF SPECIAL SEQGP INTEGER
 
 800 npoint = 3 - npoint
 IF (npoint < 0) THEN
   GO TO  1010
 ELSE IF (npoint == 0) THEN
   GO TO   830
 END IF
 810 DO  k = 1,npoint
   it = it + 1
   nt(it) = 0
 END DO
 
!     COMPUTE NUMBER
 
 830 NUMBER = 0
 IF (it == 0) THEN
   GO TO   420
 END IF
 840 DO  k = 1,it
   IF (((a67777-nt(k))/10) < NUMBER) GO TO 1040
   NUMBER = NUMBER*10 + nt(k)
 END DO
 GO TO 420
 
 
 860 WRITE  (f6,870) ufm
 870 FORMAT (a23,' 300, DATA ERROR IN FIELD UNDERLINED.')
 WRITE  (f6,880)
 880 FORMAT (10X,42HFLOATING point NUMBER out of machine range)
 GO TO  1060
 890 WRITE  (f6,870) ufm
 WRITE  (f6,900)
 900 FORMAT (10X,38HINTEGER magnitude out of machine range)
 GO TO  1060
 910 WRITE  (f6,870) ufm
 WRITE  (f6,920)
 920 FORMAT (10X,22HDATA NOT recognizeable)
 GO TO  1060
 930 CONTINUE
 940 CONTINUE
 950 CONTINUE
 960 WRITE  (f6,870) ufm
 WRITE  (f6,970)
 970 FORMAT (10X,26HPOSSIBLE error in exponent)
 GO TO  1060
 980 WRITE  (f6,870) ufm
 WRITE  (f6,990)
 990 FORMAT (10X,23HPOSSIBLE imbedded BLANK)
 GO TO  1060
 1000 CONTINUE
 1010 CONTINUE
 1020 WRITE  (f6,870) ufm
 WRITE  (f6,1030)
 1030 FORMAT (10X,30HINCORRECT dewey decimal NUMBER)
 GO TO  1060
 1040 WRITE  (f6,870) ufm
 WRITE  (f6,1050)
 1050 FORMAT (10X,49HINTERNAL conversion of dewey decimal is too large)
 1060 word = (field-1)*nwords + 2
 ASSIGN 1090 TO iretrn
 word2 = stars
 1070 line(word  ) = word2
 line(word-1) = word2
 IF (nwords == 2  .OR.  field == 1) GO TO 1080
 line(word-2) = word2
 line(word-3) = word2
 1080 GO TO iretrn,(1090,1150)
 1090 IF (nwords == 4) GO TO 1110
 WRITE  (f6,1100)
 1100 FORMAT (10X,80H.   1  ..   2  ..   3  ..   4  ..   5  ..   6  ..   &
                        &7  ..   8  ..   9  ..  10  .)
 GO TO 1130
 1110 WRITE  (f6,1120)
 1120 FORMAT (10X,80H.   1  ..   2  AND  3  ..   4  AND  5  ..   6  AND  &
     7  ..   8  AND  9  ..  10  .)
 1130 WRITE  (f6,1140) (in(i),i=1,20),line
 1140 FORMAT (10X,20A4)
 ASSIGN 1150 TO iretrn
 word2 = blanks
 GO TO 1070
 1150 iout = iout + 1
 nlines = nlines + 7
 out(iout) = 0
 ifmt = ifmt + 1
 FMT(ifmt) = -1
 nogo = .true.
 GO TO 50
END SUBROUTINE rcard
