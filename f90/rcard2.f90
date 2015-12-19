SUBROUTINE rcard2 (out,FMT,nflag,in)
!DIR$ INTEGER=64
 
!     CDIR$ IS CRAY COMPILE DIRECTIVE. 64-BIT INTEGER IS USED LOCALLY
 
!     THIS ROUTINE IS MUCH MORE EFFICIENT THAN THE OLD ROUTINE RCARD
!     IT CAN SAFELY REPLACE THE OLD RCARD ROUTINE
!     WRITTEN BY G.CHAN/UNISYS            10/1987
!     REVISED, 8/1989, IMPROVED EFFICIENCY BY REDUCING CHARACTER
!     OPERATIONS (VERY IMPORTANT FOR CDC MACHINE)
!     LAST REVISED, 8/1991, SETTING UP REAL NO. UPPER AND LOWER BOUNDS
!     FOR VARIOUS MACHINES
 
!     RCARD2 ASSUMES ALL INPUT FIELDS IN 'IN' ARE LEFT-ADJUSTED.
 
 
 INTEGER, INTENT(OUT)                     :: out(1)
 INTEGER, INTENT(OUT)                     :: FMT(1)
 INTEGER, INTENT(OUT)                     :: nflag
 INTEGER, INTENT(IN OUT)                  :: in(20)
 EXTERNAL         lshift   ,rshift   ,complf
 LOGICAL :: seqgp    ,deciml   ,minus    ,nogo     , expont   ,DOUBLE   ,blkon
 INTEGER :: TYPE(16) , nt(16)   ,outx(100),idoubl(2),value(16),  &
     num1(9)  ,chr1(16) ,a1(80),  field
 REAL :: fpt
 DOUBLE PRECISION :: ddoubl
 CHARACTER (LEN=1) :: blankc   ,starc    ,dotc     ,plusc    ,  &
     minusc   ,dc       ,ec       ,zeroc    , khr1(16) ,k1(80)
 CHARACTER (LEN=4) :: in4(40)  ,c4(1)    ,chr4(4)  ,out4(100)
 CHARACTER (LEN=5) :: d5       ,seqgp5   ,seqep5
 CHARACTER (LEN=100) :: tmp100   ,out100(4)
 CHARACTER (LEN=80) :: e80
 CHARACTER (LEN=23) :: ufm
 CHARACTER (LEN=9) :: num9
 COMMON /xmssg /  ufm
 COMMON /lhpwx /  lowpw    ,highpw
 COMMON /system/  bufsz    ,nout     ,nogo     ,dum1(8)  , nlines
 COMMON /xechox/  dum2(4)  ,xsort2
 EQUIVALENCE      (chr11,chr1(1)),   (k1(1),in4(1),d5,e80),  &
     (fpt,intgr),       (khr1(1),chr4(1)), (ddoubl,idoubl(1)),(out4(1),out100(1))
 DATA    blankc,  starc,   plusc,  minusc,  dotc,  ec,   dc    /  &
     ' ',     '*',     '+',    '-',     '.',   'E',  'D'   /
 DATA    BLANK,   stars,   seqgp5, seqep5,  zeroc, num9        /  &
     4H    ,  4H====,  'SEQGP','SEQEP', '0',   '123456789' /
 DATA    plus1 /  0    /
 
 IF (plus1 /= 0) GO TO 10
 CALL k2b (blankc,blank1,1)
 CALL k2b (starc ,star1 ,1)
 CALL k2b (plusc ,plus1 ,1)
 CALL k2b (minusc,minus1,1)
 CALL k2b (dotc  ,dot1  ,1)
 CALL k2b (ec    ,e1    ,1)
 CALL k2b (dc    ,d1    ,1)
 CALL k2b (zeroc ,zero1 ,1)
 CALL k2b (num9  ,num1  ,9)
 10 CONTINUE
 
 CALL bcdkh8 (in,e80)
 CALL k2b (e80,a1,80)
 GO TO 30
 
 
 ENTRY rcard3 (out,FMT,nflag,c4)
!     ===============================
 
!     IN RCARD2, 'IN' IS 4-BYTE BCD  AND 'OUT' IS 4-BYTE BCD
!     IN RCARD3, 'C4' IS CHARACTER*4 AND 'OUT' IS 4-BYTE BCD
!     'IN' AND 'C4' ARE INPUT, AND 'OUT' IS OUTPUT
 
 DO  i = 1,20
   in4(i) = c4(i)
 END DO
 CALL k2b (c4,a1,80)
 
 30 field  = 0
 iout   = 0
 ifmt   = 0
 word   = 0
 nwords = 2
 seqgp  = .false.
 a67777 = rshift(lshift(complf(0),1),1)/10 - 10
 n8or16 = 8
 DO  i = 1,100
   outx(i) = BLANK
 END DO
 
!     PROCESS ONE FIELD (2 OR 4 WORDS) AT A TIME,
!     GET FIRST NON-BLANK CHARATER
 
 50 IF (word == 18) GO TO 860
 field  = field + 1
 deciml =.false.
 minus  =.false.
 expont =.false.
 DOUBLE =.false.
 blkon  =.false.
 sign1  = blank1
 places = 0
 it     = 0
 power  = 0
 
!     READ 8 OR 16 CHARATERS OF ONE FIELD
!     FOR EACH CHARACTER, SET TYPE TO
!            0 IF IT IS A BLANK
!           -1 IF IT IS BCD CHARACTER, AND
!           +1 IF IT IS NUMERIC
 
 base = word*4
 word = word + nwords
 DO  n = 1,n8or16
   a1nb = a1(n+base)
   IF (a1nb == blank1) GO TO 70
   IF (a1nb == zero1 ) GO TO 80
   DO  k = 1,9
     IF (a1nb == num1(k)) GO TO 90
   END DO
   TYPE(n) = -1
   GO TO 100
   70 TYPE(n) = 0
   GO TO 100
   80 k = 0
   90 TYPE(n) = 1
   value(n)= k
   100 chr1(n) = a1nb
   khr1(n) = k1(n+base)
 END DO
 
 IF (seqgp) THEN
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
 END IF
 
!     BRANCH ON BCD, BLANK, OR NUMERIC
 
 120 IF (TYPE(1)) 150,  130,    320
!                  BCD BLANK NUMERIC
 
!     A BLANK FIELD -
!     ===============
 
 130 IF (field == 1) GO TO 180
 140 iout       = iout + 1
 outx(iout) = 0
 ifmt       = ifmt + 1
 FMT(ifmt)  = 0
 GO TO 50
 
!     BCD FIELD -
!     ===========
 
!     FIRST NON-BLANK CHARATER IS ALPHA, STAR, DOT, PLUS, OR MINUS
 
 150 IF (field == 1 .AND. chr11 == star1) GO TO 270
 IF (chr11 == plus1 ) GO TO 290
 IF (chr11 == dot1  ) GO TO 300
 IF (chr11 == minus1) GO TO 310
 
!     TRUE ALPHA BCD-CHARACTER FIELD
 
!     CHECKING FOR DOULBE-FIELD ASTERISK (*) IF WE ARE IN FIELD 1
!     SET DOUBLE FLAGS N8OR16, NWORDS, AND REMOVE THE ASTERISK
 
 IF (field /= 1) GO TO 180
 j = 8
 DO  i = 2,8
   IF (chr1(j) == star1 .AND. TYPE(j) == -1) GO TO 170
   j = j - 1
 END DO
 GO TO 180
 170 nwords = 4
 n8or16 = 16
 chr1(j) = blank1
 khr1(j) = blankc
 
 180 iout = iout + 2
 IF (TYPE(1)) 190,200,190
 190 IF (nwords == 4 .AND. field == 1) GO TO 200
 n = word - nwords
 out4(iout-1) = in4(n+1)
 out4(iout  ) = in4(n+2)
 GO TO 260
 
 200 out4(iout-1) = chr4(1)
 out4(iout  ) = chr4(2)
 260 ifmt = ifmt + 1
 FMT(ifmt) = 3
 
!     IF FIRST FIELD IS SEQGP OR SEQEP, SET SEQGP FLAG TO TRUE
 
 IF (field == 1 .AND. (d5 == seqgp5 .OR. d5 == seqep5)) seqgp = .true.
 GO TO 50
 
!     FIRST CHARATER ON CARD IS AN ASTERISK (*)
 
 270 nwords = 4
 n8or16 = 16
 280 iout = iout + 2
 outx(iout-1) = 0
 outx(iout  ) = 0
 ifmt = ifmt + 1
 FMT(ifmt) = 3
 GO TO 50
 
!     FIRST CHARATER IN FIELD IS A PLUS (+)
!     IGNOR IT AND ASSUMING REMAINING FIELD IS NUMBERIC
 
 290 IF (field-1 == 0.0) THEN
   GO TO   280
 ELSE
   GO TO   340
 END IF
 
!     FIRST CHARATER IN FIELD IS A DOT (.)
 
 300 deciml = .true.
 places = 0
 GO TO 340
 
!     FIRST CHARATER IN FIELD IS A MINUS (-)
 
 310 minus = .true.
 GO TO 340
 
!     NUMERIC -  0 TO 9
!     =================
 
 320 IF (value(1) == 0) GO TO 340
 nt(1) = value(1)
 it = 1
 
!     PROCESS REMAINING DIGITS
 
 340 DO  n = 2,n8or16
   IF (TYPE(n) > 0) GO TO 360
   
!     A NON-NUMERIC CHARACTER ENCOUNTERED
   
   IF (chr1(n) /= dot1) GO TO 430
   IF (deciml) GO TO 950
   places = 0
   deciml = .true.
   CYCLE
   
!     A NUMERIC CHARACTER, 0 TO 9, SAVE IT IN NT
   
   360 it = it + 1
   nt(it) = value(n)
   IF (deciml) places = places + 1
 END DO
 
!     IF DECIML IS .FALSE. NUMERIC IS AN INTEGER
 
 IF (deciml) GO TO 570
 
!     INTEGER FOUND.  NASTRAN INTEGER LIMIT = 10*A67777
 
 390 NUMBER = 0
 IF (it == 0) GO TO 410
 DO  i = 1,it
   IF (NUMBER > a67777) GO TO 930
   NUMBER = NUMBER*10 + nt(i)
 END DO
 410 IF (minus) NUMBER = - NUMBER
 420 iout = iout + 1
 outx(iout) = NUMBER
 ifmt = ifmt + 1
 FMT(ifmt) = 1
 GO TO 50
 
!     PROBABLY WE JUST ENCOUNTERED (E, D, +, -) EXPONENT, OR BLANK
 
 430 IF (TYPE(n) == 0) THEN
   GO TO   440
 ELSE
   GO TO   460
 END IF
 
!     IT IS A BLANK
!     THUS ONLY AN EXPONENT OR BLANKS PERMITTED FOR BALANCE OF FIELD
 
 440 IF (n == n8or16) GO TO 450
 n = n + 1
 IF (TYPE(n) < 0) THEN
   GO TO   460
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   440
 ELSE
   GO TO   970
 END IF
 
!     FALL THRU ABOVE LOOP IMPLIES BALANCE OF FIELD WAS BLANKS
 
 450 IF (deciml) GO TO 570
 GO TO 390
 
!     A NON-BLANK CHARACTER -
!     IT HAS TO BE A (+, -, D, OR E ) OF THE EXPONENT STRING
 
 460 expont = .true.
 IF (chr1(n) /= plus1) GO TO 470
 sign1  = plus1
 GO TO 500
 470 IF (chr1(n) /= e1) GO TO 480
 GO TO 500
 480 IF (chr1(n) /= minus1) GO TO 490
 sign1  = minus1
 GO TO 500
 490 IF (chr1(n) /= d1) GO TO 970
 DOUBLE = .true.
 
!     READ INTEGER POWER, WITH OR WITHOUT SIGN
 
 500 IF (n == n8or16) GO TO 970
 n = n + 1
 
 IF (TYPE(n) < 0) THEN
   GO TO   510
 ELSE IF (TYPE(n) == 0) THEN
   GO TO   500
 ELSE
   GO TO   520
 END IF
 510 IF (chr1(n) /= plus1 .AND. chr1(n) /= minus1) GO TO 520
 IF (sign1 /= blank1) GO TO 970
 sign1 = chr1(n)
 GO TO 500
 
!     FIRST DIGIT OF INTEGER POWER AT HAND NOW
 
 520 power = 0
 blkon = .false.
 
 530 IF (TYPE(n) > 0) THEN
   GO TO   540
 ELSE
   GO TO   970
 END IF
 540 power = power*10 + value(n)
 
!     GET ANY MORE DIGITS IF PRESENT
 
 550 IF (n == n8or16) GO TO 570
 n = n + 1
 IF (blkon) IF (TYPE(n)) 990,550,990
 IF (TYPE(n) == 0) THEN
   GO TO   560
 ELSE
   GO TO   530
 END IF
 
!     IS A BLANK.  BALANCE OF FIELD MUST BE BLANKS
 
 560 blkon = .true.
 GO TO 550
 
!     SINGLE OR DOUBLE PRECISION FLOATING POINT NUMBER
!     COMPLETE AND OUTPUT IT
 
!     15 SIGNIFICANT FIGURES POSSIBLE ON INPUT
!     CONSIDERED SINGLE PRECISION UNLESS D EXPONENT IS PRESENT
 
 570 IF (sign1 == minus1) power = -power
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
 620 ddoubl = DBLE(FLOAT(NUMBER))
 IF (it <= 7) GO TO 640
 NUMBER = 0
 DO  i = 8,it
   NUMBER = NUMBER*10 + nt(i)
 END DO
 ddoubl = ddoubl*10.0D0**(it-7) + DBLE(FLOAT(NUMBER))
 640 IF (minus) ddoubl = -ddoubl
 
!     CHECK FOR POWER IN RANGE OF MACHINE
 
 check = power + it
 IF (ddoubl == 0.0D0) GO TO 660
 IF (check < lowpw+1 .OR. check > highpw-1 .OR.  &
     power < lowpw+1 .OR. power > highpw-1) GO TO 900
 
 ddoubl = ddoubl*10.0D0**power
 660 ifmt = ifmt + 1
 IF (DOUBLE) GO TO 670
 fpt  = ddoubl
 iout = iout + 1
 outx(iout)= intgr
 FMT(ifmt) = 2
 GO TO 50
 670 iout = iout + 2
 outx(iout-1) = idoubl(1)
 outx(iout  ) = idoubl(2)
 FMT(ifmt) = 4
 GO TO 50
 
!     FIRST CHARATER OF FIELD 3, 5, 7,  OR 9 ON SEQGP/SEQEP CARD
!     ENCOUNTERED. IT HAS TO BE A 1 TO 9 FOR NO ERROR
 
 690 DO  n = 1,n8or16
   IF (TYPE(n) < 0) THEN
     GO TO  1020
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
 nt(it) = value(n)
 730 IF (n == n8or16) GO TO 800
 n = n + 1
 
!     GET NEXT CHARATER
 
 IF (npoint > 0 .AND. .NOT.deciml .AND. .NOT.blkon) GO TO 790
 IF (deciml)  GO TO 770
 IF (blkon )  GO TO 750
 IF (TYPE(n) > 0) THEN
   GO TO   720
 END IF
 740 IF (chr1(n) == dot1) GO TO 760
 750 IF (TYPE(n) /=    0) GO TO 1020
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
 
 790 IF (chr1(n) == dot1 .AND. TYPE(n) < 0) GO TO 760
 GO TO 750
 
!     READY TO COMPUTE INTEGER VALUE OF SPECIAL SEQGP/SEQEP INTEGER
 
 800 npoint = 3 - npoint
 IF (npoint < 0) THEN
   GO TO  1020
 ELSE IF (npoint == 0) THEN
   GO TO   830
 END IF
 810 DO  k = 1,npoint
   it = it + 1
   nt(it) = 0
 END DO
 
!     COMPUTE NUMBER.  NASTRAN INTEGER LIMIT = 10*A67777
 
 830 NUMBER = 0
 IF (it == 0) THEN
   GO TO   420
 END IF
 840 DO  k = 1,it
   IF (NUMBER > a67777) GO TO 1040
   NUMBER = NUMBER*10 + nt(k)
 END DO
 GO TO 420
 
!     ALL FIELDS PROCESSED
 
 860 nflag = iout
 FMT(ifmt+1) = -1
 
!     CONVERT CHARACTERS TO BCD, AND INSERT NUMERIC VALUES IF
!     APPLICABLE
 
 n = 1
 DO  i = 1,nflag,25
   k = i + 24
   tmp100 = out100(n)
   CALL khrbc1 (tmp100,out(i))
   DO  j = i,k
     IF (outx(j) /= BLANK) out(j)=outx(j)
   END DO
   n = n + 1
 END DO
 RETURN
 
!     ERROR
 
 900 WRITE  (nout,910) ufm
 910 FORMAT (a23,' 300, DATA ERROR IN FIELD UNDERLINED.')
 WRITE  (nout,920)
 920 FORMAT (10X,'FLOATING POINT NUMBER OUT OF MACHINE RANGE')
 WRITE  (nout,925) power,it,check,lowpw,highpw
 925 FORMAT (10X,'POWER,IT,CHECK,LOWPW,HIGHPW =',5I5)
 GO TO  1060
 930 WRITE  (nout,910) ufm
 WRITE  (nout,940)
 940 FORMAT (10X,'INTEGER MAGNITUDE OUT OF MACHINE RANGE')
 GO TO  1060
 950 IF (xsort2 == 2) GO TO 50
 WRITE  (nout,910) ufm
 WRITE  (nout,960)
 960 FORMAT (10X,'DATA NOT RECOGNIZEABLE')
 GO TO  1060
 970 expont = .false.
 IF (xsort2 == 2) GO TO 50
 WRITE  (nout,910) ufm
 WRITE  (nout,980)
 980 FORMAT (10X,'POSSIBLE ERROR IN EXPONENT')
 GO TO  1060
 990 IF (xsort2 == 2) GO TO 50
 WRITE  (nout,910) ufm
 WRITE  (nout,1000)
 1000 FORMAT (10X,'POSSIBLE IMBEDDED BLANK')
 GO TO  1060
 1020 IF (xsort2 == 2) GO TO 50
 WRITE  (nout,910) ufm
 WRITE  (nout,1030)
 1030 FORMAT (10X,'INCORRECT DEWEY DECIMAL NUMBER')
 GO TO  1060
 1040 IF (xsort2 == 2) GO TO 50
 WRITE  (nout,910) ufm
 WRITE  (nout,1050)
 1050 FORMAT (10X,'INTERNAL CONVERSION OF DEWEY DECIMAL IS TOO LARGE')
 1060 DO  j = 1,20
   IF (outx(j) /= stars) outx(j) = BLANK
 END DO
 word = (field-1)*nwords + 2
 k = stars
 1080 outx(word  ) = k
 outx(word-1) = k
 IF (nwords == 2 .OR. field == 1) GO TO 1090
 outx(word-2) = k
 outx(word-3) = k
 1090 IF (k == 0) GO TO 1150
 IF (nwords == 4) GO TO 1110
 WRITE  (nout,1100)
 1100 FORMAT (10X,'---1--- +++2+++ ---3--- +++4+++ ---5--- +++6+++ ',  &
     '---7--- +++8+++ ---9--- +++10+++')
 GO TO  1130
 1110 WRITE  (nout,1120)
 1120 FORMAT (10X,'---1--- +++++2+&+3+++++ -----4-&-5----- +++++6+&',  &
     '+7+++++ -----8-&-9----- +++10+++')
 1130 WRITE  (nout,1140) (in4(i),i=1,20),outx
 1140 FORMAT (10X,20A4)
 nlines = nlines + 7
 k = 0
 GO TO  1080
 1150 iout = iout + 1
 outx(iout) = 0
 ifmt = ifmt + 1
 FMT(ifmt) = -1
 nogo =.true.
 GO TO 50
 
END SUBROUTINE rcard2
