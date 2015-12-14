SUBROUTINE xread (*,bufx)
     
!     THIS ROUTINE MAKES FREE-FIELD INPUT PACKAGE (HANDLED BY FFREAD)
!     COMPLETELY MACHINE INDEPENDENT.
 
!     IF THE XSORT FLAG IN /XECHOX/ IS TURNED ON (XSORT=1), THIS ROUTINE
!     WILL ALSO PREPARES THE NECESSARY GROUND WORK SO THAT THE INPUT
!     CARDS CAN BE SORTED EFFICIENTLY IN XSORT2 ROUTINE. ALL FIELDS IN
!     THE INPUT CARDS ARE ALSO LEFT-ADJUSTED FOR PRINTING.
 
!     WRITTEN BY G.CHAN/UNISYS.   OCT. 1987
!     LAST REVISED, 1/1990, IMPROVED EFFICIENCY BY REDUCING CHARACTER
!     OPERATIONS (VERY IMPORTANT FOR CDC MACHINE)
 
 
 , INTENT(IN)                             :: *
 INTEGER, INTENT(IN OUT)                  :: bufx(20)
 IMPLICIT INTEGER (a-z)
 EXTERNAL         rshift,complf
 LOGICAL :: DOUBLE,bcd2,bcd3,alpha,numric
 INTEGER :: sub(2)
 INTEGER :: card1(80),khr1(43),blank1,dollr1,slash1,star1,  &
     plus1,minus1,zero1,point1,e1,d1,j1
 CHARACTER (LEN=1) :: kard1(80),khrk(43),blankk,equ1
 CHARACTER (LEN=8) :: card8(10),card81,blank8,slash8,end8(3),name8(15)
 CHARACTER (LEN=80) :: card80
 CHARACTER (LEN=43) :: khr43
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /xechox/  dummy,echou,osop(2),xsort,wasff,dum,f3long,large
 COMMON /xsortx/  ibuf(4),table(255)
 COMMON /system/  bufsz,nout,nogo
 COMMON /machin/  mach
 EQUIVALENCE      (kard1(1),card8(1), card80 ,card81  ) ,  &
     (blank1,khr1( 1)) , (khr43 ,khrk( 1)) ,  &
     (zero1 ,khr1( 2)) , (d1    ,khr1(15)) ,  &
     (e1    ,khr1(16)) , (slash1,khr1(38)) ,  &
     (dollr1,khr1(39)) , (star1 ,khr1(40)) ,  &
     (plus1 ,khr1(41)) , (minus1,khr1(42)) , (point1,khr1(43))
 DATA    blank8 , slash8 ,  blank4  , equal4 , sub            /  &
     '    ' , '/   ' ,  4H      , 4H==== , 4HXREA,4HD     /
 DATA    nname /  15 /   ,  name8                             /  &
     'SPC1 ', 'SPCS ', 'TICS '  , 'MPCS ','MPCAX', 'RELES',  &
     'GTRAN','FLUTTER','BDYC '  , 'SPCSD','SPCS1','RANDPS',  &
     'DAREAS','DELAYS', 'DPHASES'                          /
DATA    end8  /  'ENDDATA ','ENDATA ','END DATA'/, derr / -1 /

DATA    khr43 /' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/$*+-.' /
!                      2 4 6 8 1 2 4 6 8 2 2 4 6 8 3 2 4 6 8 4 2
!                              0         0         0         0
DATA    n7, n1, n2,      n3,                       n4, n5,n6 /  &
    44,  1,  2,      11,                       37, 41,43 /

DATA    plus1 , blankk, equ1 / 0, ' ', '=' /

IF (plus1 == 0) CALL k2b (khr43,khr1,43)

!     CALL FFREAD TO READ INPUT CARD
!     IF INPUT IS A COMMENT CARD, SET IBUF(1)=-1, AND RETURN
!     IF INPUT IS IN FREE-FIELD, ALL 10 BULKDATA FIELDS ARE ALREADY
!     LEFT-ADJUSTED, AND WASFF IS SET TO +1 BY FFREAD
!     IF INPUT IS IN FIXED-FIELD, ALL 10 BULKDATA FIELDS MAY NOT BE IN
!     LEFT-ADJUSTED FORMAT, AND WASFF IS SET TO -1 BY FFREAD

CALL ffread (*850,card8)
CALL k2b (card80,card1,80)
IF (card1(1) == dollr1) GO TO 770
ie = 0
IF (xsort == 0 .OR. wasff == 1) GO TO 40

!     LEFT-ADJUSTED THE BULKDATA FIELDS, FIRST 9 FIELDS
!     (FIRST 4 AND A HALF FIELDS IF DOUBLE FIELD CARDS)

ib = 1
l  = 8
IF (card1(1) == plus1 .OR. card1(1) == star1) ib = 9
IF (card1(1) == star1) l = 16
DO  i = ib,72,l
  IF (card1(i) /= blank1) CYCLE
  k  = i
  je = i + l - 1
  DO  j = i,je
    IF (card1(j) == blank1) CYCLE
    card1(k) = card1(j)
    kard1(k) = kard1(j)
    k  = k + 1
  END DO
  IF (k == i) CYCLE
  DO  j = k,je
    kard1(j) = blankk
    card1(j) = blank1
  END DO
END DO

!     CHECK COMMENT CARD WITH DOLLAR SIGN NOT IN COLUMN 1. CONVERT
!     CHARACTER STRING TO BCD STRING, AND RETURN TO CALLER IF IT IS
!     NOT CALLED BY XSORT.

40 ie = ie + 1
IF (card1(ie) == blank1) GO TO 40
IF (card1(ie) == dollr1) GO TO 760
CALL khrbcd (card80,bufx)
IF (xsort == 0) GO TO 780


!     IF THIS ROUTINE IS CALLED BY XSORT, PASS THE FIRST 3 FIELDS TO
!     IBUF ARRAY IN /XSORTX/, IN INTEGER FORMS

!     FIRST BULKDATA FIELD IS ALPHA-NUMERIC, COMPOSED OF TWO 4-CHARACTER
!     WORDS. CHECK WHETHER OR NOT THIS IS A CONTINUATION OR COMMENT CARD
!     IF IT IS NOT, WE CHANGE ALL 8 CHARACTER BYTES INTO THEIR NUMERIC
!     CODE VALUES GIVEN BY TABLE /KHR43/ AND STORE THE VALUE IN IBUF(1)
!     AND IBUF(2)

!     WE SET IBUF(1) AND (2)    IF INPUT CARD IS
!     ----------------------    -------------------
!                -1             A COMMENT CARD
!                -2             A CONTINUATION CARD
!                -3             A DELETE CARD (RANGE IN IBUF(3) AND (4))
!                -3, -4         A DIRTY DELETE CARD
!                -5             A BLANK CARD
!                -9             A ENDDATA CARD
!     AND IBUF(2) AND IBUF(3) ARE NOT SET, EXECPT -3 CASE

!     IF FIELD 2 AND/OR FIELD 3 ARE IN CHARACTERS, WE PUT THE FIRST 6
!     BYTES (OUT OF POSSIBLE 8 CHARACTER-BYTES) INTO IBUF(3) AND/OR
!     IBUFF(4) RESPECTIVELY, IN INTERNAL NUMERIC CODE QUITE SIMILAR TO
!     RADIX-50
!     IF FIELD 2 HAS MORE THAN 7 CHARACTERS, IBUF(4) IS USED TO RECEIVE
!     THE LAST 2 CHARACTERS OF FIELD 2

!     IF FIELD 2 AND/OR FIELD 3 ARE NUMERIC DATA (0-9,+,-,.,E), THEIR
!     ACTUAL INTEGER VALUES ARE STORED IBUF(3) AND/OR IBUF(4).
!     IF THEY ARE F.P. NUMBERS, THEIR EXPONENT VALUES (X100000) ARE
!     CHANGED INTO INTEGERS, AND THEN STORED IN IBUF(3) AND/OR IBUF(4)

!     NOTE - XREAD WILL HANDLE BOTH SINGLE- AND DOUBLE-FIELD BULKDATA
!     INPUT IN FIELDS 2 AND 3, AND MOVED THEM ACCORDINGLY INTO IBUF(3)
!     AND IBUF(4)


!     PRESET TABLE IF THIS IS VERY FIRST CALL TO XREAD ROUTINE
!     TABLE SETTING IS MOVED UP BY ONE IF MACHINE IS CDC (TO AVOID
!     BLANK CHARACTER WHICH IS ZERO FROM ICHAR FUNCTION)

fromy = 0
IF (xsort /= 1) GO TO 80
xsort = 2
cdc   = 0
IF (mach == 4) cdc = 1
DO  i = 1,255
  table(i) = n7
END DO
DO  i = 1,n6
  j = ICHAR(khrk(i)) + cdc
  table(j) = i
END DO
f3long = 0
large  = rshift(complf(0),1)/20

!     CHECK BLANK, ENDDATA, AND CONTINUATION CARDS

80 er = 0
j1 = card1(1)
j  = table(ICHAR(kard1(1))+cdc)
IF (j >= n7) GO TO 810
IF (card81 == blank8 .AND. card8(2) == blank8)  GO TO 90
IF (card81 == end8(1) .OR. card81 == end8(2) .OR.  &
    card81 == end8(3)) GO TO 100
IF (j1 /= plus1 .AND. j1 /= star1) GO TO 120
ibuf(1) = -2
GO TO 110
90 ibuf(1) = -5
GO TO 110
100 ibuf(1) = -9
110 ibuf(2) = ibuf(1)
GO TO 800

!     CHECK ASTERISK IN FIELD 1 (BUT NOT IN COLUMN1 1) AND SET DOUBLE-
!     FIELD FLAG. MERGE EVERY TWO SINGLE FIELDS TO ENSURE CONTINUITY OF
!     DOUBLE FIELD DATA (FIXED FIELD CARDS ONLY)

120 DOUBLE = .false.
IF (wasff == 1) GO TO 180
ie = 8
DO  j = 2,8
  IF (card1(ie) == star1) GO TO 140
  ie = ie - 1
END DO
GO TO 180
140 DOUBLE = .true.
ib = 0
DO  i = 8,71,16
  k = i
  DO  j = 1,16
    l = i + j
    IF (card1(l) == blank1) CYCLE
    k = k + 1
    IF (k == l) CYCLE
    ib = 1
    card1(k) = card1(l)
    kard1(k) = kard1(l)
  END DO
  IF (k == l) CYCLE
  k = k + 1
  DO  j = k,l
    kard1(j) = blankk
    card1(j) = blank1
  END DO
END DO
IF (ie <= 0) CALL mesage (-37,0,sub)
IF (ib == 1) CALL khrbcd (card80,bufx)
card1(ie) = blank1
kard1(ie) = blankk

!     CHECK DELETE CARD
!     SET IBUF(1)=IBUF(2)=-3 IF IT IS PRESENT, AND SET THE DELETE RANGE
!     IN IBUF(3) AND IBUF(4)
!     SET IBUF(1)=-3 AND IBUF(2)=-4 IF TRASH FOUND AFTER SLASH IN
!     FIELD 1
!     NOTE - IF FIELD 3 IS BLANK, IBUF(4) IS -3

180 IF (j1 /= slash1) GO TO 200
DO  l = 1,4
  ibuf(l) = -3
END DO
IF (card81 /= slash8) ibuf(2) = -4
l = 2
GO TO 300

!     TURN BCD2 AND BCD3 FLAGS ON IF THE 2ND AND 3RD INPUT FIELDS ARE
!     NOT NUMERIC RESPECTIVELY
!     IF 2ND FIELD HAS MORE THAN 6 CHARACTERS, REPLACE 3RD FIELD BY THE
!     7TH AND 8TH CHARACTERS OF THE 2ND FIELD
!     (FOR DMI AND DTI CARDS, MERGE 7TH AND 8TH CHARACTERS INTO 3RD
!     FIELD AND TREAT THE ORIG. 3RD FIELD AS A NEW BCD WORD)
!     IF 3RD FIELD HAS MORE THAN 6 CHARACTERS, SET F3LONG FLAG TO 1, AND
!     USER INFORMATION MESSAGE 217A WILL BE PRINTED BY XSORT
!     FIELDS 2 AND 3 SHOULD NOT START WITH A /, $, *
!     IF FIELD2 IS A BCD WORD, FIELD3 PROCESSING ACTUALLY BEGINS IN
!     CARD8(4)

200 bcd2 = .false.
IF (derr == +1) derr = 0
j = table(ICHAR(kard1(9))+cdc)
IF (j >= n7) GO TO 810
numric = (j >= n2 .AND. j <= n3) .OR. j >= n5
IF (numric) GO TO 210
bcd2 = .true.
IF (card1(15) == blank1) GO TO 210

!     SINCE THE NAME IN THE 2ND FIELD OF DMI, DTI, DMIG, DMIAX CARDS
!     ARE NOT UNIQUELY DEFINED FOR SORTING, SPECIAL CODES HERE TO MOVE
!     THE LAST PART OF A LONG NAME (7 OR 8 LETTER NAME) INTO THE 3RD
!     FIELD, AND TREAT THE NEW 3RD FIELD AS BCD WORD. THUS THE ORIGINAL
!     3RD FIELD (THE COLUMN NUMBER, RIGHT ADJUSTED WITH LEADING ZEROS)
!     IS LIMITED TO 4 DIGITS OR LESS.  IF THE NAME IN THE 2ND FIELD IS
!     SHORT (6 LETTERS OR LESS), MERGING OF THE 3RD FIELD IS NOT NEEDED.

IF (card1(1) /= d1       .OR.  card1(3) /= khr1(20) .OR.  &
    (card1(2) /= khr1(24) .AND. card1(2) /= khr1(31))) GO TO 208
bcd3 = .true.
k = 24
IF (DOUBLE) k = 32
IF (card1(k-3) == blank1) GO TO 204
IF (echou == 1) GO TO 202
IF (derr == -1) CALL page
CALL page2 (-2)
IF (DOUBLE) card1(8) = star1
WRITE  (nout,201) card8
201 FORMAT (30X,10A8)
IF (DOUBLE) card1(8) = blank1
202 CALL page2 (-2)
WRITE  (nout,203) ufm
203 FORMAT (a23,', THE 3RD INPUT FIELD OF THE ABOVE CARD IS LIMITED ',  &
    'TO 4 OR LESS DIGITS, WHEN A NAME OF 7 OR MORE', /5X,  &
    'LETTERS IS USED IN THE 2ND FIELD',/)
derr = +1
nogo =  1
204 DO  j = 1,4
  IF (card1(k-4) /= blank1) GO TO 206
  kard1(k-4) = kard1(k-5)
  kard1(k-5) = kard1(k-6)
  kard1(k-6) = kard1(k-7)
  kard1(k-7) = blankk
END DO
206 DO  j = 1,6
  kard1(k) = kard1(k-2)
  k = k-1
END DO
kard1(k  ) = kard1(16)
kard1(k-1) = kard1(15)
kard1( 15) = blankk
GO TO 215

208 kard1(17) = kard1(15)
kard1(18) = kard1(16)
DO  k = 19,24
  kard1(k) = blankk
END DO

210 bcd3 = .false.
k = 17
IF (DOUBLE) k = 25
j = table(ICHAR(kard1(k))+cdc)
alpha = j == n1 .OR. (j > n3 .AND. j < n5)
IF (alpha) bcd3 = .true.
IF (bcd3 ) GO TO 215

!     THE FIRST 3 FIELDS OF THE DMIG OR DMIAX CARDS (NOT THE 1ST HEADER
!     CARD), ARE NOT UNIQUE. MERGE THE 4TH FIELD (1 DIGIT INTEGER) INTO
!     THE 3RD FIELD (INTEGER, 8 DIGITS OR LESS) TO INCLUDE THE COMPONENT
!     FIELD FOR SORTING

IF (card1(1) /= d1 .OR. card1(2) /= khr1(24) .OR.  &
    card1(3) /= khr1(20)  .OR. (card1(4) /= khr1(18) .AND.  &
    card1(4) /= khr1(12))) GO TO 215
IF (card1(1) == khr1(2)) GO TO 215
k = 24
IF (DOUBLE) k = 32
IF (card1(k) /= blank1) GO TO 215
DO  j = 1,7
  k = k - 1
  IF (card1(k) /= blank1) EXIT
END DO
212 kard1(k+1) = kard1(25)
IF (DOUBLE) kard1(k+1) = kard1(41)


!     CHANGE ALL CHARACTERS IN FIRST 3 FIELDS TO INTEGER INTEGER CODES
!     ACCORDING TO THE TABLE ARRANGEMENT IN /KHR43/
!     MAKE SURE THE INTERNAL CODE IS NOT IN NASTRAN INTEGER RANGE (1 TO
!     8 DIGITS), AND WITHIN MACHINE INTEGER WORD LIMIT
!     IN 2ND AND 3RD FIELDS, INTERCHANGE ALPHABETS AND NUMERIC DIGITS
!     SEQUENCE TO AVOID SYSTEM INTEGER OVERFLOW

!     -------------- REMEMBER, FROM HERE DOWN,
!                    CARD1 (1-BYTE ) HOLD ONE CHARACTER, AND
!                    IBUFX (4-BYTES) HOLD AN  INTEGER -----------------
!     WE ALSO HAVE   CARD8 (8-BYTES) HOLDING 8 CHARACTERS,
!              AND   BUFX  (4-BYTES) HOLDING 4 BCD-CHARACTERS


!     MAP OF THE FIRST 3 BULKDATA FIELDS -
!     (INPUT)

!           WORD1 WORD2 WORD3 WORD4 WORD5 WORD6 WORD7 WORD8 WORD9 WORD10
!     BYTE: 1         8 9        16 17       24 25       32 33       40
!          +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
!     SF:  !<-FIELD 1->!<-FIELD 2->!<-FIELD 3->!
!     DF:  !<-FIELD 1->!<------ FIELD 2 ------>!<------ FIELD 3 ------>!


!     MAP OF IBUF -           WORD1 WORD2 WORD3 WORD4
!     (OUTPUT)         BYTE:  1         8 9  12 13 16
!                            +-----+-----+-----+-----+
!     FOR CORE SORT          !<-FIELD 1->!<--->!<--->!
!     PERFORMED IN                        FIELD FIELD
!     XSORT2                                 2     3

215 numric = .false.
l    = 0
ioo  = 100
word = 0
220 word = word + 1
GO TO 260
230 ioo  = n4
IF (.NOT.bcd2) GO TO 280
word = 3
GO TO 260
240 word = 5
IF (DOUBLE) word = 7
IF (.NOT.bcd2 .OR. kard1(15) == blankk) GO TO 250
word = 4
ioo  = 100
bcd3 = .true.
250 IF (.NOT.bcd3) GO TO 280
IF (word /= 4 .AND. kard1(word*4+3) /= blankk .AND. derr /= +1) f3long = 1
260 ie  = word*4
ib  = ie - 3
j   = table(ICHAR(kard1(ib))+cdc)
IF (j >= n7) GO TO 810
IF (MOD(word,2) == 0 .AND. .NOT.numric) GO TO 262
numric = (j >= n2 .AND. j <= n3) .OR. j >= n5
IF (numric) GO TO 280
262 IF (ioo == 100) GO TO 265
ie = ie + 2
k  = j
IF (k > n3) j = k - n3
IF (k <= n3) j = k + 25
265 sum = j
ib  = ib + 1
DO  i = ib,ie
  j   = table(ICHAR(kard1(i))+cdc)
  sum = sum*ioo + j
END DO
IF (ioo == 100) sum = sum + 200000000
ibuf(l+1) = sum
280 l = l + 1
SELECT CASE ( l )
  CASE (    1)
    GO TO 220
  CASE (    2)
    GO TO 230
  CASE (    3)
    GO TO 240
  CASE (    4)
    GO TO 290
END SELECT

!     CHECK INTEGERS ON 2ND AND 3RD FIELDS

290 IF (bcd2 .AND. bcd3) GO TO 500
l  = 2
IF (bcd2) l = 3
300 l  = l + 1
IF (l-4 < 0) THEN
  GO TO   310
ELSE IF (l-4 == 0) THEN
  GO TO   320
ELSE
  GO TO   500
END IF
310 ib = 9
GO TO 330
320 ib = 17
IF (DOUBLE) ib = 25
330 ie = ib + 7
IF (DOUBLE) ie = ib + 15
j1 = card1(ib)
IF (j1 == plus1 .OR. j1 == minus1 .OR. j1 == point1 .OR.  &
    j1 == zero1) GO  TO 340
j  = table(ICHAR(kard1(ib))+cdc)
IF (j >= n2 .AND. j <= n3) GO TO 350

!     IT IS CHARACTER FIELDS, NOTHING ELSE NEEDS TO BE DONE

GO TO 300

!     IT IS NUMERIC

340 ib  = ib + 1
350 sum = 0
fp  = 0
sigx= 1
SIGN= 1
IF (j1 == minus1) SIGN =-1
IF (j1 == point1) fp   = 1
DO  i = ib,ie
  IF (kard1(i) == blankk) GO TO 390
  j   = table(ICHAR(kard1(i))+cdc) - n2
  IF (j < 0 .OR. j > 9) GO TO 360
  IF (fp <= 0 .AND. IABS(sum) < large) sum = sum*10 + SIGN*j
  CYCLE
  
!     A NON-NUMERIC SYMBOL FOUND IN NUMERIC STRING
!     ONLY 'E', 'D', '+', '-', OR '.' ARE ACCEPTABLE HERE
  
  360 j1  = card1(i)
  IF (j1 == point1) GO TO 370
  IF (fp == 0 .OR. ibuf(3) == -3) GO TO 420
  IF (j1 /= e1 .AND. j1 /= d1 .AND. j1 /= plus1 .AND. j1 /= minus1) GO TO 420
  IF (j1 == minus1) sigx = -1
  fp  =-1
  sum = 0
  CYCLE
  370 fp  = 1
  
END DO

!     BEEF UP NUMERIC DATA BY 2,000,000,000 SO THAT THEY WILL BE
!     SORTED BEHIND ALL ALPHABETIC DATA, AND MOVE THE NUMERIC DATA, IN
!     INTEGER FORM (F.P. MAY NOT BE EXACT) INTO IBUF(3) OR IBUF(4)

390 IF (fp < 0.0) THEN
  GO TO   410
END IF
400 ibuf(l) = sum + SIGN*2000000000
GO TO 300
410 ibuf(l) = SIGN*2000000000
IF (sigx > 0 .AND. sum < 9) ibuf(l) = SIGN * (10**(sigx*sum) + 2000000000)
IF (sigx > 0 .AND. sum >= 9) ibuf(l)= 2147000000*SIGN
GO TO 300

!     ERROR IN NUMERIC FIELD

420 IF (ib == 10 .OR. ib == 18) ib = ib - 1
k = 1
IF (echou == 0 .AND. er /= -9) k = 2
CALL page2 (-k)
IF (echou == 0 .AND. er /= -9) WRITE (nout,430) card80
430 FORMAT (1H ,29X,a80)
k = 2
IF (.NOT.DOUBLE) GO TO 440
k = 4
IF (l /= 4) word = word + 1
440 IF (l == 4) word = word + 2
WRITE  (nout,450) (blank4,i=1,word),(equal4,i=1,k)
450 FORMAT (7X,'*** ERROR -',24A4)
nogo = 1
er   =-9
GO TO 500

!     BOTH FIELDS 2 AND 3 (OF BULK DATA CARD) DONE.


!     FOR MOST BULK DATA CARDS, EXCEPT THE ONES IN NAME8, THE FIRST
!     3 FIELDS, IN INTERNAL CODES AND SAVED IN THE IBUF 4-WORD ARRAY,
!     ARE SUFFICIENT FOR ALPHA-NUMERIC SORT (BY XSORT2)

!     THOSE SPECIAL ONES IN NAME8 ADDITIONAL FIELDS FOR SORTING

500 DO  TYPE = 1,nname
  IF (card81 == name8(TYPE)) GO TO  &
      (520,   520,   520,   520,   600,   520,    520,   520,  &
      520,   560,   570,   580,   560,   560,    560),  TYPE
  
!    1   SPC1   SPCS   TICS   MPCS  MPCAX  RELES   GTRAN  FLUTTER
!    2   BDYC  SPCSD   SPCS1 RANDPS DELAYS DAREAS  DPHASES
  
END DO
GO TO 700

!     SPC1,SPCS,TICS,MPCS,RELES,GTRAN,FLUTTER,BDYC CARDS -
!     ADD 4TH INTEGER FIELD TO IBUF ARRAY

520 ibuf(2) = ibuf(3)
ibuf(3) = ibuf(4)
530 sum = 0
DO  i = 25,32
  j1  = card1(i)
  IF (j1 == blank1) EXIT
  j   = table(ICHAR(kard1(i))+cdc) - n2
  IF (j >= 0 .AND. j <= 9) sum = sum*10 + j
END DO
550 ibuf(4) = sum
IF (TYPE == 12) GO TO 590
GO TO 700

!     DAREAS,DELAYS,DPHASES,SPCSD CARDS -
!     ADD ONE TO IBUF(1), THUS CREATE DARF,DELB,DPHB,OR SPCT IN
!     IBUF(1), THEN ADD 4TH INTEGER FIELD INTO IBUF ARRAY

560 ibuf(1) = ibuf(1) + 1
GO TO 520

!     SPCS1 CARD -
!     ADD TWO TO IBUF(1), THUS CREATE SPCU IN IBUF(1), THEN ADD
!     4TH INTEGER FIELD INTO IBUF ARRAY

570 ibuf(1) = ibuf(1) + 2
GO TO 520

!     RANDPS -
!     MERGE FIELDS 3 AND 4 IF SUBCASE NUMBERS ARE NOT TOO BIG

580 IF (ibuf(4) >= 10000 .OR. bufx(8) /= blank4) GO TO 700
ioooo = ibuf(4)*10000
GO TO 530
590 ibuf(4) = ibuf(4) + ioooo
GO TO 700

!     MPCAX -
!     MOVE THE 6TH FIELD INTO IBUF(4)

600 j = 41
DO  i = 25,32
  card1(i) = card1(j)
  kard1(i) = kard1(j)
  j = j+1
END DO
GO TO 530

!     CHECK NUMERIC ERROR IN 4TH TO 9TH FIELDS IF NO ERROR IN FIRST
!     3 FIELDS (NEW BULK DATA CARDS ONLY)

700 IF (fromy == 1 .OR. er == -9) GO TO 800
word = 5
IF (DOUBLE) word = 7
710 word = word + 2
IF (DOUBLE) word = word + 2
IF (word >= 19) GO TO 800
ib = word*4 - 3
j  = table(ICHAR(kard1(ib))+cdc)
IF (j >= n7) GO TO 710
alpha = j == n1 .OR. (j > n3 .AND. j < n5)
IF (alpha) GO TO 710
ie = ib + 7
IF (DOUBLE) ie = ib + 15
l  = ib + 1
DO  i = l,ie
  j1 = card1(i)
  IF (j1 == blank1) GO TO 710
  j  = table(ICHAR(kard1(i))+cdc)
  numric = (j >= n2 .AND. j <= n3) .OR. (j >= n5 .AND. j <= n6)
  IF (numric .OR. j == 15 .OR. j == 16) CYCLE
!                           D            E
  k = 1
  IF (echou == 0 .AND. er /= -9) k = 2
  CALL page2 (-k)
  IF (echou == 0 .AND. er /= -9) WRITE (nout,430) card80
  word = word + 2
  k = 2
  IF (.NOT. DOUBLE) GO TO 730
  k = 4
  730 WRITE (nout,450) (blank4,j=1,word),(equal4,j=1,k)
  nogo = 1
  GO TO 800
END DO
GO TO 800

760 IF (xsort == 0) kard1(ie) = khrk( 1)
770 IF (xsort == 0) kard1( 1) = khrk(39)
ibuf(1) = -1
CALL khrbcd (card80,bufx)
GO TO 800

780 ibuf(1) = 0

800 RETURN

810 IF (xsort == 2) GO TO 830
WRITE  (nout,820) xsort
820 FORMAT (//,' *** TABLE IN XREAD HAS NOT BEEN INITIALIZED.',  &
    /5X,'XSORT=',i4)
CALL mesage (-37,0,sub)
830 WRITE  (nout,840) card8
840 FORMAT (/,' *** ILLEGAL CHARACTER ENCOUNTERED IN INPUT CARD',  &
    /4X,1H',10A8,1H' )
nogo = 1
850 RETURN 1


ENTRY yread (*,bufx)
!     ====================

!     YREAD IS CALLED ONLY BY XSORT TO RE-PROCESS CARD IMAGES FROM
!     THE OPTP FILE

CALL bcdkh8 (bufx,card80)
CALL k2b (card80,card1,80)
fromy = 1
GO TO 80


ENTRY rmveq (bufx)
!     ==================

!     RMVEQ, CALLED ONLY BY XCSA, REMOVES AN EQUAL SIGN FROM TEXT.
!     THUS, 1 EQUAL SIGN BEFORE COLUMN 36 IS ALLOWED ON ONE EXECUTIVE
!     CONTROL LINE

!     AT THIS POINT, THE DATA IN KARD1 IS STILL GOOD

DO  i = 1,36
  IF (kard1(i) == equ1) GO TO 910
END DO
GO TO 920
910 kard1(i) = blankk
CALL khrbcd (card80,bufx)
920 RETURN
END SUBROUTINE xread
