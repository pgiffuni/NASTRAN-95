SUBROUTINE matgen
     
!     THE PURPOSE OF THIS MODULE IS TO GENERATE CERTAIN KINDS OF
!     MATRICES ACCORDING TO ONE OF SEVERAL SIMPLE USER SELECTED OPTIONS
 
!     MATGEN      TAB/OUT/P1/P2/P3/P4/P5/P6/P7/P8/P9/P10/P11  $
 
!     TAB - INPUT TABLE - (OPTIONAL) FOR USE IN GENERATING THE MATIRX
!                 (THIS DATA MAY BE ASSUMED TO BE INPUT VIA DTI CARDS.)
!             = EQEXIN TABLE  FOR P1 =  9
!             = USET   TABLE  FOR P1 = 11
!             = ANY GINO FILE FOR P1 = 10
 
!     OUT - OUTPUT MATRIX - IF PURGED AND P1 IS NOT 10, P1 WILL BE SET
!               TO -1 AND RETURN
 
!     P1      - INPUT - INTEGER, OPTION SELECTION. (DEFULAT P1 = 3)
!             = 1, GENERATE A RSP IDENTITY MATRIX OF ORDER P2.
!             = 2, GENERATE AN IDENTITY MATRIX OF ORDER P2, FORM 8
!             = 3, GENERATE A DIAGONAL MATRIX FORM INPUT FILE T
!             = 4, GENERATE A PARTERN MATRIX
!             = 5, GENERATE A MATRIX OF PSEUDO-RANDOM NUMBERS.
!             = 6, GENERATE PARTITION VECTOR OF ORDER P2, WITH P3 ZERO'S
!                  CLOOWED BY P4 ONE'S FOLLOWED BY P5 ZERO'S ETC.
!                  REMAINER IS ALWAYS AERO. TOO MANY DEFINITIONS IS AN
!                  ERROR.
!             = 7, GENERATE A NULL MATRIX
!             = 8, GENERATE A MATRIX FROM EQUATIONS BASED ON ITS INDICES
!             = 9, GENERATE A TRANSFORMATION BETWEEN EXTERNAL AND
!                  INTERANL MATRICES, OF G-SET SIZE.
!                  P2 = 0, OUTPUTS INT-EXT (DEFAULT)
!                  P2 = 1, OUTPUTS TRANSPOSE EXT-INT
!                  P3 = NO. OF TERMS IN G-SET (REQUIRED). USE LUSET IN
!                       MOST SOLUTION SEQUENCES.
!             =10, ALLOW USER TO ALTER DATA BLOCK TRAILER.
!             =11, GENERATE A RECTANGULAR MATRIX, DRIVEN BY USET TABLE
 
!     P2 - P11 -   OPTION PARAMETERS - INTEGER - INPUT AND OUTPUT
!                  INPUT  AS OPTION VALUE (1 THRU NP)
!                  OUTPUT AS -1 IF AND ONLY IF OUTPUT DATA BLOCK IS PRE-
!                  PURGED
!                  DEFAULT VALUES FOR P2 THRU P11 ARE ZEROS
 
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL         lshift
 INTEGER :: mcb(7),nam(2),p(11),code(2),ix(12)
 REAL :: val,rx(7),tmp(7)
 DOUBLE PRECISION :: d(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
 COMMON /system/  sysbuf,nout,dum37(37),nbpw,dum14(14),iprec
 COMMON /machin/  machx
 COMMON /zblpkx/  val(4),row
 COMMON /packx /  ita,itb,i2,j2,incr2
!ZZ   COMMON /ZZMGEN/  X(1)
 COMMON /zzzzzz/  x(20000)
 EQUIVALENCE      (val(1),d(1)),(x(1),ix(1),rx(1)),(p(1),p1)
 DATA    nam   /  4HMATG,4HEN  /,  np /  11   /
 DATA    eqe   ,  xin / 4HEQEX, 4HIN  /, code / 6,1 /
 DATA    out   ,  t      / 201   ,  101    /
 
!     IF OUTPUT DATA BLOCK IS PRE-PURGED, SET P1 = -1 AND RETURN
 
 IF (p1 == 10) GO TO 30
 mcb(1) = out
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 1580
 
!     CHECK INPUT FILE REQUIREMENT
 
 30 ix(1) = t
 CALL rdtrl (ix(1))
 IF (p1 == 10) GO TO 1000
 IF (p1 == 3 .OR. p1 == 9 .OR. p1 == 10 .OR. p1 == 11) GO TO 50
 IF (ix(1) == 0) GO TO 50
 CALL fname (ix(1),rx(2))
 WRITE  (nout,40) uwm,rx(2),rx(3),p1
 40 FORMAT (a25,' FROM MODULE MATGEN. INPUT DATA BLOCK ',2A4,' IS ',  &
     'NOT NEEDED FOR OPTION',i3)
 
!     CHECK OPEN CORE AND OPEN OUTPUT DATA BLOCK
 
 50 lcor = korsz(ix(1))
 IF (p1 == 2) GO TO 200
 IF (lcor < sysbuf) GO TO 1500
 ibuf1 = lcor  - sysbuf - 1
 ibuf2 = ibuf1 - sysbuf
 CALL gopen (out,ix(ibuf1),1)
 lcor = lcor - sysbuf
 
!     TEST FOR VALID OPTION AND BRANCH ON OPTION.
 
 IF (p1 == 0) p1 = 3
 IF (p1 < 0 .OR. p1 > np) GO TO 1510
 SELECT CASE ( p1 )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 300
   CASE (    4)
     GO TO 400
   CASE (    5)
     GO TO 500
   CASE (    6)
     GO TO 600
   CASE (    7)
     GO TO 700
   CASE (    8)
     GO TO 800
   CASE (    9)
     GO TO 900
   CASE (   10)
     GO TO 1000
   CASE (   11)
     GO TO 1100
 END SELECT
 
!     OPTION 1 - GENERATE A RSP IDENTITY MATRIX OF ORDER P2, AND TRAILER
!     ========   P2 = ORDER OF MATRIX
!                P3 = SKEW FLAG, IF NONZERO, GENERATE A SKEW-DIAGONAL
!                     MATRIX
!                P4 = PRECISION (1 OR 2). IF ZERO, USE MACHINE PRECISION
 
 100 IF (p2 > 0) GO TO 110
 ipx = 2
 px  = p2
 GO TO 1530
 110 ita = 1
 itb = p4
 IF (p4 == 0) itb = iprec
 incr2 = 1
 CALL makmcb (mcb,out,p2,6,itb)
 DO  i = 1,p2
   rx(i) = 1.0
   i2 = i
   j2 = i
   CALL pack (ix,out,mcb)
 END DO
!WKBI SPR93023 12/93
 CALL CLOSE ( out, 1 )
 GO TO 210
 
!     OPTION 2 - GENERATE AN IDENTITY TRAILER (FORM = 8)
!     ========   P2 = ORDER OF MATRIX
 
!                ** CAUTION ** FORM = 8 MATRICES DO NOT REALLY EXIST
!                ONLY CERTAIN  ROUTINES CAN PROCESS THEM
!                e.g. FBS, MPYAD, CEAD etc.
 
 200 mcb(1) = out
 mcb(2) = p2
 mcb(3) = p2
 mcb(4) = 8
 mcb(5) = 1
 mcb(6) = 1
!     MCB(7) = LSHIFT(1,NBPW-2) + P2
 mcb(7) = lshift(1,nbpw-2 - (nbpw-32)) + p2
 
!     ADD (NBPW-32) TO MCB(7) SO THAT CRAY, WITH 48-BIT INTEGER WILL
!     NOT GET INTO TROUBLE. (SEE SDCOMP AND WRTTRL)
 
 210 CALL wrttrl (mcb)
 GO TO 1700
 
!     OPTION 3 - GENERATE A DIAGONAL MATRIX FROM INPUT TABLE T
!     ========   P2 = DATA TYPE OF T
!                P3 = 0, FORM 6 MATRIX IS GENERATED
!                   = 1, FORM 3 MATRIX IS GENERATED
 
!     THIS OPTION IS THE ORIGINAL MATGEN IN COSMIC MATGEN
!     SKIP HEADER RECORD, AND BEGINNING RECORD ON T
!     PICKUP DATA IN ARRAY OF 7 WORDS. DIAGONAL VAULE ON THE 3RD
 
 300 lcor = lcor - sysbuf
 IF (lcor < sysbuf) GO TO 1500
 IF (p2 == 0) p2 = 1
 CALL OPEN (*1550,t,ix(ibuf2),0)
 CALL skprec (t,2)
 ita = p2
 itb = iprec
 incr2 = 1
 FORM = 6
 IF (p3 == 1) FORM = 3
 CALL makmcb (mcb,out,0,FORM,iprec)
 m  = 0
 310 CALL READ (*1560,*330,t,tmp,7,0,0)
 m  = m + 1
 IF (p3 == 1) GO TO 320
 i2 = m
 j2 = m
 CALL pack (tmp(3),out,mcb)
 GO TO 310
 320 rx(m) = tmp(3)
 GO TO 310
 330 IF (p3 == 1) GO TO 340
 mcb(3) = mcb(2)
 GO TO 350
 340 i2 = 1
 j2 = m
 CALL pack (rx,out,mcb)
 mcb(2) = 1
 mcb(3) = m
 350 CALL CLOSE (t,1)
 GO TO 210
 
!     OPTION 4 - GENERATE A PATTERN MATRIX
!     ========   P2 = NUNBER OF COLUMNS
!                P3 = NUMBER OF ROWS
!                P4 = PRECISION (1 OR 2). IF 0, USE MACHINE PRECISION
!                P5 = NUMBER OF TERMS PER STRING. IF 0, USE 1
!                P6 = INCREMENT BETWEEN STRINGS. IF 0, USE 1
!                P7 = ROW NUMBER OF 1ST STRING IN COLUMN 1. IF 0, USE 1
!                P8 = INCREMENT TO 1ST ROW FOR SUBSEQUENT COLUMNS.
!                P9 = NUMBER OF COLS BEFORE RETURNING TO P7.
 
!                THE VALUE OF EACH NON-ZERO TERM IN THE MATRIX WILL BE
!                THE COLUMN NUMBER
!                e.g. TO GENERATE A 10x10 DIAGONAL MATRIX WITH THE COL.
!                NUMBER IN EACH DIAGONAL POSITION, CODE
 
!                MATGEN   ,/DIAG/4/10/10/0/1/10/1/1/10  $
 
 400 p2 = MAX0(p2,1)
 p3 = MAX0(p3,1)
 IF (p4 /= 1 .AND. p4 /= 2) p4 = 0
 IF (p4 == 0) p4 = iprec
 p5 = MAX0(p5,1)
 p6 = MAX0(p6,1)
 p7 = MAX0(p7,1)
 p8 = MAX0(p8,0)
 p9 = MAX0(p9,1)
 irow1 = p7
 l = 1
 CALL makmcb (mcb,out,p3,2,p4)
 
 DO  j = 1,p2
   IF (p4 == 1) val(1) = j
   IF (p4 == 2) d(  1) = j
   row = irow1
   CALL bldpk (p4,p4,out,0,0)
   410 CONTINUE
   DO  k = 1,p5
     IF (row > p3) GO TO 430
     CALL zblpki
     row = row + 1
   END DO
   row = row + p6 - 1
   GO TO 410
   430 CALL bldpkn (out,0,mcb)
!WKBI 9/93
   435 CONTINUE
   l = l + 1
   irow1 = irow1 + p8
!WKBR 9/93 IF (L .LE. P9) GO TO 430
   IF (l <= p9) GO TO 435
   l = 1
   irow1 = p7
 END DO
 GO TO 1400
 
!     OPTION 5 - GENERATE A MATRIX OF PSEUDO-RANDOM NUMBERS. THE NUMBERS
!     ========   SPAN THE RANGE 1. TO 1.0 WITH A NORMAL DISTRIBUTION
!                P2 = NUMBER OF COLUMNS
!                P3 = NUMBER OF ROWS
!                P4 = PRECISION (1 OR 2).  IF 0, USED MACHINE PRECISION
!                P5 = SEED FOR RANDOM NUMBER GENERATION.  IF P5.LE.0,
!                     THE TIME OF DAY (SECONDS PAST MIDNIGHT) WILL BE
!                     USED
 
!     OPTION 5 WAS WRITTEN BY G.CHAN/UNISYS 2/93
 
 500 ita = 1
 itb = p4
 IF (p4 == 0) itb = iprec
 FORM = 2
 IF (p2 == p3) FORM = 1
 i2 = 1
 j2 = p3
 incr2 = 1
 CALL makmcb (mcb,out,p2,FORM,itb)
 k  = p5
 IF (machx == 4) GO TO 560
!                   CDC
 IF (machx == 9) GO TO 530
!                    HP
 
 DO  i = 1,p2
   IF (p5 == 0) CALL cputim (k,k,0)
   k  = (k/2)*2 + 1
   DO  j = 1,p3
!WKBR 5/95 SUN     RX(J) = RAN(K)
     rx(j) = rand(k)
   END DO
   CALL pack (rx(1),out,mcb)
 END DO
 GO TO 590
 
!     HP ONLY
!     ACTIVATE SRAND AND RAND() BELOW, AND COMMENT OUT RAN(K) ABOVE
 
 530 CONTINUE
 WRITE  (nout,535) sfm
 535 FORMAT (a25,'. MATGEN NEEDS TO ACTIVATE SRAND AND RAND() FOR HP')
 CALL mesage (-61,0,0)
 DO  i = 1,p2
   IF (p5 == 0) CALL cputim (k,k,0)
!     CALL SRAND (K)
   DO  j = 1,p3
!     RX(J) = RAND()
   END DO
   CALL pack (rx(1),out,mcb)
 END DO
 GO TO 590
 
!     CDC ONLY
!     ACTIVATE SRAND AND RAND() BELOW, AND COMMENT OUT RAN(K) ABOVE
 
 560 CONTINUE
 WRITE  (nout,565) sfm
 565 FORMAT (a25,'. MATGEN NEEDS TO ACTIVATE RANSET AND RANF() FOR CDC' )
 CALL mesage (-61,0,0)
 DO  i = 1,p2
   IF (p5 == 0) CALL cputim (k,k,0)
!     CALL RANSET (K)
   DO  j = 1,p3
!     RX(J) = RANF()
   END DO
   CALL pack (rx(1),out,mcb)
 END DO
 
 590 CALL CLOSE  (out,1)
 CALL wrttrl (mcb(1))
 GO TO 1700
 
!     OPTION 6 - GENERATE A PARTITIONING VECTOR FOR USE IN PARTN OR
!     ========   MERGE
!                P2 = NUMBER OF ROWS
!                P3,P5,P7,P9  = NUMBER OF ROWS WITH ZERO COEFFICIENTS
!                P4,P6,P8,P10 = NUMBER OF ROWS WITH UNIT COEFFICIENTS
 
!                IF SUM OF P3 THRU P10 IS .LT. P2, THE REMAINING TERMS
!                CONTAIN ZEROS
!                IF SUM OF P3 THRU P10 IS .GT. P2, THE TERMS ARE IGNORED
!                AFTER P2
!                e.g. GENERATE A VECTOR OF 5 UNIT TERMS FOLLOWED BY 7
!                ZEROS, FOLLOWED BY TWO UNIT TERMS
 
!                MATGEN,   ,/UPART/6/14/0/5/7/2   $
 
!     OPTION 6 WAS ORIGINALLY WRITTEN BY P.KIRCHMAN/SWALES 1/92
!     RE-CODED BY G.CHAN/UNISYS FOR ALL COMPILERS,  2/93
 
 600 ipx = 2
 px  = p2
 IF (p2 <= 0) GO TO 1530
 incr2 = 1
 i2  = 1
 j2  = p2
 ita = 1
 itb = 1
 CALL makmcb (mcb,out,p2,2,itb)
 tot = 0
 DO  i = 3,11
   tot = tot + p(i)
 END DO
 IF (tot > p2) WRITE (nout,620) ufm,p1,p2
 IF (tot < p2) WRITE (nout,630) uwm,p1
 620 FORMAT (a23,' FROM MATGEN, OPTION',i3,'. TOO MANY ENTRIES FOR ',  &
     'SPECIFIED SIZE',i7)
 630 FORMAT (a25,' FORM MATGEN, OPTION',i3,'. THE NUMBER OF ENTRIES ',  &
     'SPECIFIED BY PARAMETERS IS LESS THAN THE TOTAL SIZE', /5X,  &
     'OF THE PARTITION. THE REMAINING RENTRIES ARE ZERO FILLED')
 k = 1
 DO  i = 3,9,2
   pi = p(i)
   DO  j = 1,pi
!WKBR SPR 93023 12/93      IX(K) = 0
     rx(k) = 0.
     k = k + 1
   END DO
   pi = p(i+1)
   DO  j = 1,pi
!WKBR SPR 93024 12/93      IX(K) = 1
     rx(k) = 1.0
     k = k + 1
   END DO
 END DO
 IF (k >= p2) GO TO 680
 DO  i = k,p2
   ix(i) = 0
 END DO
 680 CALL pack (ix,out,mcb)
 CALL CLOSE (out,1)
 CALL wrttrl (mcb)
 GO TO 1700
 
!     OPTION 7 - GENERATE A NULL MATRIX
!     ========   P2 = NUMBER OF ROWS
!                P3 = NUMBER OF COLUMNS
!                P4 = FORM; IF P4 = 0, AND P2 = P3, FORM WILL BE 6
!                     (SYMMETRIC). OTHERWISE P4 IS 2 (RECTANGULAR)
!                P5 = TYPE: IF P5 = 0, TYPE IS MACHINE PRECISION
 
 700 d(1) = 0.0D0
 d(2) = 0.0D0
 ita = 1
 itb = p5
 IF (p5 == 0) itb = iprec
 FORM = p4
 IF (p4 == 0 .AND. p2 == p3) FORM = 6
 IF (p4 == 0 .AND. p2 /= p3) FORM = 2
 i2 = 1
 j2 = 1
 incr2 = 1
 CALL makmcb (mcb,out,p2,FORM,itb)
 DO  i = 1,p3
   CALL pack (val,out,mcb)
 END DO
 CALL CLOSE  (out,1)
 CALL wrttrl (mcb(1))
 GO TO 1700
 
!     OPTION 8 - GENERATE A MATRIX FROM EQUATIONS BASED ON IT INDICES
!     ========   P2 =  0, GENERATE ALL TERMS
!                  .NE.0, GENERATE ONLY DIAGONAL TERMS
!                P3 =     NUMBER OF ROWS
!                P4 =     NUMBER OF COLUMNS
!                P5 =     NUMBER OF THE RECORD IN THE INPUT DTI TABLE
!                         USED TO DEFINE REAL COEFFICIENTS
!                  .LT.0, COEFFICIENT TAKEN FROM DTI TRAILER
!                         COEFF(TRAILER1) = FLOAT(TRAILER2)   TRAILER
!                         COEFF(TRAILER3) = FLOAT(TRAILER4)   ITEMS ARE
!                         COEFF(TRAILER5) = FLOAT(TRAILER6)   INTEGERS
 
!                   = 0,  DATA PAIRS FROM RECORD 0 (DATA BLOCK HEADER
!                         RECORD) ARE INTERPRETED AS DFINING
!                         COEFF(V1) = V2     V1 IS INTEGER, V2 IS REAL
!                  .GT.0, DATA PAIRS FROM RECORD P5 INTERPRETED AS ABOVE
!                P6 =     NUMBER OF THE RECORD IN THE INPUT DTI TABLE
!                         USED TO DEFINE IMAGINARY DOEFFICIENTS D(I)
!                  .LE.0, NO DOEFFICIENTS DEFINED
!                  .GT.0, DATA PAIRS FROM RECORD P6 INTERPRETED AS ABOVE
!                         WHERE D(V1) = V2
!                P7 =     FORM OF OUTPUT MATRIX
!                  .LE.0, FORM = 1 OR 2, DEPENDING ON P3 AND P4
!                  .GT.0, FORM SET TO P7
!                P8 =     COEFFICIENT PRINT FLAG
!                   = 0,  DO NOT PRINT COEFFICIENT LISTS
!                  .NE.0, PRINT COEFFICIENTS LISTS C(L) AND D(L) FROM
!                         DTI INPUT. (PRINT D(L) LIST ONLY IF P6.GT.0)
 
!                SEE USER MANUAL FOR THE EQUATION USED TO DETERMINE THE
!                COEFFICIENT OF THE (I,J)TH TERM OF THE OUTPUT MATRIX
 
 800 WRITE (nout,1200) uwm,p1
 GO TO 1700
 
!     OPTION 9 - GENERATE A TRANSFORMATION BETWEEN EXTERNAL AND INTERNAL
!     ========   MATRICES FOR G-SET SIZE MATRICES
!                P2 = 0, OUTPUT NON-TRANSPOSED FACTOR, UEXT = MAT*UINT
!                   = 1, OUTPUT TRANSPOSED FACTOR, UEXT = MAT*UINT
!                P3 = NUMBER OF TERMS IN G-SET. THE PARAMETER LUSET
!                     CONTAINS THIS NUMBER IN MOST SOLUTION SEQUENCES
 
!                EXAMPLES -
!                1. TRANSFORM A g-SET SIZE VECTOR TO EXTERNAL SEQUENCE
!                ALTER XX  $ AFTER SDR1, ALL SDR1 OUTPUTS ARE IN
!                            INTERNAL SEQUENCE
!                MATGEN  EQEXIN/INTEXT/9//LUSET $
!                MPYAD   INTEXT,UGV,/UGVEXT/1 $
 
!                2. TRANSFORM AN a-SET SIZE MATRIX TO EXTERNAL SEQUENCE
!                ALTER XX  $ AFTER KAA IS GENERATED, ALL MATRICES ARE IN
!                            INTERNAL SEQUENCE
!                MATGET  EQEXIN/INTEXT/9/0/LUSET $
!                SMPYAD  INTEXT,KAGG,INTEXT,,/KAAGEXT/3////1////6 $
!                $ (KAAGEXT) = TRANSPOSE(INTEXT)*(KAAG)*(INTEXT)
!                $ ITS FORM IS 6 (SYMMETRIC)
 
!     OPTION 9 WAS ORIGINALLY WRITTEN BY P.KIRCHMAN/SWALES 1/92
!     RE-CODED BY G.CHAN/UNISYS FOR ALL COMPILERS,  2/93
 
 900 ipx = 3
 px  = p(3)
 IF (px <= 0) GO TO 1530
 nuset = px
 l     = 2
 nval  = ix(l)
 CALL fname (t,tmp(1))
 IF (tmp(1) /= eqe .OR. tmp(2) /= xin) GO TO 1600
 CALL OPEN (*1550,t,ix(ibuf2),0)
 CALL fwdrec (*1560,t)
 CALL fwdrec (*1560,t)
 CALL READ (*1560,*910,t,ix(1),ibuf2-1,1,l)
 910 CALL CLOSE (t,1)
 IF (l /= nval*2) GO TO 1620
 ita = iprec
 itb = ita
 CALL makmcb (mcb,out,nuset,2,itb)
 incr2 = 1
 val(1)= 1.0
 IF (ita == 2) d(1) = 1.0D+0
 tot   = 0
 IF (p2 > 0) GO TO 930
 
!     NO TRANSPOSE
 
 DO  i = 1,nval
   is2 = i*2
   a = ix(is2)/10
   b = MOD(ix(is2),10)
   c = code(b)
   DO  j = 1,c
     i2 = a
     j2 = a
     CALL pack (val,out,mcb)
     a = a + 1
   END DO
 END DO
 tot = tot + c
 GO TO 980
 
!     TRANSPOSE
 
 930 nval2 = nval*2
 pos = 1
 DO  i = 1,nval
   is2 = i*2
   a = ix(is2)/10
   b = MOD(ix(is2),10)
   ix(is2-1) = pos
   pos = pos + code(b)
 END DO
 DO  i = 4,nval2,2
   j = nval2
   flag = 0
   950 IF (ix(j) >= ix(j-2)) GO TO 960
   flag = 1
   k = ix(j  )
   l = ix(j-1)
   ix(j  ) = ix(j-2)
   ix(j-1) = ix(j-3)
   ix(j-2) = k
   ix(j-3) = l
   960 j = j - 2
   IF (j >= i) GO TO 950
   IF (flag == 0) GO TO 980
 END DO
 
 980 DO  i = 1,nval
   is2 = i*2
   a = ix(is2)/10
   b = MOD(ix(is2),10)
   a = ix(is2-1)
   c = code(b)
   DO  j = 1,c
     i2 = a
     j2 = a
     CALL pack (val,out,mcb)
     a  = a + 1
   END DO
 END DO
 tot = nval*c
 IF (nuset /= tot) GO TO 1640
 CALL wrttrl (mcb)
 CALL CLOSE  (out,1)
 GO TO 1700
 
!     OPTION 10 - ALLOW USER TO ALTER DATA BLOCK TRAILER
!     =========
!     IF PI IS NEGATIVE, THE CORRESPONDING TRAILER WORD (I) IS SET TO
!     ZERO
 
 1000 IF (ix(1) == 0) GO TO 1050
 CALL fname (ix(1),ix(11))
 WRITE  (nout,1010) uim,ix(11),ix(12),(ix(i),i=2,7)
 1010 FORMAT (a29,' FROM MATGEN MODULE, OPTION 10. TRAILER OF ',2A4,2H -  &
     ,      /5X,'OLD - ',6I7)
 DO  i = 2,7
   IF (p(i) /= 0) ix(i) = p(i)
   IF (p(i) < 0) ix(i) = 0
 END DO
 WRITE  (nout,1030) (ix(i),i=2,7)
 1030 FORMAT (5X,'NEW - ',6I7)
 IF (ix(2) == ix(3) .AND. ix(4) == 2 .AND. ix(7) /= 0) WRITE (nout,1040) uim
 1040 FORMAT (a29,'. SINCE ROW = COLUMN, RECTANGULAR FORM 2 WILL BE ',  &
     'CHANGED TO SQUARE FORM 1 AUTOMATICALLY')
 ix(1) = 199
 CALL wrttrl (ix(1))
 GO TO 1700
 
 1050 WRITE  (nout,1060) uwm
 1060 FORMAT (a25,' FROM MATGEN, OPTION 10. INPUT FILE MISSING')
 GO TO 1700
 
!     OPTION 11 - GENERATE A RECTANGULAR MATRIX, DRIVEN BY USET TABLE
!     =========   P2 = 1, GENERATE A NULL MATRIX
!                  .NE.1, GENERATE A NULL MATRIX WITH AN IDENTITY MATRIX
!                         STORED IN IT
!                 P3 = NUMBER OF COLUMNS OF OUTPUT MATRIX, IF P2 = 1
!                    = BIT POSITION OF SET THAT DEFINES NUMBER OF ROW,
!                      IF P2.NE.1. SEE SECTION 1.4.10 FOR BIT POSITION
!                      LIST. DEFAULT IS A-SET SIZE.
!                 P4 = NOT USED IF P2 = 1. THE OUTPUT MATRIX WILL BE
!                      NULL AND HAVE P3 COLUMNS AND A-SET SIZE ROWS
!                    = BIT POSITION OF SET THAT DEFINES NUMB OF COLUMNS
!                      IF P2.NE.1. DEFAULT IS L-SET SIZE
 
!                 IF P2.NE.1, AND ONE OR BOTH OF THE SETS REQUESTED IN
!                 P3 AND P4 DOES NOT EXIST, THEN MAT IS RETURNED PURGED,
!                 AND P5 IS RETURNED WITH THE VALUE OF -1. IF MAT DOES
!                 EXISTS, P5 IS RETURNED WITH THE VALUE 0
 
 1100 WRITE  (nout,1200) uwm,p1
 1200 FORMAT (a25,' FROM MATGEN MODULE, OPTION',i3,' IS NOT AVAILABLE')
 GO TO 1700
 
!     WRAP-UP AND RETURN TO EXECUTIVE SYSTEM
 
 1400 CALL CLOSE  (out,1)
 CALL wrttrl (mcb)
 GO TO 1700
 
!     ERROR MESSAGES
 
 1500 CONTINUE
 lcor = sysbuf - lcor
 CALL mesage (-8,lcor,nam)
 GO TO 1700
 
 1510 WRITE  (nout,1520) ufm,p1
 1520 FORMAT (a23,' IN MATGEN, ILLEGAL VALUE FOR OPTION PARAMETER =',i5)
 GO TO 1690
 
 1530 WRITE  (nout,1540) ufm,ipx,px
 1540 FORMAT (a23,' IN MATGEN, ILLEGAL VALUE FOR PARAMETER ',i1,3H = , i5)
 
 1550 j = -1
 GO TO 1570
 1560 j = -2
 1570 CALL mesage (j,t,nam)
 
 1580 WRITE  (nout,1590) ufm,p1
 1590 FORMAT (a23,'. OPTION',i3,' OUTPUT DATA BLOCK IS MISSING')
 p1 = -1
 GO TO  1700
 1600 WRITE  (nout,1610) ufm,tmp(1),tmp(2)
 1610 FORMAT (a23,'. OPTION 9. INPUT FILE IS ',2A4,', NOT EQEXIN')
 GO TO  1690
 1620 WRITE  (nout,1630) ufm,l,nval
 1630 FORMAT (a23,'. EQEXIN RECORD LENGTH NOT MATCH TWICE TRAIL(2)',2I9)
 GO TO  1690
 1640 WRITE  (nout,1650) ufm,nuset,tot
 1650 FORMAT (a23,'. OPTION 9, LUSET OF',i9,' DOES NOT AGREE WITH SIZE',  &
     ' OF EQEXIN',i9)
 1690 CALL mesage (-61,0,nam)
 1700 RETURN
END SUBROUTINE matgen
