SUBROUTINE inptt2
     
!     READ DATA BLOCK(S) FROM A FORTRAN UNIT.
 
!     CALL TO THIS MODULE IS
 
!     INPUTT2   /O1,O2,O3,O4,O5/V,N,P1/V,N,P2/V,N,P3/V,N,P4/V,N,P5/
!                               V,N,P6 $
 
!     PARAMETERS P1, P2, P4, AND P5 ARE INTEGER INPUT. P3 AND P6 ARE BCD
 
!            P1 =+N, SKIP FORWARD N DATA BLOCKS BEFORE READ
!               = 0, NO ACTION TAKEN BEFORE READ (DEFAULT)
!               =-1, BEFORE READ, FORTRAN TAPE IS REWOUND AND TAPE
!                    HEADER RECORD (RECORD NUMBER ZERO) IS CHECKED
!               =-3, THE NAMES OF ALL DATA BLOCKS ON FORTRAN TAPE
!                    ARE PRINTED AND READ OCCURS AT BEGINNING OF TAPE
!               =-5, SEARCH FORTRAN TAPE FOR FIRST VERSION OF DATA
!                    BLOCKS REQUESTED.
!                    IF ANY ARE NOT FOUND, A FATAL TERMINATION OCCURS.
!               =-6, SEARCH FORTRAN TAPE FOR FINAL VERSION OF DATA
!                    BLOCKS REQUESTED.
!                    IF ANY ARE NOT FOUND, A FATAL TERMINATION OCCURS.
!               =-7, SEARCH FORTRAN TAPE FOR FIRST VERSION OF DATA
!                    BLOCKS REQUESTED.
!                    IF ANY ARE NOT FOUND, A WARNING OCCURS.
!               =-8, SEARCH FORTRAN TAPE FOR FINAL VERSION OF DATA
!                    BLOCKS REQUESTED.
!                    IF ANY ARE NOT FOUND, A WARNING OCCURS.
 
!            P2 =    THE FORTRAN UNIT FROM WHICH THE DATA BLOCK(S)
!                    WILL BE READ. (DEFAULT P2 = 11, OR 14)
 
!            P3 =    TAPE ID CODE FOR FORTRAN TAPE, AN ALPHANUMERIC
!                    VARIABLE WHOSE VALUE MUST MATCH A CORRESPONDING
!                    VALUE ON THE FORTRAN TAPE.
!                    THIS CHECK IS DEPENDENT ON THE VALUE OF P1 AS
!                    FOLLOWS..
 
!                    *P1*             *TAPE ID CHECKED*
!                     +N                     NO
!                      0                     NO
!                     -1                    YES
!                     -3                    YES (WARNING CHECK)
!                     -5                    YES
!                     -6                    YES
!                     -7                    YES
!                     -8                    YES
!                    THE MPL DEFAULT VALUE FOR P3 IS XXXXXXXX .
 
!            P4 =    NOT USED IN INPUTT2.
!                    (USED ONLY IN OUTPUT2 FOR MAXIMUM RECORD SIZE)
 
!            P5 = 0, NON-SPARSE MATRIX IF INPUT IS A MATRIX DATA BLOCK
!               = NON-0, SPARSE MATRIX IF INPUT IS A MATRIX DATA BLOCK
!                    (P4 IS IGNORED IF INPUT IS A TALBE DATA BLOCK.
!                     P4 IS EQUIVALENT TO P5 IN OUTPUT2 MODULE)
 
!            P6 = BLANK, (DEFAULT)
!               = 'MSC', THE INPUT TAPE WAS WRITTEN IN MSC/OUTPUT2
!                     COMPATIBEL RECORD FORMAT.
 
!     OUTPT2 DOES NOT AUTOMATICALLY OUTPUT THE MATRIX IN STRING OR
!     SPARSE FORM. UNLESS P5 IS REQUESTED.
!     SIMILARILY, INPUT2 DOES NOT AUTOMATICALLY PROCESS MATRIX IN SPARSE
!     MATRIX FORM, UNLESS P5 IS REQUESTED).
 
!     REVISED  11/90 BY G.CHAN/UNISYS
!              (1) TO ACCEPT MSC/OUTPUT2 DATA (CALLED FROM INPTT4, 11/90
!                  OR INPTT2, 2/93)
!              (2) TO ACCEPT SPARSE MATRIX COMING FORM COSMIC/OUTPT2
!                  (SEE P5 PARAMETER IN INPTT2 AND P5 IN OUTPT2)
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: sparse,dp
 INTEGER :: trl(8),NAME(2),subnam(2),typin,mcb(7),dx(3),  &
     namex(2),idhdr(7),idhdrx(7),p3x(2),nt(5,3),  &
     tapcod(2),NONE(2),bcdbin(4),blk(20)
 REAL :: core(1)
 DOUBLE PRECISION :: dcore(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON /dsname/ dsnames(80)
!WKBNE
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3(2),p4,p5,p6(2) /system/  ksystm(65)  &
     /packx /  typin,typout,irow,nrow,incr /TYPE  /  prec(2),nwds(4)  &
     /zzzzzz/  x(1)
 EQUIVALENCE      (core(1),x(1))
 EQUIVALENCE      (ksystm(1),nb    ), (ksystm( 2),nout ),  &
     (ksystm(9),nlpp  ), (ksystm(12),line ),  &
     (blk( 1)  ,bname ), (blk( 2)   ,btyp ),  &
     (blk( 3)  ,bform ), (blk( 4)   ,brow ),  &
     (blk( 5)  ,bpoint), (blk( 6)   ,brav ),  &
     (blk( 7)  ,bwrt  ), (blk( 8)   ,bflag),  &
     (blk(12)  ,bcol  ), (dcore(1),core(1))
!WKBI
 DATA    ifirst / 0 /
 DATA    subnam/  4HINPT, 4HT2  /  ,   NONE / 4H (NO,4HNE)       /
 DATA    zero  ,  mone,mtwo,mtre,mfor /0,-1,-2,-3,-4  /, i3 / 3  /  &
         mfiv  ,  msix,mete /-5,-6,-8 /,    ipt2,ipt4 / 1H2, 1H4 /
 DATA    idhdr /  4HNAST,4HRAN ,4HFORT,4H TAP,4HE ID,4H COD,3HE -/
 DATA    bcdbin/  4HBCD ,4H    ,4HBINA,4HRY    /, msc / 4HMSC    /
 
 
 iptx  = ipt2
 ntrl  = 8
 IF (p4 == 0) GO TO 20
 IF (p5 /= 0) GO TO 20
 WRITE  (nout,10) uwm
 10 FORMAT (a25,'. THE 4TH PARAMETER IN INPUTT2 MODULE IS NO LONGER ',  &
     'USED.', /5X,'SPARSE MATRIX FLAG IS NOW THE 5TH PARAMETER',  &
     ', A MOVE TO SYNCHRONIZE THE PARAMETERS USED IN OUTPUT2')
 p5 = p4
 20 sparse = .false.
 IF (p5 /= 0) sparse = .true.
 IF (p6(1) /= msc) GO TO 100
 GO TO 50
 
 
 ENTRY input2
!     ============
 
!     INPUT2 IS CALLED TO HANDLE MSC/OUTPUT2 DATA.
!     IT IS CALLED FROM INPTT2 WITH P6 PARAMETER = 'MSC', OR
!     FROM INPTT4
 
 iptx  = ipt4
 50 WRITE  (nout,60) uim,iptx
 60 FORMAT (a29,' FROM INPUTT',a1,'. USER INPUT TAPE IN MSC/OUTPUT2',  &
     ' COMPATIBLE RECORDS')
 iptx  = ipt4
 ntrl  = 7
 irecf = 0
 sparse= .false.
 IF (p3(1) == bcdbin(1) .AND. p3(2) == bcdbin(2)) GO TO 1580
 IF (p3(1) == bcdbin(3) .AND. p3(2) == bcdbin(4)) GO TO 1580
 
 100 lcor = korsz(x) - nb
 IF (lcor <= 0) CALL mesage (-8,lcor,subnam)
!WKBNB
 IF ( ifirst /= 0) GO TO 61
 CLOSE ( UNIT=p2 )
 OPEN ( UNIT=p2, FILE=dsnames(p2), FORM='UNFORMATTED', STATUS='UNKNOWN' )
 ifirst = 1
 61    CONTINUE
!WKBNE
 oubuf = lcor + 1
 tapcod(1) = p3(1)
 tapcod(2) = p3(2)
 in = p2
 IF (p1 < mete .OR. p1 == mtwo .OR. p1 == mfor) GO TO 1420
 
 IF (p1 < mfor) GO TO 700
 IF (p1 == mtre) GO TO 500
 IF (p1 <= zero) GO TO 130
 
 i = 1
 110 READ (in) key
 keyx = 2
 IF (key /= keyx) GO TO 1530
 READ (in) namex
 READ (in) key
 imhere = 115
 IF (key >= 0) GO TO 1560
 ASSIGN 120 TO ret
 nskip = 1
 GO TO 1300
 
 120 i = i + 1
 IF (i <= p1) GO TO 110
 GO TO 160
 
!     OPEN FORTRAN TAPE TO READ TAPE-LABEL WITHOUT REWIND.
 
 130 IF (p1 /= mone) GO TO 160
 REWIND in
 READ (in) key
 keyx = 3
 IF (key /= keyx) GO TO 1530
 READ (in) dx
 READ (in) key
 keyx = 7
 IF (key /= keyx) GO TO 1530
 READ (in) idhdrx
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 1460
 END DO
 READ (in) key
 keyx = 2
 IF (key /= keyx) GO TO 1530
 READ (in) p3x
 READ (in) key
 imhere = 145
 IF (key >= 0) GO TO 1560
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 1440
 ASSIGN 150 TO ret
 nskip = 1
 GO TO 1300
 150 CONTINUE
 
 160 DO  i = 1,5
   
   output = 200 + i
   trl(1) = output
   CALL rdtrl (trl)
   IF (trl(1) <= 0) CYCLE
   CALL fname (output,NAME)
   IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) CYCLE
   
!     READ FILE NAME HEADER RECORD.
   
   READ (in) key
   IF (key == 0) GO TO 440
   keyx = 2
   IF (key /= keyx) GO TO 1530
   READ (in) namex
   READ (in) key
   imhere = 163
   IF (key >= 0) GO TO 1560
   
!     READ TRAILER RECORD.
   
   READ (in) key
   keyx = ntrl
   IF (key /= keyx) GO TO 1530
   READ (in) (trl(l),l=1,ntrl)
   IF (iptx == ipt2) irecf = trl(8)
   READ (in) key
   imhere = 165
   IF (key >= 0) GO TO 1560
   
!     OPEN OUTPUT DATA BLOCK TO WRITE WITH REWIND.
   
   CALL OPEN (*1400,output,x(oubuf),1)
   
!     COPY CONTENTS OF FORTRAN TAPE ONTO OUTPUT DATA BLOCK.
   
!     TRL(8) = 0, DATA BLOCK IS A TALBE
!            = 1, DATA BLOCK IS A MATRIX, WRITTEN IN STRING FORMAT
!            = 2, DATA BLOCK IS A VECTOR (1ST RECORD IS REGULAR, 2ND
!                 RECORD IS A STRING)
   
   INDEX = 0
   READ (in) key
   IF (iptx == ipt2) GO TO 180
   bname = output
   keyx  = 1
   IF (key /= keyx) GO TO 1530
   READ (in) krec
   imhere = 170
   IF (krec /= 0) GO TO 1560
   
   READ (in) key
   180 keyx  = 2
   IF (key < keyx) GO TO 1530
   IF (key > lcor) GO TO 1510
   READ (in) (x(l),l=1,key)
   CALL WRITE (output,NAME,2,0)
   IF (key == keyx) GO TO 200
   CALL WRITE (output,x(i3),key-2,0)
   
   200 IF (iptx == ipt2) GO TO 220
   READ (in) key
   imhere = 205
   IF (key >= 0) GO TO 1560
   btyp  = trl(5)
   bform = 0
   bcol  = 0
   nwd   = nwds(btyp)
   dp    = btyp == 2 .OR. btyp == 4
   CALL WRITE (output,x,0,1)
   210 READ (in) key
   keyx = 1
   IF (key /= keyx) GO TO 1530
   READ (in) krec
   IF (krec /= 0) GO TO 350
   
!     TABLE DATA BLOCK(S)
   
   220 READ (in) key
   IF (key < 0) THEN
     GO TO   240
   ELSE IF (key == 0) THEN
     GO TO   400
   END IF
!              EOR, EOF, KEY
   
   230 IF (key > lcor) GO TO 1510
   READ (in) (x(l),l=1,key)
   CALL WRITE (output,x,key,0)
   GO TO 220
   240 CALL WRITE (output,x,0,1)
   IF (iptx  == ipt4) GO TO 210
   IF (irecf ==    0) GO TO 200
   IF (irecf == 1 .OR. INDEX > 0) GO TO 250
   INDEX = 1
   GO TO 220
   
!     READ STRING FORMATTED MATRIX
   
   250 IF (irecf == 2 .AND. INDEX == 2) GO TO 260
   INDEX = 2
   CALL makmcb (mcb(1),output,trl(3),trl(4),trl(5))
   irow  = 1
   nrow  = trl(3)
   typin = trl(5)
   typout= trl(5)
   nwdsx = nwds(typout)
   ncol  = trl(2)
   
!     CHECK FOR NULL MATRIX
   
   IF (nrow == 0 .OR. ncol == 0) GO TO 400
   IF (irecf == 2) ncol = 1
   incr  = 1
   nwdsx = nrow*nwdsx
   260 keyx  = nwdsx
   
!     NWDSX IS NUMBER OF WORDS NEEDED PER COLUMN
   
   IF (sparse) GO TO 300
   DO  l = 1,ncol
     READ (in) key
     IF (key /= keyx) GO TO 1530
     IF (key > lcor) GO TO 1510
     READ (in) (x(k),k=1,nwdsx)
     CALL pack (x,output,mcb)
     READ (in) key
     imhere = 265
     IF (key > 0) GO TO 1560
   END DO
   280 IF (irecf == 2) GO TO 200
   keyx = 0
   READ (in) key
   imhere = 285
   IF (key /= keyx) GO TO 1530
   GO TO 400
   
!     SPARSE MATRIX INPUT (P5 = NON-ZERO)
!     (NOT CALLING FROM INPTT4 (IPTX=IPT2)
   
   300 DO  l = 1,ncol
     DO  k = 1,nwdsx
       x(k) = 0.0
     END DO
     320 READ (in) key,base
     IF (key < 0) GO TO 330
     READ (in) (x(k+base),k=1,key)
     GO TO 320
     330 CALL pack (x,output,mcb)
   END DO
   GO TO 280
   
!     MATRIX DATA BLOCK - MSC/STRING RECORD. (IPTX=IPT4)
   
   
   350 bflag = -1
   bcol  = bcol + 1
   360 READ (in) key
   CALL putstr (blk)
   imhere = 360
   IF (key < 0) THEN
     GO TO   390
   ELSE IF (key == 0) THEN
     GO TO  1560
   END IF
!       NULL or EOR, ERR, KEY
   
   370 bwrt = key/nwd
   imhere = 370
   IF (bwrt > brav) GO TO 1560
   
!     COMMENTS FROM G.C./UNISYS   3/93
!     UNLESS MSC/PUTSTR IS DIFFERENT FROM COSMIC/PUTSTR, THE FOLLOWING
!     3 LINES, ORIGINATED FROM MSC SOURCE CODE, DO NOT WORK FOR D.P.
!     DATA ON VAX, AND POSSIBLY SILICON-GRAHPICS. THEY ARE REPLACED BY
!     NEXT 17 LINES BELOW.
!     (I TRIED  SETTING L1=(BPOINT-1)*NWD+1, AND STILL DID NOT WORK.)
!     THE PROBLEM HERE IS D.P. DATA MAY FALL SHORT ON DOUBLE WORD
!     BOUNADRY, AND THEREFORE BECOME GARBAGE, WHICH MAY CAUSE FATAL
!     ERROR IN PRINTING.
   
!     L1 = BPOINT
!     L2 = L1 - 1 + KEY
!     READ (IN) BROW,(CORE(L),L=L1,L2)
   
   l1 = bpoint*nwd
   l2 = l1 - 1 + key
   IF (dp) GO TO 380
!     L  = 375
!     WRITE  (NOUT,375) L,L1,L2,KEY,BROW,BTYP,BPOINT
! 375 FORMAT (' /@',I3,'  L1,L2,KEY,BROW,BTYP,BPOINT =',4I7,4I4)
   READ   (in) brow,(core(l),l=l1,l2)
!     WRITE  (NOUT,376,ERR=388) (CORE(L),L=L1,L2)
! 376 FORMAT (10X,' CORE =',/,(1X,11E11.3))
   GO TO 385
   380 l1 = l1/2
   l2 = l2/2
!     L  = 382
!     WRITE  (NOUT,375) L,L1,L2,KEY,BROW,BTYP,BPOINT
   READ   (in) brow,(dcore(l),l=l1,l2)
!     WRITE  (NOUT,382,ERR=388) (DCORE(L),L=L1,L2)
! 382 FORMAT (10X,'DCORE =',/,(1X,11D11.3))
   385 CALL endput (blk)
   GO TO 360
   390 bflag = +1
   bwrt  =  0
   CALL endput (blk)
   GO TO 210
   
!     CLOSE OUTPUT DATA BLOCK WITH REWIND AND EOF.
   
   400 CALL CLOSE (output,1)
   
!     WRITE TRAILER.
   
   trl(1) = output
   CALL wrttrl (trl)
   CALL page2 (-3)
   WRITE  (nout,410) uim,NAME,in,namex
   410 FORMAT (a29,' 4105, DATA BLOCK ',2A4,' RETRIEVED FROM FORTRAN ',  &
       'TAPE ',i2, /5X,'ORIGINAL NAME OF DATA BLOCK WAS ',2A4)
   IF (sparse .AND. ntrl == 8 .AND. trl(8) /= 0)  &
       WRITE (nout,420) trl(2),trl(3)
   420 FORMAT (1H+,55X,'(A SPARSE MATRIX',i6,2H x,i6,')')
   
 END DO
 
!     CLOSE FORTRAN TAPE WITHOUT REWIND.
 
 440 CONTINUE
 RETURN
 
!     OBTAIN LIST OF DATA BLOCKS ON FORTRAN TAPE.
 
 500 REWIND in
 READ (in) key
 keyx = 3
 IF (key /= keyx) GO TO 1530
 READ (in) dx
 READ (in) key
 keyx = 7
 IF (key /= keyx) GO TO 1530
 READ (in) idhdrx
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 1460
 END DO
 READ (in) key
 keyx = 2
 IF (key /= keyx) GO TO 1530
 READ (in) p3x
 READ (in) key
 imhere = 515
 IF (key >= 0) GO TO 1560
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 1480
 520 ASSIGN 530 TO ret
 nskip = 1
 GO TO 1300
 530 kf   = 0
 540 CALL page1
 line = line + 8
 WRITE  (nout,550) in
 550 FORMAT (//50X,'FILE CONTENTS ON FORTRAN UNIT ',i2, /51X,32(1H-),  &
     //54X,4HFILE,18X,4HNAME,//)
 560 READ (in) key
 IF (key == 0) GO TO 590
!     KEYX = 2
!     IF (KEY .NE. KEYX) GO TO 9918
 READ (in) namex
!     READ (IN) KEY
!     IF (KEY .GE. 0) GO TO 9919
 ASSIGN 570 TO ret
 nskip = 1
 GO TO 1300
 570 kf   = kf + 1
 line = line + 1
 WRITE  (nout,580) kf,namex
 580 FORMAT (53X,i5,18X,2A4)
 IF (line - nlpp < 0) THEN
   GO TO   560
 ELSE
   GO TO   540
 END IF
 590 REWIND in
 ASSIGN 600 TO ret
 nskip = 1
 GO TO 1300
 600 CONTINUE
 GO TO 160
 
!     SEARCH MODE
 
!     EXAMINE OUTPUT REQUESTS AND FILL NAME TABLE.
 
 700 nnt = 0
 DO  i = 1,5
   output = 200 + i
   trl(1) = output
   CALL rdtrl (trl)
   IF (trl(1) <= 0) GO TO 710
   CALL fname (output,NAME)
   IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) GO TO 710
   nt(i,1) = 0
   nt(i,2) = NAME(1)
   nt(i,3) = NAME(2)
   nnt = nnt + 1
   CYCLE
   710 nt(i,1) = -1
!     IF (IPTX .NE. IPT2) GO TO 3050
   nt(i,2) = NONE(1)
   nt(i,3) = NONE(2)
 END DO
 
 IF (nnt > 0) GO TO 800
 CALL page2 (-2)
 WRITE  (nout,730) uwm,iptx
 730 FORMAT (a25,' 4137, ALL OUTPUT DATA BLOCKS FOR INPUTT',a1,  &
     ' ARE PURGED.')
!     CLOSE (UNIT=IN)
 RETURN
 
!     CHECK TAPE ID LABEL.
 
 800 REWIND in
 READ (in) key
 keyx = 3
 IF (key /= keyx) GO TO 1530
 READ (in) dx
 READ (in) key
 keyx = 7
 IF (key /= keyx) GO TO 1530
 READ (in) idhdrx
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 1460
 END DO
 READ (in) key
 keyx = 2
 IF (key /= keyx) GO TO 1530
 READ (in) p3x
 READ (in) key
 imhere = 815
 IF (key >= 0) GO TO 1560
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 1440
 ASSIGN 820 TO ret
 nskip = 1
 GO TO 1300
 820 CONTINUE
 
!     BEGIN SEARCH OF TAPE.
 
 kf = 0
 830 READ (in) key
 IF (key == 0) GO TO 1140
!     KEYX = 2
!     IF (KEY .NE. KEYX) GO TO 9918
 READ (in) namex
 READ (in) key
 imhere = 835
 IF (key >= 0) GO TO 1560
 kf = kf + 1
 
 DO  i = 1,5
   NAME(1) = nt(i,2)
   NAME(2) = nt(i,3)
   IF (nt(i,1) < 0) CYCLE
   IF (NAME(1) /= namex(1) .OR. NAME(2) /= namex(2)) CYCLE
   nt(i,1) = nt(i,1) + 1
   IF (nt(i,1) == 1 .OR. p1 == msix .OR. p1 == mete) GO TO 850
   CALL page2 (-3)
   WRITE  (nout,840) uwm,NAME,kf,in
   840 FORMAT (a25,' 4138, DATA BLOCK ,',2A4,' (DATA BLOCK COUNT =',i6,  &
       ')  HAS PREVIOUSLY BEEN RETRIEVED FROM ', /36X,  &
       'FORTRAN TAPE ',i2,' AND WILL BE IGNORED.')
   GO TO 1110
   850 READ (in) key
   keyx = ntrl
   IF (key /= keyx) GO TO 1530
   READ (in) (trl(l),l=1,ntrl)
   IF (iptx == ipt2) irecf = trl(8)
   READ (in) key
   imhere = 855
   IF (key >= 0) GO TO 1560
   
   output = 200 + i
   CALL OPEN (*1400,output,x(oubuf),1)
   
   INDEX = 0
!     IF (IPTX .EQ. IPT4) GO TO 890        ! FROM MSC/INPTT4
   READ (in) key
   IF (iptx == ipt2) GO TO 860
   keyx = 1
   IF (key == keyx) GO TO 1530
   READ (in) krec
   imhere = 857
   IF (krec < 0) GO TO 1560
   READ (in) key
   860 keyx = 2
   IF (key < keyx) GO TO 1530
   IF (key > lcor) GO TO 1510
   READ (in) (x(l),l=1,key)
   CALL WRITE (output,NAME,2,0)
   IF (key == keyx) GO TO 870
   CALL WRITE (output,x(i3),key-2,0)
   
   870 IF (iptx == ipt2) GO TO 890
   READ (in) key
   imhere = 875
   IF (key > 0) GO TO 1560
   btyp  = trl(5)
   bform = 0
   bcol  = 0
   nwd   = nwds(btyp)
   dp    = btyp == 2 .OR. btyp == 4
   CALL WRITE (output,0,0,1)
   880 READ (in) key
   keyx = 1
   IF (key /= keyx) GO TO 1530
   READ (in) krec
   IF (krec /= 0) GO TO 1010
   
!     TABLE DATA BLOCK(S)
   
   890 READ (in) key
   IF (key < 0) THEN
     GO TO   910
   ELSE IF (key == 0) THEN
     GO TO  1060
   END IF
!              EOR, EOF, KEY
   
   900 IF (key > lcor) GO TO 1510
   READ (in) (x(l),l=1,key)
   CALL WRITE (output,x,key,0)
   GO TO 890
   910 CALL WRITE (output,x,0,1)
!     IF (IPTX .EQ. IPT4) GO TO 890
   IF (iptx == ipt4) GO TO 880
   IF (irecf == 0) GO TO 870
   IF (irecf == 1) GO TO 920
   IF (INDEX > 0) GO TO 920
   INDEX = 1
   GO TO 870
   
!     READ STRING FORMATTED MATRIX
   
   920 IF (irecf == 2 .AND. INDEX == 2) GO TO 930
   INDEX = 2
   CALL makmcb (mcb(1),output,trl(3),trl(4),trl(5))
   irow  = 1
   nrow  = trl(3)
   typin = trl(5)
   typout= trl(5)
   nwdsx = nwds(typout)
   ncol  = trl(2)
   
!     CHECK FOR NULL MATRIX
   
   IF (nrow == 0 .OR. ncol == 0) GO TO 1060
   IF (irecf == 2) ncol = 1
   incr  = 1
   nwdsx = nrow*nwdsx
   930 keyx  = nwdsx
   
!     NWDSX IS NUMBER OF WORDS NEEDED PER COLUMN
   
   IF (sparse) GO TO 960
   DO  l = 1,ncol
     READ (in) key
     IF (key /= keyx) GO TO 1530
     IF (key > lcor) GO TO 1510
     READ (in) (x(k),k=1,nwdsx)
     CALL pack (x,output,mcb)
     READ (in) key
     imhere = 935
     IF (key > 0) GO TO 1560
   END DO
   950 IF (irecf == 2) GO TO 870
   keyx = 0
   READ (in) key
   IF (key /= keyx) GO TO 1530
   GO TO 1060
   
!     SPARSE MATRIX INPUT (P4 = NON-ZERO)
!     (NOT CALLING FROM INPTT4 (IPTX=IPT2)
   
   960 DO  l = 1,ncol
     DO  k  = 1,nwdsx
       x(k) = 0.0
     END DO
     980 READ (in) key,base
     IF (key < 0) GO TO 990
     READ (in) (x(k+base),k=1,key)
     GO TO 980
     990 CALL pack (x,output,mcb)
   END DO
   GO TO 950
   
!     MATRIX DATA BLOCK, IPTX = IPT4.  MSC/STRING RECORD
   
   
   1010 bflag = -1
   bcol  = bcol + 1
   1020 READ (in) key
   CALL putstr (blk)
   imhere = 1025
   IF (key < 0) THEN
     GO TO  1050
   ELSE IF (key == 0) THEN
     GO TO  1560
   END IF
!       NULL or EOR,  ERR, KEY
   
   1030 bwrt = key/nwd
   imhere = 1030
   IF (bwrt > brav) GO TO 1560
   
!     L1 = BPOINT
!     L2 = L1 - 1 + KEY
!     READ (IN) BROW,(CORE(L),L=L1,L2)
   
   l1 = bpoint*nwd
   l2 = l1 - 1 + key
   IF (dp) GO TO 1035
   READ (in) brow,(core(l),l=l1,l2)
   GO TO 1040
   1035 l1 = l1/2
   l2 = l2/2
   READ (in) brow,(dcore(l),l=l1,l2)
   1040 CALL endput (blk)
   GO TO 1020
   1050 bflag = +1
   bwrt  =  0
   CALL endput (blk)
   GO TO 880
   
!     CLOSE OUTPUT DATA BLOCK WITH REWIND AND EOF
   
   1060 CALL CLOSE (output,1)
   
!     WRITE TRAILER
   
   trl(1) = output
   CALL wrttrl (trl)
   CALL page2 (-2)
   WRITE  (nout,1070) uim,NAME,in,kf
   1070 FORMAT (a29,' 4139, DATA BLOCK ',2A4,' RETRIEVED FROM FORTRAN ',  &
       'TAPE ',i2,' (DATA BLOCK COUNT =',i6,1H))
   IF (nt(i,1) > 1) GO TO 1080
   nnt = nnt - 1
   GO TO 1130
   1080 WRITE  (nout,1090) uwm
   1090 FORMAT (a25,' 4140, SECONDARY VERSION OF DATA BLOCK HAS REPLACED',  &
       ' EARLIER ONE.')
   CALL page2 (-2)
   GO TO 1130
 END DO
 
 1110 ASSIGN 1120 TO ret
 nskip = 1
 GO TO 1300
 1120 CONTINUE
 1130 IF (nnt > 0 .OR. p1 == msix .OR. p1 == mete) GO TO 830
 GO TO 1200
 
 1140 IF (nnt <= 0) GO TO 1200
 CALL page2 (-7)
 IF (p1 == mfiv .OR. p1 == msix) GO TO 1160
 WRITE  (nout,1150) uwm
 1150 FORMAT (a25,' 4141, ONE OR MORE DATA BLOCKS NOT FOUND ON FORTRAN',  &
     ' TAPE.')
 GO TO 1170
 1160 WRITE (nout,1500) ufm
 1170 DO  i = 1,5
   IF (nt(i,1) /= 0) GO TO 1190
   WRITE  (nout,1180) nt(i,2),nt(i,3)
   1180 FORMAT (20X,21HNAME of DATA BLOCK = ,2A4)
   1190 IF (iptx == ipt4) GO TO 1200
 END DO
 IF (p1 == mfiv .OR. p1 == msix) GO TO 1600
 
 1200 ASSIGN 1210 TO ret
 nskip = -1
 GO TO 1300
 1210 CONTINUE
 RETURN
 
!     SIMULATION OF SKPFIL (IN,NSKIP)
 
 1300 IF (nskip < 0) THEN
   GO TO  1320
 ELSE IF (nskip == 0) THEN
   GO TO  1310
 ELSE
   GO TO  1330
 END IF
 1310 GO TO ret, (120,150,530,570,600,820,1120,1210)
 1320 REWIND in
 
!     NSKIP = COMPLEMENT OF NSKIP.
 
 1330 DO  ns = 1,nskip
   1340 READ (in) key
   IF (key < 0) THEN
     GO TO  1340
   ELSE IF (key == 0) THEN
     GO TO  1360
   END IF
!               EOR, EOF, KEY
   
   1350 IF (key > lcor) GO TO 1510
   READ (in) (x(l),l=1,key)
   GO TO 1340
   1360 CONTINUE
 END DO
 GO TO 1310
 
!     ERRORS
 
 1400 WRITE  (nout,1410) ufm,iptx,output
 1410 FORMAT (a23,' 4108, SUBROUTINE INPTT',a1,' UNABLE TO OPEN OUTPUT',  &
     ' DATA BLOCK',i6)
 GO TO  1600
 1420 WRITE  (nout,1430) ufm,iptx,p1
 1430 FORMAT (a23,' 4113, MODULE INPUTT',a1,' - ILLEGAL VALUE FOR ',  &
     'FIRST PARAMETER =',i20)
 GO TO  1600
 1440 WRITE  (nout,1450) ufm,p3x,iptx,p3
 1450 FORMAT (a23,' 4136, USER TAPE ID CODE -',2A4,'- DOES NOT MATCH ',  &
     'THIRD INPUTT',a1,' DMAP PARAMETER -',2A4,2H-.)
 line = line + 1
 GO TO  1600
 1460 WRITE  (nout,1470) ufm,iptx,idhdrx
 1470 FORMAT (a23,' 4134, MODULE INPUTT',a1,' - ILLEGAL TAPE CODE ',  &
     'HEADER = ',7A4)
 GO TO  1600
 1480 WRITE  (nout,1490) uwm,p3x,iptx,p3
 1490 FORMAT (a25,' 4135, USER TAPE ID CODE -',2A4,'- DOES NOT MATCH ',  &
     'THIRD INPUTT',a1,' DMAP PARAMETER -',2A4,2H-.)
 GO TO  520
 1500 FORMAT (a23,' 4142, ONE OR MORE DATA BLOCKS NOT FOUND ON USER ', 'TAPE')
 1510 WRITE  (nout,1520) ufm,lcor,key
 1520 FORMAT (a23,' 2187, INSUFFICIENT WORKING CORE TO HOLD FORTRAN ',  &
     'LOGICAL RECORD.', /5X,'LENGTH OF WORKING CORE =',i11,  &
     ',  LENGTH OF FORTRAN LOGICAL RECORD =',i11,1H.)
 line = line + 1
 GO TO  1600
 1530 WRITE  (nout,1540) sfm,key,keyx
 1540 FORMAT (a25,' 2190, ILLEGAL VALUE FOR KEY =',i10,  &
     ',   EXPECTED VALUE =',i11,1H.)
 IF (key == 2 .AND. keyx == 3) WRITE (nout,1550)
 1550 FORMAT (5X,'POSSIBLY DUE TO IMPROPER TAPE GENERATION PROCEDURE')
 GO TO  1600
 1560 WRITE  (nout,1570) sfm,key,imhere
 1570 FORMAT (a25,' 2190, ILLEGAL VALUE FOR KEY =',i10,'.  IMHERE =',i4)
 GO TO  1600
 1580 WRITE  (nout,1590) ufm,p3
 1590 FORMAT (a23,', ILLEGAL TAPE LABEL NAME -',2A4,'-  POSSIBLY ',  &
     'THE 4TH PARAMETER OF INPTT4 IS IN ERROR')
 GO TO  1600
 
 1600 line = line + 2
 CALL mesage (-61,lcor,subnam)
 RETURN
 
END SUBROUTINE inptt2
