SUBROUTINE xsort
     
!     SORT READS BULK DATA CARDS FROM THE INPUT TAPE, ADJUSTS THE
!     FIELDS, PERFORMS AN ALPHA-NUMERIC SORT ON THE CARD IMAGES FROM
!     LEFT TO RIGHT, INSERTS CONTINUATION CARDS IN THEIR PROPER
!     POSITION, AND PLACES THE RESULTING SORTED IMAGES ON THE NEW
!     PROBLEM TAPE.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: dec
 DIMENSION       headu(32),heads(32),headn(32),iblkda(2),cdcnt(3),  &
     bk(4),mk(4),ibuf1(20),ibuf2(20),ibuf3(2),  &
     kparnt(2),ibuf1a(2),ibuf2a(2),nsort(2),iiend(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /system/ ibufsz,outtap,nogo,intape,d1(14),iecho,d,iapprc,  &
     dum1(2),iuedit,dum44(44),isubs,dum12(12),icpflg, dum8(8),lpch
 COMMON /output/ dum2(96),head1(32),head2(32),head3(32)
 COMMON /zzzzzz/ skip1,buf(1)
 COMMON /xsrtcm/ bimsk1(6),bimsk2(5),bimsk3(4),bimsk4(4),bimsk5(2),  &
     bimsk6,bkmsk1(8),bkmsk2,shifts(4),  &
     icon1,icon2,star,plus,dollar,starl,slash,sftm, mask,BLANK,mka,is,mbit4
 COMMON /stapid/ krap(12),kumf
 COMMON /xechox/ ffflag,echou,echos,echop
 EQUIVALENCE     (bk(1),bkmsk1(5)),(mk(1),bimsk2(2)),  &
     (mkb  ,bimsk5(1)),(inf  ,bimsk2(1)), (sfta ,shifts(2)),(mkd  ,bimsk2(2)),  &
     (mke  ,bimsk5(2)),(mkc  ,bimsk4(1))
 EQUIVALENCE     (blanx,bkmsk1(8))
 DATA headu/10*4H    ,4H i n,4H p u,4H t  ,4H b u,4H l k,4H   d,  &
     4H a t,4H a  ,4H d e,4H c k,4H   e,4H c h,4H o  ,9*4H    /
 DATA heads/11*4H    ,4H s o,4H r t,4H e d,4H   b,4H u l,4H k  ,  &
     4H d a,4H t a,4H   e,4H c h,4H o  ,10*4H    /
 DATA headn/ 3*4H    ,4H    ,4H    ,4H    ,4H .  ,4H 1  ,4H..  ,  &
     4H 2  ,4H..  ,4H 3  ,4H..  ,4H 4  ,4H..  ,4H 5  ,4H..  ,  &
     4H 6  ,4H..  ,4H 7  ,4H..  ,4H 8  ,4H..  ,4H 9  ,4H..  ,  &
     4H10  ,4H.   ,5*4H    /
 DATA cdcnt/4HCARD,4HCOUN,4HT   /,nsort/4HXSOR,4HT   /
!     DATA BK/4H000 ,4H00  ,4H0   ,4H    /
!     DATA (MK(I),I=1,4)/O777777007777,O777700007777,O770000007777,O0/
!     DATA MKA,MKB,INF,SFTA/O000000777777,O377777777777,O777777777777,6/
!     DATA MKC/O007777777777/,MKD/O777777007777/,MKE/O377777007777/
 DATA iend1,iend2/4HENDD,4HATA /
 DATA iend3,iend4/4HENDA,4HTA  /
DATA iend5,iend6/4HEND ,4HDATA/
!     DATA STAR,PLUS,DOLLAR,STARL/4H000*,4H+000,4H$000,4H*000/
DATA iblkda/4HBULK,4HDATA/, iptp/4HOPTP/, nptp/4HNPTP/
DATA itape1,itape2,itape3,itape4,itape5/301,302,303,304,305/
DATA umf/4HUMF /
DATA idup,iok/4HDUPL, 4HOK  /


!     XSORT MAY NOT WORK PROPERLY IN ALL UNIX MACHINES, WHICH FOLLOW
!     THE VAX LINE.

dec    = mach == 5 .OR. mach == 6 .OR. mach == 21
IF (dec .AND. lpch /= 77) WRITE (outtap,5) uwm
5 FORMAT (a25,', SWITCHING TO OLD XSORT VIA DIAG 42 HAS NOT BEEN ',  &
    'THOROUGHLY TESTED', /5X,'FOR THE UNIX MACHINES.')

!     INITIALIZE XSORT AND TURN ON FREE-FIELD FLAG FOR XREAD

ffflag = 1234
echou  = 0
echos  = 0
echop  = 0
iend   = 0
iseq   = 0
iccbrk = 0
notsor = 0
optp   = iptp
kin    = 0
irestr = -iapprc
IF (kumf <= 0) GO TO 90
kin    = 1
CALL OPEN (*50,umf,buf(1),2)

!     FIND PARTICULAR BULK DATA FILE ON UMF AS REQUESTED BY USER

10 CALL READ (*30,*60,umf,pid,1,1,iflg)
IF (kumf-pid < 0) THEN
  GO TO    30
ELSE IF (kumf-pid == 0) THEN
  GO TO    80
END IF
20 CALL skpfil (umf,1)
GO TO 10
30 WRITE  (outtap,35) ufm,kumf
35 FORMAT (a23,' 201, REQUESTED BULK DATA DECK',i8,' NOT ON USER ',  &
    'MASTER FILE.')
CALL page2 (2)
nogo =-1
CALL CLOSE (umf,1)
RETURN

50 WRITE  (outtap,55) sfm
55 FORMAT (a25,' 202, UMF COULD NOT BE OPENED')
GO TO  1800
60 WRITE  (outtap,65) sfm
65 FORMAT (a25,' 203, ILLEGAL EOR ON UMF')
GO TO  1800
80 CALL CLOSE (umf,2)

90 CALL initco
IF (iecho < 0) GO TO 110
echou = andf(iecho,1)
echos = andf(iecho,2)
echop = andf(iecho,4)
IF (icpflg /= 0) echos = 1
110 ASSIGN 1260 TO ibrana
ASSIGN 810  TO ibranb
ASSIGN 1220 TO ibranf

!     SET ASSIGN GO TO SWITCHES FOR MACHINE CONFIGURATIONS
!     THE 8 BIT CHARACTER BYTE OF THE 360 WILL HOLD THE INTERNAL
!     CHARACTER CODE (MAX=37) WITHOUT USE OF THE 1ST BIT POSITION -
!     THE OTHER 3 MACHINES HAVE 6 BIT CHARACTERS THEREFORE A SHIFT RIGHT
!     OF ONE MUST BE DONE TO REMOVE A POSSIBLE BIT FROM THE SIGN
!     POSITION THE FOLLOWING ASSIGNS SET THOSE BRANCHES BASED ON MACHINE

IF (mach == 2 .OR. dec) GO TO 120
ASSIGN 350  TO mx3
ASSIGN 790  TO my1
ASSIGN 820  TO my2
ASSIGN 960  TO my3
ASSIGN 990  TO my4
ASSIGN 1030 TO my5
ASSIGN 840  TO my6
ASSIGN 730  TO mz1
linf   = 0
nshift = 1

!     SET NSHIFT TO ZERO FOR UNIVAC ASCII VERSION ONLY (NOT FORTRAN 5)

IF (mach == 3) nshift = 0
GO TO 130
120 CONTINUE
ASSIGN 360  TO mx3
ASSIGN 800  TO my1
ASSIGN 970  TO my2
ASSIGN 970  TO my3
ASSIGN 1040 TO my4
ASSIGN 1040 TO my5
ASSIGN 850  TO my6
ASSIGN 740  TO mz1
linf = orf(is,1)
130 CONTINUE

!     START WORKING SORT BUFFER BELOW GINO I/O BUFFERS

ii     = 5*ibufsz + 1
ibufbg = ii + 42
ibuflg = korsz(buf) - 21
IF (ibuflg-ibufbg < 210) CALL mesage (-8,ibufbg+210-ibuflg,nsort)
itape  = itape1
jtape  = itape2

!     OPEN ITAPE4 AND ITAPE5
!     (4 CONTAINS CONTINUATIONS, 5 CONTAINS ALTERS)

nbuf3 = 3*ibufsz + 1
CALL OPEN (*1700,itape4,buf(nbuf3),1)
nbuf4 = 4*ibufsz + 1
CALL OPEN (*1700,itape5,buf(nbuf4),1)

!     A BUFFER LINE IS 20 WORDS OF CARD IMAGE PLUS A 1 WORD POINTER TO
!     THE NEXT IMAGE IN THE SORT SEQUENCE - A ZERO POINTER INDICATES
!     THE LAST IMAGE (LARGEST IN SORT)
!     INITIALIZE WORKING BUFFER - 1ST LINE ZEROS, 2ND LINE ALL BITS

k = ii + 19
DO  j = ii,k
  buf(j) = linf
  buf(j+ 21) = inf
END DO
buf(ii+41) = 0

!     SET UP UNSORTED HEADING

DO  j = 1,32
  head1(j) = headu(j)
  head3(j) = headn(j)
END DO
head2(4) = headn(1)
iccnt    = 0
IF (echou == 0) GO TO 160
CALL page

!     OPEN ITAPE (LOCATION FOR EACH SORTED CORE LOAD AS ITS FORCED TO
!     EMPTY

160 CALL OPEN (*1700,itape,buf(1),1)
170 buf(ii+20) = 1
k    = ii
ncnt = 2

!     LOOP TO INPUT AND SORT CARD IMAGES - USES OPEN CORE FOR SORTED
!     IMAGES

DO  n1 = ibufbg,ibuflg,21
  n2 = n1 + 19
  n3 = n2 + 1
  180 CALL xread (*1770,buf(n1))
  iccnt = iccnt + 1
  IF (echou == 0) GO TO 220
  CALL page2 (-1)
  WRITE  (outtap,200)(buf(i),i=n1,n2)
  200 FORMAT (30X,20A4)
  210 FORMAT (13X,i8,1H-,8X,20A4)
  
!     IGNORE BLANK CARDS
  
  220 IF (buf(n1) == blanx .AND. buf(n1+1) == blanx) GO TO 180
  
!     LEFT ADJUST FIELD 1
  
  CALL xfadj1 (buf(n1),lshift,0)
  
!     TEST FOR END OF INPUT DATA STREAM (ENDDATA)
  
  iiend(1) = iend1
  iiend(2) = iend2
  IF (buf(n1) == iend1 .AND. buf(n1+1) == iend2) GO TO 560
  iiend(1) = iend3
  iiend(2) = iend4
  IF (buf(n1) == iend3 .AND. buf(n1+1) == iend4) GO TO 560
  iiend(1) = iend5
  iiend(2) = iend6
  IF (buf(n1) == iend5 .AND. buf(n1+1) == iend6) GO TO 560
  
!     IS THIS A CONTINUATION, COMMENT, OR DELETE CARD
  
  IF (.NOT.dec) tst = andf(mk(3),buf(n1))
  IF (     dec) tst = khrfn1(bkmsk2,1,buf(n1),1)
  
!     WRITE CONTINUATIONS ON ITAPE4
  
  IF (tst == starl .OR. tst == plus) GO TO 530
  
!     IGNORE COMMENT CARDS
  
  IF (tst == dollar) GO TO 180
  
!     WRITE DELETES ON ITAPE5
  
  IF (tst == slash) GO TO 540
  
!     IF A STAR IS FOUND IN FIELD 1, MOVE IT TO COLUMN 8
  
  ny  = 4
  DO  j = 1,2
    nx  = n1 + 2 - j
    tst = buf(nx)
    DO  i = 1,ny
      IF (.NOT.dec) ptst = andf(mka,tst)
      IF (     dec) ptst = khrfn1(bkmsk2,4,tst,4)
      IF (ptst /= bk(1)) GO TO 250
      IF (.NOT.dec) tst = rshift(tst,sfta)
      IF (     dec) tst = khrfn3(bkmsk2,tst,1,0)
    END DO
    ny = 3
  END DO
  GO TO 260
  
!     STARSW = 0 FOR A SINGLE FIELD CARD (NO STAR)
!            = 1 FOR A DOUBLE FIELD CARD (W/ STAR)
  
  250 starsw = 0
  IF (ptst /= star) GO TO 260
  starsw = 1
  IF (j == 1 .AND. i == 1) GO TO 260
  IF (dec) GO TO 258
  buf(nx  ) = orf(andf(mk(i),buf(nx)),bk(i))
  buf(n1+1) = orf(andf(mk(1),buf(n1+1)),star)
  GO TO 260
  258 buf(nx  ) = khrfn1(buf(nx),5-i,bk(i),5-i)
  buf(n1+1) = khrfn1(buf(n1+1),4,star,4)
  260 CONTINUE
  CALL xfadj (buf(n1+2),starsw,ny)
  CALL extint (buf(n1))
  
  
!     START SORT LOGIC
  
!     WITHOUT THE FOLLOWING CARD, XSORT WILL ASSUME SOME DEGREE OF SORT
!     EXISTS (I.E.,THE NEXT CARD WILL FOLLOW THE PREVIOUS CARD, MORE
!     OFTEN THAN NOT)
!     K  = II  (THIS CARD WILL FORCE SORT TO BEGINNING OF CHAIN)
  
  kp = 0
  
!     K TYPE SUBSCRIPTS REFER TO POSITIONS AND ITEMS IN THE SORTED
!     TABLE CURRENTLY BEING BUILT
!     N TYPE SUBSCRIPTS REFER TO ITEMS ABOUT THE NEWEST CARD IN
  
  270 fcnt = 1
  ni = 0
  ki = 0
  nx = n1
  
!     THE RIGHT SHIFT IN THE FOLLOWING CODE IS USED TO AVOID THE
!     NEGATIVE SIGN PROBLEM WHICH WOULD REVERSE THE SORT ORDER ON SOME
!     MACHINES.
!     (NOTE THAT THE SORT COMPARES CAN BE MADE BOTH WITH OR WITHOUT
!     THE SIGN SHIFT DEPENDING ON THE MACHINES CHARACTER CONFIG)
  
  300 kx = k
  GO TO 340
  330 IF (buf(nx) == buf(kx)) GO TO 400
  IF (buf(nx) == bk(4)  ) GO TO 380
  IF (buf(kx) == bk(4)  ) GO TO 370
  340 GO TO mx3, (350,360)
  350 IF (rshift(buf(nx),nshift)-rshift(buf(kx),nshift) < 0.0) THEN
    GO TO   380
  ELSE IF (rshift(buf(nx),nshift)-rshift(buf(kx),nshift) == 0.0) THEN
    GO TO   400
  ELSE
    GO TO   370
  END IF
  360 IF (dec) GO TO 365
  IF (buf(nx) < buf(kx)) GO TO 380
  IF (buf(nx) > buf(kx)) GO TO 370
  GO TO 400
  365 IF (rshift(khrfn4(buf(nx)),1)-rshift(khrfn4(buf(kx)),1)) 380,366,370
  366 IF (rshift(lshift(khrfn4(buf(nx)),1),1)-  &
      rshift(lshift(khrfn4(buf(kx)),1),1)) 380,400,370
  
!     GO ON, LOOK AT NEXT ITEM IN THE SORTED TABLE
  
  370 kp = k
  k  = buf(k+20)*21 + ii
  IF (nx == n1) GO TO 300
  GO TO 270
  
!     CARD POSITION FOUND IN SORT, SET THE CHAINING POINTER
  
  380 IF (kp == 0) GO TO 390
  buf(n3   ) = buf(kp+20)
  buf(kp+20) = ncnt
  k    = kp
  ncnt = ncnt + 1
  CYCLE
  390 k = ii
  GO TO 270
  
!     TWO FIELDS EQUAL - SLIDE TO NEXT FIELD ON CARD
  
  400 fcnt = fcnt + 1
  nx   = nx + 1
  kx   = kx + 1
  SELECT CASE ( fcnt )
    CASE (    1)
      GO TO 1760
    CASE (    2)
      GO TO 410
    CASE (    3)
      GO TO 470
    CASE (    4)
      GO TO 330
    CASE (    5)
      GO TO 510
    CASE (    6)
      GO TO 330
    CASE (    7)
      GO TO 430
    CASE (    8)
      GO TO 330
    CASE (    9)
      GO TO 510
    CASE (   10)
      GO TO 330
    CASE (   11)
      GO TO 520
    CASE (   12)
      GO TO 330
    CASE (   13)
      GO TO 510
    CASE (   14)
      GO TO 330
    CASE (   15)
      GO TO 430
    CASE (   16)
      GO TO 330
    CASE (   17)
      GO TO 510
    CASE (   18)
      GO TO 330
    CASE (   19)
      GO TO 380
  END SELECT
  410 ktarsw = 0
  IF (.NOT.dec) itst = andf(mka,buf(k+1))
  IF (     dec) itst = khrfn1(bkmsk2,4,buf(k+1),4)
  IF (itst == star) ktarsw = 1
  IF (starsw == ktarsw) GO TO 340
  
!     IF ONE MEMBER OF THE 2ND FIELD HAS A STAR AND THE OTHER DOES NOT,
!     DELETE STARS FOR THE COMPARE
  
  IF (dec) GO TO 415
  in1 = rshift(andf(mkd,buf(nx)),1)
  ik2 = rshift(andf(mkd,buf(kx)),1)
  GO TO 418
  415 in1 = rshift(khrfn4(khrfn1(buf(nx),4,bkmsk2,1)),1)
  ik2 = rshift(khrfn4(khrfn1(buf(kx),4,bkmsk2,1)),1)
  418 IF (in1 /= ik2) GO TO 428
  IF (dec) GO TO 420
  in1 = andf(mke,buf(nx))
  ik2 = andf(mke,buf(kx))
  GO TO 425
  420 in1 = rshift(lshift(khrfn4(khrfn1(buf(nx),4,bkmsk2,1)),1),1)
  ik2 = rshift(lshift(khrfn4(khrfn1(buf(kx),4,bkmsk2,1)),1),1)
  425 IF (in1 == ik2) GO TO 400
  428 IF (in1 < ik2) GO TO 380
  GO TO 370
  
!     INCREMENT FIELD LOCATIONS IF FIELD TYPES DID NOT MATCH
  
  430 IF (ni-ki < 0) THEN
    GO TO   450
  ELSE IF (ni-ki == 0) THEN
    GO TO   460
  END IF
  440 nx = nx + ni
  ni = 0
  GO TO 460
  450 kx = kx + ki
  ki = 0
  
!     ADJUST FIELDS RIGHT OR LEFT AS REQUIRED
  
  460 CALL xfadj (buf(nx),starsw,k1)
  CALL xfadj (buf(kx),ktarsw,k2)
  GO TO 480
  470 IF (starsw == ktarsw) GO TO 330
  k1 = 0
  k2 = 0
  IF (dec) GO TO 472
  IF (andf(mk(3),buf(nx)) /= bkmsk1(4)) k1 = 1
  IF (andf(mk(3),buf(kx)) /= bkmsk1(4)) k2 = 1
  GO TO 480
  472 IF (khrfn1(bkmsk2,1,buf(nx),1) /= bkmsk1(4)) k1 = 1
  IF (khrfn1(bkmsk2,1,buf(kx),1) /= bkmsk1(4)) k2 = 1
  480 IF (starsw-ktarsw < 0.0) THEN
    GO TO   500
  ELSE IF (starsw-ktarsw == 0.0) THEN
    GO TO   330
  END IF
  490 ni = 2
  IF (k1+k2 == 2) GO TO 330
  nx = nx + 2
  ni = 0
  GO TO 330
  500 ki = 2
  IF (k1+k2 == 2) GO TO 330
  kx = kx + 2
  ki = 0
  GO TO 330
  510 IF (starsw /= ktarsw) GO TO 430
  IF (starsw ==      0) GO TO 430
  GO TO 330
  520 IF (starsw == ktarsw) GO TO 430
  GO TO 380
  
!     CONTINUATION CARD - PUT ON ITAPE4
  
  530 CALL WRITE (itape4,buf(n1),20,1)
  GO TO 180
  
!     BULK DATA DELETE CARD - PUT ON ITAPE5
  
!     TEST FOR EXTRANEOUS DATA IN FIELD 1 OF DELETE CARD
!     AND WRITE OUT TO SCRATCH FILE
  
  540 IF (.NOT.dec) itst1 = andf(buf(n1),bimsk1(6))
  IF (     dec) itst1 = andf(buf(n1),bimsk1(1))
  itst2 = andf(buf(n1+1),mbit4)
  ibk3  = andf(bk(3),mbit4)
  ibk4  = andf(bk(4),mbit4)
  IF (itst1 == ibk3 .AND. itst2 == ibk4) GO TO 545
  CALL page2 (2)
  WRITE  (outtap,541) ufm
  541 FORMAT (a23,' 221, EXTRANEOUS DATA IN FIELD 1 OF BULK DATA ',  &
      'DELETE CARD.')
  nogo = -2
  545 CALL xfadj1 (buf(n1+2),rshift,0)
  CALL xbcdbi (buf(n1+2))
  CALL xfadj1 (buf(n1+4),rshift,0)
  CALL xbcdbi (buf(n1+4))
  buf(n1+4) = buf(n1+5)
  CALL WRITE (itape5,buf(n1+3),2,1)
  GO TO 180
  
!     END OF BIG SORT LOOP
  
END DO
GO TO 590


!     SET (ENDDATA) CARD FOUND FLAG

560 iend = -1
IF (echou /= 1) GO TO 572
CALL page2 (2)
WRITE  (outtap,570) iccnt
570 FORMAT (//24X,12HTOTAL count=,i5)
572 CONTINUE

!     TEST FOR COLD-START WITH NO BULK DATA

IF (iccnt > 1 .OR. irestr > 0 .OR. kumf > 0) GO TO 590
IF (iapprc == 1) GO TO 590
IF (isubs  /= 0) GO TO 590
CALL page2 (2)
WRITE  (outtap,580) ufm
580 FORMAT (a23,' 204, COLD START NO BULK DATA.')
nogo = -2
RETURN


!     IF MODIFIED RESTART - TURN ON SORT ECHO

590 CONTINUE

!     THIS SECTION UNCHAINS THE SORTED TABLE AND WRITES A CORE LOAD,
!     IN ITS ACTUAL ORDER, ONTO A MERGE SCRATCH TAPE.

j = buf(ii+20)
j1st = j*21 + ii
KEEP = 1
610 j  = j*21 + ii
j1 = buf(j+20)
IF (j1 == 0) GO TO 620

!     ITAPE IS PRIMARY CORE UNLOAD TAPE

CALL WRITE (itape,buf(j),20,1)
IF (j < KEEP) notsor = 1
KEEP = j
j    = j1
GO TO 610
620 iseq = iseq + 1
IF (iseq == 2) GO TO 640
IF (iseq > 2) GO TO 650
IF (iend /= 0) GO TO 630
itape = itape2
CALL OPEN (*1700,itape,buf(ibufsz+1),1)
GO TO 170

!     NO MERGING IS REQUIRED, ALL CARDS FIT WITHIN ONE WORKING BUFFER
!     LOAD

630 ktape = itape
CALL CLOSE (ktape,1)
GO TO 1260

!     SET UP 1ST MERGE

640 CALL CLOSE (itape,1)
itape = itape1
kop   = 0

!     SET UP SUBSEQUENT MERGE TAPES

650 IF (MOD(iseq,2) == 0) GO TO 660
jtape = itape3
ktape = itape2
GO TO 670
660 jtape = itape2
ktape = itape3
670 CALL CLOSE (itape,1)

!     SPECIAL LOGIC TO AVOID MERGE IF NEW CORE LOAD FOLLOWS ALL PREVIOUS

IF (kop-1 < 0) THEN
  GO TO   760
ELSE IF (kop-1 == 0) THEN
  GO TO   700
END IF
680 DO  j = 1,18
  ibuf1(j) = ibuf2(j)
END DO
700 DO  j = 1,18
  IF (buf(j1st) /= ibuf1(j)) GO TO 720
  j1st = j1st + 1
END DO
GO TO 760
720 GO TO mz1, (730,740)
730 IF (rshift(buf(j1st),nshift) < rshift(ibuf1(j),nshift)) GO TO 755
GO TO 750
740 IF (dec) GO TO  745
IF (buf(j1st) < ibuf1(j)) GO TO 755
GO TO 750
745 IF (khrfn4(buf(j1st)) < khrfn4(ibuf1(j))) GO TO 755
750 trial = ktape
ktape = jtape
jtape = trial
iseq  = iseq - 1
CALL OPEN (*1700,itape,buf(1),0)
CALL OPEN (*1700,ktape,buf(ibufsz+1),3)
GO TO 1210

!     THIS SECTION PERFORMS A 2 TAPE ALPHANUMERIC MERGE
!     (ITAPE+JTAPE=KTAPE)
!     SAME BASIC LOGIC AS ORIGINAL SORT COMPARES (COMMENT CARDS OMITTED)

755 notsor = 1
760 CALL OPEN (*1700,itape,buf(1),0)
770 CALL OPEN (*1700,jtape,buf(ibufsz+1),0)
nbuf2 = 2*ibufsz + 1
CALL OPEN (*1700,ktape,buf(nbuf2),1)
ccnt = 0
780 CALL READ (*1190,*1710,jtape,ibuf2,20,1,iflg)
IF (mach == 2 .AND. (jtape == umf .OR. jtape == iptp)) CALL umftrn (ibuf2)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
ldup = 0
IF (itape == optp) CALL crdflg (ibuf2)
ktarsw = 0
IF (.NOT.dec) itst = andf(mka,ibuf2(2))
IF (     dec) itst = khrfn1(bkmsk2,4,ibuf2(2),4)
IF (itst == star) ktarsw = 1
GO TO my1, (790,800)
790 ibuf2a(1) = rshift(ibuf2(1),nshift)
ibuf2a(2) = rshift(ibuf2(2),nshift)
800 CALL READ (*1240,*1710,itape,ibuf1,20,1,iflg)
IF (mach == 2 .AND. (itape == umf .OR. itape == iptp)) CALL umftrn (ibuf1)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf1)
starsw = 0
IF (.NOT.dec) itst = andf(mka,ibuf1(2))
IF (     dec) itst = khrfn1(bkmsk2,4,ibuf1(2),4)
IF (itst == star) starsw = 1
GO TO ibranb, (830,810)
810 GO TO my2, (820,970)
820 ibuf1a(1) = rshift(ibuf1(1),nshift)
ibuf1a(2) = rshift(ibuf1(2),nshift)
GO TO 970

!     TEST IF CARD IS TO BE DELETED

830 ccnt = ccnt + 1
IF (.NOT.dec) tst = andf(mk(3),ibuf1(1))
IF (     dec) tst = khrfn1(bkmsk2,1,ibuf1(1),1)
iccflg = -1
IF (tst == plus .OR. tst == starl) GO TO 860
CALL extint (ibuf1(1))
GO TO my6, (840,850)
840 ibuf1a(1) = rshift(ibuf1(1),nshift)
ibuf1a(2) = rshift(ibuf1(2),nshift)
850 iccflg    = 0
kparnt(1) = ibuf1(1)
kparnt(2) = ibuf1(2)
860 GO TO ibranc, (870,880,900)
870 CALL READ (*920,*1710,itape5,ibuf3,2,1,iflg)
IF (ibuf3(1) == 0) GO TO 870
ASSIGN 880 TO ibranc
880 IF (ibuf3(2) /=    0) GO TO 890
IF (ibuf3(1) /= ccnt) GO TO 900
ASSIGN 870 TO ibranc
GO TO 930
890 IF (ibuf3(2) == ccnt) ASSIGN 870 toibranc
IF (ibuf3(1) <= ccnt .AND. ibuf3(2) >= ccnt) GO TO 930

!     REMOVE ANY UNDELETED CONTINUATION CARDS DURING RESTART MERGE

900 IF (iccflg == 0) GO TO ibrane, (970,1220)
CALL WRITE (itape4,ibuf1(1),20,1)
910 GO TO ibrand, (800,1210)
920 ASSIGN 900 TO ibranc
CALL CLOSE (itape5,1)
GO TO 900

!     IF CONTINUATION WAS DELETED, FLAG PARENT

930 IF (iccflg == 0) GO TO 940
CALL crdflg (kparnt)
GO TO 910
940 CALL crdflg (ibuf1)
GO TO 910
950 CALL READ (*1190,*1710,jtape,ibuf2,20,1,iflg)
IF (mach == 2 .AND. (jtape == umf .OR. jtape == iptp)) CALL umftrn (ibuf2)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
IF (itape == optp) CALL crdflg (ibuf2)
ktarsw = 0
IF (.NOT.dec) itst = andf(mka,ibuf2(2))
IF (     dec) itst = khrfn1(bkmsk2,4,ibuf2(2),4)
IF (itst == star) ktarsw = 1
GO TO my3, (960,970)
960 ibuf2a(1) = rshift(ibuf2(1),nshift)
ibuf2a(2) = rshift(ibuf2(2),nshift)
970 j  = 1
j1 = 1
j2 = 1
ni = 0
ki = 0
980 GO TO my4, (990,1040)
990 IF (ibuf1a(j1)-ibuf2a(j2) < 0) THEN
  GO TO  1050
ELSE IF (ibuf1a(j1)-ibuf2a(j2) == 0) THEN
  GO TO  1070
ELSE
  GO TO  1060
END IF
1000 IF (starsw == ktarsw) GO TO 980
IF (dec) GO TO 1005
in1 = rshift(andf(mkd,ibuf1(j1)),1)
ik2 = rshift(andf(mkd,ibuf2(j2)),1)
GO TO 1008
1005 in1 = rshift(khrfn4(khrfn1(ibuf1(j1),4,bkmsk2,1)),1)
ik2 = rshift(khrfn4(khrfn1(ibuf2(j2),4,bkmsk2,1)),1)
1008 IF (in1 /= ik2) GO TO 1018
IF (dec) GO TO 1010
in1 = andf(mke,ibuf1(j1))
ik2 = andf(mke,ibuf2(j2))
GO TO 1015
1010 in1 = rshift(lshift(khrfn4(khrfn1(ibuf1(j1),4,bkmsk2,1)),1),1)
ik2 = rshift(lshift(khrfn4(khrfn1(ibuf2(j2),4,bkmsk2,1)),1),1)
1015 IF (in1 == ik2) GO TO 1070
1018 IF (in1 < ik2) GO TO 1050
GO TO 1060
1020 IF (ibuf1(j1) == ibuf2(j2)) GO TO 1070
IF (ibuf1(j1) == bk(4)    ) GO TO 1050
IF (ibuf2(j2) == bk(4)    ) GO TO 1060
GO TO my5, (1030,1040)
1030 IF (rshift(ibuf1(j1),1)-rshift(ibuf2(j2),1) < 0.0) THEN
  GO TO  1050
ELSE IF (rshift(ibuf1(j1),1)-rshift(ibuf2(j2),1) == 0.0) THEN
  GO TO  1070
ELSE
  GO TO  1060
END IF
1040 IF (dec) GO TO 1045
IF (ibuf1(j1) < ibuf2(j2)) GO TO 1050
IF (ibuf1(j1) > ibuf2(j2)) GO TO 1060
GO TO 1070
1045 IF (khrfn4(ibuf1(j1))-khrfn4(ibuf2(j2)) < 0) THEN
  GO TO  1050
ELSE IF (khrfn4(ibuf1(j1))-khrfn4(ibuf2(j2)) == 0) THEN
  GO TO  1070
ELSE
  GO TO  1060
END IF
1050 CALL WRITE (ktape,ibuf1,20,1)
kop = 1
GO TO 800
1060 CALL WRITE (ktape,ibuf2,20,1)
kop = 2
GO TO 950
1070 j  = j  + 1
j1 = j1 + 1
j2 = j2 + 1
SELECT CASE ( j )
  CASE (    1)
    GO TO 1760
  CASE (    2)
    GO TO 1000
  CASE (    3)
    GO TO 1120
  CASE (    4)
    GO TO 1020
  CASE (    5)
    GO TO 1160
  CASE (    6)
    GO TO 1020
  CASE (    7)
    GO TO 1080
  CASE (    8)
    GO TO 1020
  CASE (    9)
    GO TO 1160
  CASE (   10)
    GO TO 1020
  CASE (   11)
    GO TO 1170
  CASE (   12)
    GO TO 1020
  CASE (   13)
    GO TO 1160
  CASE (   14)
    GO TO 1020
  CASE (   15)
    GO TO 1080
  CASE (   16)
    GO TO 1020
  CASE (   17)
    GO TO 1160
  CASE (   18)
    GO TO 1020
  CASE (   19)
    GO TO 1180
END SELECT
1080 IF (ni-ki < 0) THEN
  GO TO  1100
ELSE IF (ni-ki == 0) THEN
  GO TO  1110
END IF
1090 j1 = j1 + ni
ni = 0
GO TO 1110
1100 j2 = j2 + ki
ki = 0
1110 CALL xfadj (ibuf1(j1),starsw,k1)
CALL xfadj (ibuf2(j2),ktarsw,k2)
GO TO 1130
1120 IF (starsw == ktarsw) GO TO 1020
k1 = 0
k2 = 0
IF (dec) GO TO 1122
IF (andf(mk(3),ibuf1(j1)) /= bkmsk1(4)) k1 = 1
IF (andf(mk(3),ibuf2(j2)) /= bkmsk1(4)) k2 = 1
GO TO 1130
1122 IF (khrfn1(bkmsk2,1,ibuf1(j1),1) /= bkmsk1(4)) k1 = 1
IF (khrfn1(bkmsk2,1,ibuf2(j2),1) /= bkmsk1(4)) k2 = 1
1130 IF (starsw-ktarsw < 0.0) THEN
  GO TO  1150
ELSE IF (starsw-ktarsw == 0.0) THEN
  GO TO  1020
END IF
1140 ni = 2
IF (k1+k2 == 2) GO TO 1020
j1 = j1 + 2
ni = 0
GO TO 1020
1150 ki = 2
IF (k1+k2 == 2) GO TO 1020
j2 = j2 + 2
ki = 0
GO TO 1020
1160 IF (starsw /= ktarsw) GO TO 1080
IF (starsw ==      0) GO TO 1080
GO TO 1020

!     DUPLICATE CARD

1170 IF (starsw == ktarsw) GO TO 1080
1180 CALL WRITE (ktape,ibuf1,20,1)
CALL WRITE (ktape,ibuf2,20,1)
ldup = -1
GO TO 780

!     ONE OF TWO TAPES BEING MERGED IS EXHAUSTED, OTHER TAPE IS COPIED
!     ONTO THE MERGE TAPE

1190 IF (itape /= optp) GO TO 1200
ASSIGN 1210 TO ibrand
ASSIGN 1220 TO ibrane
ASSIGN 830 TO ibranf
IF (ccnt == 0) GO TO 1210
1200 IF (ldup < 0) GO TO 1210
GO TO 1220
1210 CALL READ (*1250,*1710,itape,ibuf1,20,1,iflg)
IF (mach == 2 .AND. (itape == umf .OR. itape == iptp)) CALL umftrn (ibuf1)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf1)
GO TO ibranf, (830,1220)
1220 CALL WRITE (ktape,ibuf1,20,1)
kop = 1
GO TO 1210
1230 CALL READ (*1250,*1710,jtape,ibuf2,20,1,iflg)
IF (mach == 2 .AND. (jtape == umf .OR. jtape == iptp)) CALL umftrn (ibuf2)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
1240 CALL WRITE (ktape,ibuf2,20,1)
kop = 2
GO TO 1230

!     CLOSE TAPES ENVOLVED IN MERGE

1249 CALL CLOSE (itape,2)
GO TO 1251
1250 IF (iuedit ==   1) GO TO 1249
IF (itape == optp) GO TO 1249
CALL CLOSE (itape,1)
1251 CALL CLOSE (jtape,3)
CALL CLOSE (ktape,3)
GO TO ibrana, (1260,1440)

!     WAS THIS THE FINAL MERGE (LAST CORE LOAD OF CARDS)

1260 IF (iend == 0) GO TO 160
CALL page2 (2)
WRITE (outtap,1660)
CALL CLOSE (itape5,1)

!     PROCESS DELETE CARDS (IF ANY)

nbuf4 = 4*ibufsz + 1
CALL OPEN (*1700,itape5,buf(nbuf4),0)

!     IF NOT RESTART - NO DELETES SHOULD EXIST

IF (irestr > 0 .OR. kin > 0) GO TO 1280

CALL READ (*1440,*1710,itape5,ibuf3,1,1,iflg)

!     NOT RESTART AND DELETES DO EXIST - WARNING

CALL CLOSE (itape5,1)
CALL page2 (2)
WRITE  (outtap,1270) uwm
1270 FORMAT (a25,' 205, COLD START,DELETE CARDS IGNORED.')
GO TO 1440

!     FORM DELETE CARD LIST

1280 ibuf3(1)  = inf
buf(ii  ) = mkb
buf(ii+1) = mkb
DO  j = ii,ibuflg,2
  CALL READ (*1330,*1710,itape5,ibuf3,2,1,iflg)
  DO  i = ii,j,2
    IF (ibuf3(1) <= buf(i)) GO TO 1300
  END DO
  
!     PUSH DOWN LIST - MAKE DOUBLE WORD SLOT
  
  1300 kk = j + 2
  k1 = (j-i)/2 + 1
  DO  k = 1,k1
    buf(kk+1) = buf(kk-1)
    buf(kk  ) = buf(kk-2)
    kk = kk - 2
  END DO
  buf(i  ) = ibuf3(1)
  buf(i+1) = ibuf3(2)
END DO

!     IF DELETE CARD LIST WILL NOT FIT

CALL mesage (-8,0,nsort)

!     EOF ON ITAPE5, IF IBUF3(1)= INF, THERE ARE NO DELETE CARDS

1330 IF (ibuf3(1) == inf) GO TO 1400
j = j - 1

!     CHECK FOR AND ELIMINATE OVERLAPS AND REDUNDANCYS IN DELETES

imin = 0
DO  i = ii,j,2
  IF (buf(i) == 0) CYCLE
  IF (buf(i) < buf(i+1)) GO TO 1340
  buf(i+1) = 0
  IF (buf(i) == buf(i+2)) GO TO 1350
  IF (imin == 0) CYCLE
  IF (buf(i) > buf(imax)) GO TO 1370
  GO TO 1350
  1340 IF (imin == 0) GO TO 1360
  IF (buf(i  ) > buf(imax)) GO TO 1360
  IF (buf(i+1) < buf(imax)) GO TO 1350
  buf(imax) = buf(i+1)
  1350 buf(i) = 0
  CYCLE
  1360 imin = i
  imax = i + 1
  CYCLE
  1370 imin = 0
END DO
CALL CLOSE (itape5,1)

!     PUT OUT SORTED DELETE CARD LIST

nbuf4 = 4*ibufsz + 1
CALL OPEN (*1700,itape5,buf(nbuf4),1)
DO  i = ii,j,2
  IF (buf(i) == 0) CYCLE
  CALL WRITE (itape5,buf(i),2,1)
END DO
1400 CALL CLOSE (itape5,1)

!     AT THIS POINT, IF THIS IS A RESTART, MERGE OPTP, FINAL KTAPE,
!     + DELETE

ASSIGN 1440 TO ibrana
ASSIGN 830  TO ibranb
ASSIGN 870  TO ibranc
ASSIGN 800  TO ibrand
ASSIGN 970  TO ibrane
nbuf4 = 4*ibufsz + 1
CALL OPEN (*1700,itape5,buf(nbuf4),0)

IF (kin > 0) GO TO 1430

CALL OPEN (*1740,optp,buf(1),0)
1410 CALL READ (*1730,*1710,optp,ibuf3,2,1,iflg)
IF (ibuf3(1) == iblkda(1) .AND. ibuf3(2) == iblkda(2)) GO TO 1420
CALL skpfil (optp,+1)
GO TO 1410
1420 itape = optp
trial = jtape
jtape = ktape
ktape = trial
IF (iccnt == 1) GO TO 1440
CALL WRITE (itape4,mkb,20,1)
GO TO 770

1430 optp = umf
CALL OPEN (*50,umf,buf(1),2)
GO TO 1420

!     PROCESS CONTINUATION CARDS (IF ANY)

1440 CALL CLOSE (itape4,1)
nbuf3 = 3*ibufsz + 1
CALL OPEN (*1700,itape4,buf(nbuf3),0)
IF (iccnt == 1 .AND. (irestr > 0 .OR. kin > 0)) ktape = optp
IF (iccnt == 1 .AND. (irestr > 0 .OR. kin > 0)) GO TO 1441
nbuf2 = 2*ibufsz + 1
CALL OPEN (*1700,ktape,buf(nbuf2),0)

!     FORM CONTINUATION CARD DICTIONARY

1441 CONTINUE
ibuf1(1) = 0
DO  j = ii,ibuflg,4
  CALL READ (*1480,*1710,itape4,ibuf1,20,1,iflg)
  IF (mach == 2 .AND. ibuf1(1) /= mkb) CALL umftrn (ibuf1)
  IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf1)
  IF (ibuf1(1) /= mkb) GO TO 1460
  IF (j == ii) GO TO 1450
  iccbrk = j
  1450 buf(j) = dollar
  GO TO 1465
  1460 IF (.NOT.dec) buf(j) = andf(mkc,ibuf1(1))
  IF (     dec) buf(j) = khrfn1(ibuf1(1),1,bkmsk2,1)
  1465 buf(j+1) = ibuf1(2)
  IF (.NOT.dec) buf(j+2) = andf(mkc,ibuf1(19))
  IF (     dec) buf(j+2) = khrfn1(ibuf1(19),1,bkmsk2,1)
  buf(j+3) = ibuf1(20)
END DO


!     CORE INSUFFICIENT TO ACCOMMODATE 4-WORD PER CARD DICTIONARY
!     OF CONTINUATION CARDS

CALL mesage (-8,0,nsort)

!     EOF ON ITAPE4, IF IBUF1(1)= 0, THERE ARE NO CONTINUATION CARDS

1480 IF (ibuf1(1) == 0) GO TO 1510
CALL REWIND (itape4)
jo = 1
iconlg = j - 1

!     CHECK AND SET FLAGS FOR DUPLICATE CONTINUATION CARDS

k = iconlg - 4
IF (k <= ii) GO TO 1510
DO  j = ii,k,4
  IF (buf(j) == idup) CYCLE
  INDEX = 0
  m = j + 4
  DO  jj = m,iconlg,4
    IF (buf(jj ) == idup     ) CYCLE
    IF (buf(j  ) /= buf(jj)  ) CYCLE
    IF (buf(j+1) /= buf(jj+1)) CYCLE
    buf(jj) = idup
    INDEX = 1
  END DO
  IF (INDEX == 1) buf(j) = idup
END DO

!     SET UP AND PUT OUT SORTED HEADING

1510 IF (notsor == 0) GO TO 1515
CALL page2 (2)
WRITE  (outtap,1511) uim
1511 FORMAT (a29,' 207, BULK DATA NOT SORTED, XSORT WILL RE-ORDER ', 'DECK.')
1515 IF (echos == 0) GO TO 1530
DO  j = 1,32
  head1(j) = heads(j)
END DO
head2(4) = cdcnt(1)
head3(4) = cdcnt(2)
head3(5) = cdcnt(3)
CALL page
ccnt = 0
1530 CALL CLOSE (itape5,1)
j = ii
nbuf4 = 4*ibufsz + 1
CALL OPEN (*1750,nptp,buf(nbuf4),3)
CALL WRITE (nptp,iblkda,2,1)
IF (ibuf1(1) == 0) GO TO 1630

!     MERGE CONTINUATION CARDS - PRODUCE DATA ON NPTP

1540 CALL READ (*1640,*1710,ktape,ibuf1,20,1,iflg)
IF (iccbrk == 0) GO TO 1550
kparnt(1) = ibuf1(1)
kparnt(2) = ibuf1(2)
1550 CALL intext (ibuf1(1))
IF (mach == 2) CALL umftrn (ibuf1)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf1)
CALL WRITE (nptp,ibuf1,20,1)
IF (echos == 0) GO TO 1551
CALL page2 (-1)
ccnt = ccnt + 1
CALL xprety (ibuf1)
WRITE (outtap,210) ccnt,ibuf1

!      PUNCH OUT DECK

1551 IF (echop == 0) GO TO 1554
IF (echos /= 0) GO TO 1552
CALL xprety (ibuf1)
1552 WRITE  (lpch,1553) ibuf1
1553 FORMAT (20A4)
1554 CONTINUE

!     SEE IF PREVIOUS CARD HAS A CONTINUATION
!     IF CONTINUATION FIELD BLANK - CONTINUATION NOT POSSIBLE

IF (ibuf1(19) == bk(4) .AND. ibuf1(20) == bk(4)) GO TO 1540
IF (.NOT.dec) trial = andf(mkc,ibuf1(19))
IF (     dec) trial = khrfn1(ibuf1(19),1,bkmsk2,1)
jn = 0
1571 CONTINUE

DO  k = ii,iconlg,4
  
!     IGNORE DUPLICATE CONTINUATION CARDS
  
  IF (buf(j) == idup) GO TO 1600
  IF (ibuf1(20) /= buf(j+1)) GO TO 1600
  IF (.NOT.dec) itst = andf(mkc,buf(j))
  IF (     dec) itst = khrfn1(buf(j),1,bkmsk2,1)
  IF (itst /= trial) GO TO 1600
  
!     A CONTINUATION EXISTS, HAS IT ALREADY BEEN USED
  
  IF (.NOT.dec) itst = andf(mk(3),buf(j))
  IF (     dec) itst = khrfn1(bkmsk2,1,buf(j),1)
  IF (itst == dollar) GO TO 1610
  IF (j > iccbrk) GO TO 1580
  CALL crdflg (kparnt)
  1580 IF (.NOT.dec) buf(j) = orf(buf(j),dollar)
  IF (     dec) buf(j) = khrfn1(buf(j),1,dollar,1)
  jn = (j-ii)/4 + 1
  CALL xrecps (jn,jo)
  CALL READ (*1720,*1710,itape4,ibuf1,20,1,iflg)
  IF (mach == 2) CALL umftrn (ibuf1)
  IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf1)
  CALL WRITE (nptp,ibuf1,20,1)
  IF (echos == 0) GO TO 1581
  CALL page2 (-1)
  ccnt = ccnt+ 1
  WRITE (outtap,210) ccnt,ibuf1
  1581 IF (echop == 0) GO TO 1584
  WRITE (lpch,1553) ibuf1
  1584 CONTINUE
  IF (.NOT.dec) trial = andf(mkc,ibuf1(19))
  IF (     dec) trial = khrfn1(ibuf1(19),1,bkmsk2,1)
  IF (ibuf1(19) == bk(4) .AND. ibuf1(20) == bk(4)) GO TO 1540
  GO TO 1571
  1600 j = j + 4
  IF (j > iconlg) j = ii
END DO
GO TO 1540

!     DUPLICATE PARENT - ERROR

1610 nl = 0
IF (echos /= 0) GO TO 1612
nl = 1
WRITE (outtap,200) ibuf1
1612 nl = nl +2
CALL page2 (-nl)
WRITE  (outtap,1620) ufm
1620 FORMAT (a23,' 208, PREVIOUS CARD IS A DUPLICATE PARENT.')
nogo = -1
GO TO 1540

!     NO CONTINUATION CARDS

1630 CALL READ (*1640,*1710,ktape,ibuf2,20,1,iflg)
IF (iccnt == 1) GO TO 1631
CALL intext (ibuf2(1))
1631 CONTINUE
IF (mach == 2) CALL umftrn (ibuf2)
IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
CALL WRITE (nptp,ibuf2,20,1)
IF (echos == 0) GO TO 16311
CALL page2 (-1)
ccnt = ccnt + 1
CALL xprety (ibuf2)
WRITE (outtap,210) ccnt,ibuf2
16311 IF (echop == 0) GO TO 1630
IF (echos /= 0) GO TO 1632
CALL xprety (ibuf2)
1632 WRITE (lpch,1553) ibuf2
GO TO 1630

!     CLOSE KTAPE AND WRITE (ENDDATA)

1640 CALL CLOSE (ktape,2)
CALL eof   (nptp)
CALL CLOSE (nptp,1)
IF (echos == 0) GO TO 1650
CALL page2 (-1)
WRITE (outtap,200) iiend
1650 IF (ibuf1(1) == 0) GO TO 1690
CALL page2 (2)
WRITE  (outtap,1660)
1660 FORMAT (1H0)

!     IDENTIFY DUPLICATE OR PARENTLESS CONTINUATION CARDS

ncnt = 0
DO  j = ii,iconlg,4
  IF (.NOT.dec) itst = andf(mk(3),buf(j))
  IF (     dec) itst = khrfn1(bkmsk2,1,buf(j),1)
  IF (itst == dollar) CYCLE
  
!     CHECK FOR DUPLICATE CONTINUATION CARDS
  
  IF (buf(j) == idup) GO TO 1666
  
!     CHECK FOR PARENTLESS CONTINUATION CARDS
  
  DO  jj = ii,iconlg,4
    IF (j == jj) CYCLE
    IF (buf(j) == buf(jj+2) .AND. buf(j+1) == buf(jj+3)) GO TO 1668
  END DO
  1666 ncnt = ncnt + 1
  jn   = (j-ii)/4 + 1
  CALL xrecps (jn,jo)
  CALL READ (*1720,*1710,itape4,ibuf2,20,1,iflg)
  IF (mach == 2) CALL umftrn (ibuf2)
  IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
  CALL page2 (-1)
  WRITE (outtap,200) ibuf2
  CYCLE
  1668 buf(j) = iok
END DO
IF (ncnt == 0) GO TO 1690
CALL page2 (3)
WRITE  (outtap,1680) ufm,ncnt
1680 FORMAT (a23,' 209, PREVIOUS',i7,' CONTINUATION MNEMONICS HAVE NO',  &
    ' PARENTS AND/OR ARE DUPLICATES.',/)
nogo = -1

!     IDENTIFY THOSE CONTINUATION CARDS THAT ARE VALID, BUT YET CANNOT
!     BE PROCESSED BECAUSE OF ERRORS ON OTHER RELATED CONTINUATION CARDS

ncnt = 0
DO  j = ii,iconlg,4
  IF (buf(j) /= iok) CYCLE
  ncnt = ncnt + 1
  jn   = (j-ii)/4 + 1
  CALL xrecps (jn,jo)
  CALL READ (*1720,*1710,itape4,ibuf2,20,1,iflg)
  IF (mach == 2) CALL umftrn (ibuf2)
  IF (mach == 3 .AND. kin == 1) CALL umffd (ibuf2)
  CALL page2 (-1)
  WRITE (outtap,200) ibuf2
END DO
IF (ncnt == 0) GO TO 1690
CALL page2 (4)
WRITE  (outtap,1686) ufm,ncnt
1686 FORMAT (a23,' 206, PREVIOUS',i7,' CONTINUATION CARDS, THOUGH ',  &
    'VALID, CANNOT BE PROCESSED', /5X,  &
    'BECAUSE OF ERRORS ON OTHER RELATED CONTINUATION CARDS.',/)
1690 CALL CLOSE (itape4,1)

!     REACTIVE DIAG 47 TO PRINT THE CONTENTS OF NTPT

l47 = 0
IF (l47 == 0) GO TO 1699
CALL OPEN (*1750,nptp,buf(1),0)
1691 CALL skpfil (nptp,+1)
CALL READ (*1697,*1697,nptp,ibuf1(1),2,1,j)
IF (ibuf1(1) /= iblkda(1) .OR. ibuf1(2) /= iblkda(2)) GO TO 1691
1693 CALL READ (*1697,*1697,nptp,ibuf1(1),20,1,j)
WRITE  (outtap,1695) (ibuf1(j),j=1,10),(ibuf1(j),j=17,20)
1695 FORMAT (' ==NPTP==>',5(1X,2A4),'...',2(1X,2A4))
GO TO 1693
1697 CALL CLOSE (nptp,1)
1699 CONTINUE

!     DISABLE FREE-FIELD INPUT OPTION IN XREAD.

ffflag = 0
RETURN

!     ERROR MESSAGES

1700 WRITE  (outtap,1701) sfm
1701 FORMAT (a25,' 210, SCRATCH COULD NOT BE OPENED')
GO TO  1800
1710 WRITE  (outtap,1711) sfm
1711 FORMAT (a25,' 211, ILLEGAL EOR ON SCRATCH')
GO TO  1800
1720 WRITE  (outtap,1721) sfm
1721 FORMAT (a25,' 212, ILLEGAL EOF ON ITAPE4')
GO TO  1800
1730 WRITE  (outtap,1731) sfm
1731 FORMAT (a25,' 213, ILLEGAL EOF ON OPTP')
GO TO  1800
1740 WRITE  (outtap,1741) sfm
1741 FORMAT (a25,' 214, OPTP COULD NOT BE OPENED')
GO TO  1800
1750 WRITE  (outtap,1751) sfm
1751 FORMAT (a25,' 215, NPTP COULD NOT BE OPENED')
GO TO  1800
1760 WRITE  (outtap,1761) sfm
1761 FORMAT (a25,' 216, ILLEGAL INDEX')
GO TO  1800
1770 WRITE  (outtap,1771) sfm
1771 FORMAT (a25,' 219, MISSING ENDDATA CARD.')
1800 CALL page2 (2)
CALL mesage (-37,0,nsort)
RETURN
END SUBROUTINE xsort
