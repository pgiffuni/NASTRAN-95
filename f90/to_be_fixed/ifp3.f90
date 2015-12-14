SUBROUTINE ifp3
     
!     DATA PROCESSING AND GENERATION OF THE AXIS-SYMETRIC-CONICAL SHELL
 
!          CARDS          TYPE       REC.ID-BIT CARDS-FILE,  CARDS-FILE
!     ===  =======      ===========  ========== ===========  ==========
!      1   AXIC     --  AX.SY.SHELL     515- 5
!      2   CCONEAX  --  AX.SY.SHELL    8515-85  CCONE-GEOM2,
!      3   FORCEAX  --  AX.SY.SHELL    2115-21  FORCE-GEOM3,
!      4   FORCE    --  STANDARD       4201-42  FORCE-GEOM3,
!      5   GRAV     --  STANDARD       4401-44   GRAV-GEOM3,
!      6   LOAD     --  STANDARD       4551-61   LOAD-GEOM3,
!      7   MOMAX    --  AX.SY.SHELL    3815-38  MOMNT-GEOM3,
!      8   MOMENT   --  STANDARD       4801-48  MOMNT-GEOM3,
!      9   MPCADD   --  STANDARD       4891-60 MPCADD-GEOM4,
!     10   MPCAX    --  AX.SY.SHELL    4015-40    MPC-GEOM4,
!     11   OMITAX   --  AX.SY.SHELL    4315-43   OMIT-GEOM4,
!     12   POINTAX  --  AX.SY.SHELL    4915-49    MPC-GEOM4, GRID-GEOM1
!     13+  RFORCE   --  STANDARD       5509-55 RFORCE-GEOM3,
!     14   RINGAX   --  AX.SY.SHELL    5615-56    SPC-GEOM4, GRID-GEOM1
!     15   SECTAX   --  AX.SY.SHELL    6315-63    MPC-GEOM4, GRID-GEOM1
!     16   SEQGP    --  STANDARD       5301-53  SEQGP-GEOM1,
!     17   SPCADD   --  STANDARD       5491-59 SPCADD-GEOM4,
!     18   SPCAX    --  AX.SY.SHELL    6215-62    SPC-GEOM4,
!     19   SUPAX    --  AX.SY.SHELL    6415-64 SUPORT-GEOM4,
!     20   TEMPAX   --  AX.SY.SHELL    6815-68   TEMP-GEOM3,
!     21   TEMPD    --  STANDARD       5641-65  TEMPD-GEOM3,
!     22   CTRIAAX  --  AX.TR.CR       7012-70  CTRIA-GEOM2
!     23   CTRAPAX  --  AX.TRA.CR      7042-74  CTRAP-GEOM2
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift    ,rshift    ,andf      ,orf       , complf
 LOGICAL :: secd      ,nogo      ,recoff    ,piez
 REAL :: nphi      ,nphi1     ,nisq      ,ni        ,  &
     a1        ,a2        ,a3        ,a4        ,  &
     angle     ,raddeg    ,pi        ,difphi    ,  &
     rz        ,t1        ,t2        ,coef      , consts    ,sum       ,twopi
 DIMENSION       geom(4)   ,z(8)      ,num(11)   ,inum(11)  ,  &
     msg1(2)   ,msg2(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm       ,uwm       ,uim       ,sfm
 COMMON /system/ ibufsz    ,nout      ,noflag    ,dumdum(8) ,  &
     nlines    ,dum1(26)  ,nbpc      ,nbpw      , dum37(37) ,ipiez
 COMMON /two   / two(32)
 COMMON /condas/ consts(5)
 COMMON /ifp3lv/            recid(3)  ,recid1(3) ,recidx(3) ,  &
     iend      ,REC(3)    ,rec1(3)   ,trail(7)  ,  &
     it        ,axtrl(7)  ,openfl(6) ,n         ,  &
     a1        ,csid      ,ni        ,nisq      ,  &
     a2        ,ibuff1    ,ibuff2    ,ibuff3    ,  &
     a3        ,buff      ,nogo      ,op        ,  &
     a4        ,iheadr    ,ibitr     ,ifile     ,  &
     noreg     ,last      ,ierrtn    ,icont     ,  &
     noaxic    ,ringid    ,outbuf    ,veor      ,  &
     istart    ,iretrn    ,flag      ,iamt      ,  &
     sum       ,ibit      ,setid     ,sorc      ,  &
     ibegin    ,mpcon     ,nwords    ,nnn       ,  &
     angle     ,k3or6     ,nphi1     ,zpt       ,  &
     nmove     ,csset     ,nopont    ,non       ,  &
     iphi      ,recoff    ,nphi      ,n3or5     ,  &
     ion       ,nplus1    ,nosect    ,coef      ,  &
     ipt       ,compon    ,icore     ,iscrat    ,  &
     icore1    ,ncards    ,i1        ,iat       ,  &
     i2        ,t1        ,t2        ,nfile     , nadd      ,ncard
 COMMON /ifp3cm/ FILE(6)   ,iname(12) ,cdtype(50),axic1(3)  ,  &
     cconex(3) ,forcex(3) ,force(3)  ,grav(3)   ,  &
     load(3)   ,momax(3)  ,moment(3) ,mpcadd(3) ,  &
     mpcax(3)  ,omitax(3) ,pointx(3) ,presax(3) ,  &
     ringax(3) ,sectax(3) ,seqgp(3)  ,spcax(3)  ,  &
     supax(3)  ,tempax(3) ,tempd(3)  ,pload(3)  ,  &
     mpc(3)    ,spc(3)    ,grid(3)   ,suport(3) ,  &
     neg111(3) ,t65535(3) ,temp(3)   ,omit(3)   ,  &
     spcadd(3) ,one       ,zero      ,iheadb(96),  &
     ctriaa(3) ,ctrapa(3) ,iconso    ,rforce(3)
 COMMON /output/ dummy(96) ,ihead(96)
 COMMON /zzzzzz/ rz(1)
 EQUIVALENCE     (consts(1),pi     )  ,(consts(2),twopi  )  ,  &
     (consts(4),raddeg )  ,(z(1)     ,rz(1)  )  ,  &
     (geom(1)  ,FILE(1))  ,(scrtch   ,FILE(5))  ,  &
     (axic     ,FILE(6))  ,(num(11)  ,b      )  ,  &
     (noeor    ,inprwd    ,zero)                ,  &
     (eor      ,clorwd    ,outrwd    ,one    )
 DATA    inum  / 1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H /
 DATA    ifiat / 4HFIAT/, ifist  /  4HFIST/,i5,i6 / 5,6     /
DATA    msg1  / 4HIFP3 , 4HBEGN /, msg2  / 4HIFP3, 4HEND   /

CALL conmsg (msg1,2,0)

!     RIGHT-JUSTIFY INUM AND CALL IT NUM

DO  i = 1,11
  num(i) = rshift(inum(i),nbpw-nbpc)
END DO

!     INITIAL CHECK TO MAKE SURE TRAILER BITS ARE ALL OFF FOR GEOM1,
!     GEOM2, GEOM3, GEOM4.

DO  i = 1,96
  ihead(i) = iheadb(i)
END DO

IF (noflag == 0) THEN
  GO TO    20
ELSE
  GO TO    30
END IF
20 nogo = .false.
GO TO 40
30 nogo = .true.

40 openfl(1) = 0
openfl(2) = 0
openfl(3) = 0
openfl(4) = 0
openfl(5) = 0
openfl(6) = 0
loop110:  DO  i = 1,4
  trail(1) = geom(i)
  CALL rdtrl (trail(1))
  IF (trail(1) > 0.0) THEN
    GO TO    70
  END IF
  50 CALL page2 (3)
  imsg = 1061
  WRITE  (nout,55) sfm,imsg
  55 FORMAT (a25,i5)
  WRITE  (nout,60) geom(i),iname(2*i-1),iname(2*i),ifiat
  60 FORMAT (5X,11HFILE NUMBER,i4,3H ( ,2A4,12H) is NOT in ,a4,1H.)
  nogo = .true.
  CYCLE loop110
  
  70 DO  j = 2,7
    IF (trail(j) == 0.0) THEN
      GO TO   100
    END IF
    80 CALL page2 (3)
    imsg = 1062
    WRITE  (nout,55) sfm,imsg
    WRITE  (nout,90) geom(i),iname(2*i-1),iname(2*i)
    90 FORMAT (5X,'FILE NUMBER',i4,3H ( ,2A4,') HAS TRAILER BIT ON.  ',  &
        'FILE SHOULD BE CLEAN AT ENTRY TO IFP3.')
    nogo = .true.
    CYCLE loop110
  END DO
END DO loop110

!     PROCEED TO SETUP CORE AND OPEN AXIC FILE
!     ICORE1 WILL ALWAYS EQUAL THE GROSS OPEN CORE TO IFP3 AT START

icore1 = korsz(z)
ibuff1 = icore1 - ibufsz - 2
ibuff2 = ibuff1 - ibufsz
ibuff3 = ibuff2 - ibufsz
icore  = ibuff3 - 1
icrq   = 100 - icore
IF (icore < 100) GO TO 1310

!     OPEN  AXIC FILE

CALL preloc (*1330,z(ibuff1),axic)
openfl(6) = 1
axtrl(1)  = axic
CALL rdtrl (axtrl(1))

!     READ AXIC CARD

CALL locate (*130,z(ibuff1),axic1(1),flag)
CALL READ (*1600,*130,axic,z(1),2,eor,flag)
n    = z(1)
csid = z(2)
nnn  = n
ncard= 1
ASSIGN 140 TO ierrtn
GO TO 1420

!     MISSING REQUIRED AXIC CARD

130 ASSIGN 140 TO ierrtn
nnn   = 0
ncard = 1
GO TO 1510
140 n = nnn
nplus1 = n + 1


!     GEOM2  PROCESSING
!     =================

!     OPEN GEOM2

ifile= geom(2)
i    = 2
op   = outrwd
buff = ibuff2
ASSIGN 150 TO iretrn
GO TO 1350

!     CCONEAX CARDS

150 REC(1) = cconex(1)
REC(2) = cconex(2)
REC(3) = cconex(3)
ncard  = 2

!     IF THERE IS NO CCONEAX CARD, THEN GO TO 1750 AND LOOK FOR
!     CTRAPAX OR CTRIAAX CARDS

iconb  = 0
iconso = 0
CALL locate (*1750,z(ibuff1),REC(1),flag)

!     INPUT IS IN 4-WORD CARDS
!     OUTPUT IS N+1 4-WORD CARDS FOR EACH CARD INPUT

!     RECORD HEADER FOR CCONES

ASSIGN 160 TO iheadr
GO TO 1470

160 CALL READ (*1600,*220,axic,z(1),4,noeor,iamt)

!     CHECK RING ID-S FOR SIZE

nnn = z(3)
ASSIGN 170 TO ierrtn
GO TO 1440
170 nnn = z(4)
ASSIGN 180 TO ierrtn
GO TO 1440

!     CHECK CCONEAX ID FOR 1-9999 ALLOWABLE RANGE

180 IF (z(1) > 0 .AND. z(1) < 10000) GO TO 200
CALL page2 (3)
imsg = 361
WRITE  (nout,185) ufm,imsg
185 FORMAT (a23,i4)
WRITE  (nout,190) z(1)
190 FORMAT (5X,'CCONEAX ID =',i10,'.  OUT OF 1 TO 9999 PERMISSIBLE ', 'RANGE')
nogo = .true.

200 z(1) = z(1)*1000
DO  i = 1,nplus1
  z(1) = z(1) + 1
  z(3) = z(3) + 1000000
  z(4) = z(4) + 1000000
  IF (nogo) CYCLE
  CALL WRITE (geom(2),z(1),4,noeor)
END DO
GO TO 160

!     OUT OF CCONEAX CARDS

220 IF (iamt == 0) THEN
  GO TO   240
END IF

!     GO TO 356 FOR RECORD ERROR

230 ASSIGN 260 TO ierrtn
GO TO 1490

!     WRITE EOR AND PUT BITS IN TRAILER

240 ASSIGN 250 TO iretrn
iconso = 1
GO TO 1270
250 iconb = 1
GO TO 1750

!     CLOSE GEOM2

260 i = 2
ASSIGN 270 TO iretrn
GO TO 1380

!     GEOM3 PROCESSING
!     ================

!     OPEN GEOM3

270 ifile= geom(3)
i    = 3
op   = outrwd
buff = ibuff2
ASSIGN 280 TO iretrn
GO TO 1350

!     FORCE, FORCEAX, MOMNT, AND MOMNTAX CARDS

280 recid(1) = force(1)
recid(2) = force(2)
recid(3) = force(3)
recidx(1)= forcex(1)
recidx(2)= forcex(2)
recidx(3)= forcex(3)
ncard    = 3
ASSIGN 620 TO icont

!     SET NOREG = 0 OR 1, DEPENDING ON PRESSENCE OF RECID
!     SET NOAXIC= 0 OR 1, DEPENDING ON PRESSENCE OF RECIDX

290 ibit = recidx(2)
ASSIGN 300 TO ibitr
GO TO 1460
300 noaxic = non
ibit = recid(2)
ASSIGN 310 TO ibitr
GO TO 1460
310 noreg = non

REC(1) = recid(1)
REC(2) = recid(2)
REC(3) = recid(3)

IF (noaxic == 0) THEN
  GO TO   320
ELSE
  GO TO   340
END IF
320 IF (noreg  == 0) THEN
  GO TO   610
END IF

!     TRANSFER FORCE OR MOMENT RECORD DIRECTLY.
!     THERE ARE NO FORCEAX OR MOMAX CARDS RESPECTIVELY.

330 ASSIGN 610 TO iretrn
GO TO 1230

!     AT 410 READ IN ALL FORCEAX OR MOMNTAX CARDS AND PUT OUT ON GEOM(3)
!     IF NOREG=0,AND ON SCRTCH IF NOREG NON-ZERO.FIRST WRITE 3-WORD-
!     REC ID ON GEOM3.

340 ASSIGN 350 TO iheadr
GO TO 1470

!     OPEN SCRATCH IF NEEDED

350 IF (noreg == 0) THEN
  GO TO   370
END IF
360 i    = 5
op   = outrwd
buff = ibuff3
ASSIGN 370 TO iretrn
GO TO 1350
370 CALL locate (*1530,z(ibuff1),recidx(1),flag)
380 CALL READ (*1600,*440,axic,z(1),8,noeor,iamt)

!     CHECK RING ID

ASSIGN 390 TO ierrtn
nnn = z(2)
GO TO 1440

!     CHECK HARMONIC NUMBER AND FOR A SEQUENCE OF HARMONICS

390 IF (z(4) == 0) GO TO 396
ii   = 1
nh1  = 0
nh2  = 0
secd = .true.
word = 4
DO  ij = 1,2
  DO  ix = 1,4
    chr = rshift(lshift(z(word),nbpc*IABS(ix-4)),nbpw-nbpc)
    IF (chr == b) CYCLE
    DO  i = 1,10
      k = i-1
      IF (num(i) == chr) GO TO 394
    END DO
    secd = .false.
    ii   = 1
    CYCLE
    394 IF (secd) GO TO 395
    nh1 = nh1 + ii*k
    ii  = ii*10
    CYCLE
    395 nh2 = nh2 + ii*k
    ii  = ii*10
  END DO
  word = word -1
END DO
IF (nh1 <= nh2) GO TO 398
word = nh1
nh1  = nh2
nh2  = word
398 nnn  = nh1
ASSIGN 397 TO ierrtn
GO TO 1420
396 nh1  = z(3)
nh2  = z(3)
397 nnn  = nh2
ASSIGN 400 TO ierrtn
GO TO 1420
400 z(4) = z(5)
z(5) = z(6)
z(6) = z(7)
z(7) = z(8)
nh1  = nh1 + 1
nh2  = nh2 + 1
sum  = z(2)
mus  = z(2)
DO  i = nh1,nh2
  z(2) = mus + i*1000000
  z(3) = 0
  
!     OUTPUT TO GEOM(3) IF NOREG = 0
!     OUTPUT TO SCRTCH  IF NOREG = NON-ZERO
  
  IF (nogo ) GO TO 380
  IF (noreg == 0) THEN
    GO TO   410
  ELSE
    GO TO   420
  END IF
  410 nfile = geom(3)
  GO TO 430
  420 nfile = scrtch
  430 CALL WRITE (nfile,z(1),7,noeor)
END DO
GO TO 380

!     OUT OF CARDS

440 IF (iamt == 0) THEN
  GO TO   460
END IF

!     CHECK FOR RECORD INCONSISTANCY ERROR.

450 REC(1) = recidx(1)
REC(2) = recidx(2)
REC(3) = recidx(3)
ASSIGN 460 TO ierrtn
GO TO 1490

460 IF (noreg == 0) THEN
  GO TO   590
END IF

!     CLOSE THE SCRTCH FILE AND THEN MERGE SCRTCH WITH AXIC
!     ON TO GEOM3

470 i = 5
ASSIGN 480 TO iretrn
GO TO 1380

!     OPEN SCRTCH FILE FOR INPUT AND LOCATE FORCE OR MOMENT CARDS ON
!     AXIC FILE.

480 ASSIGN 490 TO iretrn
op = inprwd
GO TO 1350
490 CALL locate (*1560,z(ibuff1),recid(1),flag)
IF (nogo) GO TO 610

CALL READ (*1600,*600,axic,z(1),7,noeor,iamt)
CALL READ (*1610,*1610,scrtch,z(8),7,noeor,iamt)

500 IF (z(1) <= z(8)) GO TO 510

nfile  = scrtch
outbuf = 8
GO TO 520

510 nfile  = axic
outbuf = 1

520 IF (nogo) GO TO 610
CALL WRITE (geom(3),z(outbuf),7,noeor)
CALL READ (*1620,*540,nfile,z(outbuf),7,noeor,iamt)
GO TO 500

!     OK ALL WORDS PROCESSED FOR FILE-NFILE

540 IF (nfile == axic) GO TO 550
nfile  = axic
outbuf = 1
GO TO 560
550 nfile  = scrtch
outbuf = 8
560 IF (nogo) GO TO 610
CALL WRITE (geom(3),z(outbuf),7,noeor)
CALL READ (*1620,*580,nfile,z(outbuf),7,noeor,iamt)
GO TO 560

!     CLOSE SCRTCH, WRITE EOR, AND PUT BITS IN TRAILER.

580 i = 5
ASSIGN 590 TO iretrn
GO TO 1380
590 ASSIGN 610 TO iretrn
GO TO 1270

!     RECORD LENGTH ERROR

600 REC(1) = recid(1)
REC(2) = recid(2)
REC(3) = recid(3)
ASSIGN 610 TO ierrtn
GO TO 1490

610 GO TO icont, (620,650)

!     GRAV CARD

620 REC(1) = grav(1)
REC(2) = grav(2)
REC(3) = grav(3)
ASSIGN 630 TO iretrn
GO TO 1230

!     LOAD CARD

630 REC(1) = load(1)
REC(2) = load(2)
REC(3) = load(3)
ASSIGN 640 TO iretrn
GO TO 1230

!     MOMENT AND MOMAX CARDS

640 recid(1)  = moment(1)
recid(2)  = moment(2)
recid(3)  = moment(3)
recidx(1) = momax(1)
recidx(2) = momax(2)
recidx(3) = momax(3)
ncard = 7
ASSIGN 650 TO icont
GO TO 290

!     PRESAX CARD

650 CALL locate (*722,z(ibuff1),presax(1),flag)

!     RECORD HEADER FOR PRESAX CARDS IS FORMED HERE

REC(1) = presax(1)
REC(2) = presax(2)
REC(3) = presax(3)
ncard  = 13
ASSIGN 660 TO iheadr
GO TO 1470

660 CALL READ (*1600,*700,axic,z(1),6,noeor,iamt)

!     CREATE N+1 CARDS OF SAME LENGTH AS INPUT CARD.

!     CHECK RING ID-S IN FIELDS 3 AND 4 FOR PROPER SIZE.

!     CHECK FOR PIEZOELECTRIC

piez = .false.
IF (ipiez == 1 .AND. z(3) < 0) piez = .true.
IF (.NOT. piez) GO TO 661
z(3) = -z(3)
661 CONTINUE
nnn = z(3)
ASSIGN 670 TO ierrtn
GO TO 1440
670 nnn = z(4)
ASSIGN 680 TO ierrtn
GO TO 1440

680 difphi = ABS(rz(i6) - rz(i5))
DO  i = 1,nplus1
  z(7) = i - 1
  z(3) = z(3) + 1000000
  IF (piez) z(3) = -z(3)
  z(4) = z(4) + 1000000
  IF (nogo) CYCLE
  IF (difphi == 0.0) CYCLE
  IF (i > 1 .AND. ABS(difphi-360.) < 1.e-6) CYCLE
  CALL WRITE (geom(3),z(1),7,noeor)
  IF (piez) z(3) = -z(3)
END DO
GO TO 660

!     OUT OF PRESAX CARDS

700 IF (iamt == 0) THEN
  GO TO   720
END IF

!     CHECK FOR RECORD INCONSISTANCY ERROR.

710 ASSIGN 722 TO ierrtn
REC(1) = presax(1)
REC(2) = presax(2)
REC(3) = presax(3)
GO TO 1490

!     WRITE EOR AND PUT BITS IN TRAILER

720 ASSIGN 722 TO iretrn
GO TO 1270

!     RFORCE CARD

722  CALL locate (*730,z(ibuff1),rforce(1),flag)
REC(1) = rforce(1)
REC(2) = rforce(2)
REC(3) = rforce(3)
ncard  = 24
ASSIGN 723 TO iheadr
GO TO 1470

!     PROCESS RFORCE DATA

723 CALL READ (*1600,*725,axic,z(1),7,noeor,iamt)
IF (z(2) == 0 .AND. z(3) == 0 .AND. z(5) == 0 .AND. z(6) == 0) GO TO 7240
WRITE  (nout,724) ufm,z(1)
724 FORMAT (a23,' 336, RFORCE DATA IN SET NO.',i8,  &
    ' CONTAINS ILLEGAL DIRECTION FOR AXISYMMETRIC PROBLEM')
nogo = .true.
GO TO 723
7240 z(2) = 0
z(3) = 0
z(5) = 0
z(6) = z(7)
z(7) = 0
CALL WRITE (geom(3),z(1),7,noeor)
GO TO 723

!     END OF RFORCE CARDS

725 IF (iamt == 0) THEN
  GO TO   727
END IF

!     RECORD INCONSISTENCY ERROR

726 ASSIGN 730 TO ierrtn
REC(1) = rforce(1)
REC(2) = rforce(2)
REC(3) = rforce(3)
GO TO 1490

!     WRITE EOR AND BITS IN TRAILER

727 ASSIGN 730 TO iretrn
GO TO 1270

!     TEMPD CARD

730 REC(1) = tempd(1)
REC(2) = tempd(2)
REC(3) = tempd(3)
ASSIGN 740 TO iretrn
IF (nogo) GO TO 740
CALL locate (*740,z(ibuff1),REC(1),flag)
CALL WRITE (ifile,REC(1),3,noeor)
veor = 0
735 CALL READ (*1600,*738,axic,z(1),icore,noeor,iamt)
iamt = icore
736 DO  i = 1,iamt,2
  z(i) = z(i) + 100000000
END DO
CALL WRITE (ifile,z(1),iamt,0)
DO  i = 1,iamt,2
  z(i) = z(i) + 100000000
END DO
CALL WRITE (ifile,z(1),iamt,veor)
IF (veor == 0.0) THEN
  GO TO   735
ELSE
  GO TO  1290
END IF
738 veor = 1
GO TO 736

!     TEMPAX CARD

740 CALL locate (*1210,z(ibuff1),tempax(1),flag)

!     RECORD HEADER ON GEOM3 FOR TEMP CARDS

REC(1) = temp(1)
REC(2) = temp(2)
REC(3) = temp(3)
ncard  = 20
ASSIGN 750 TO iheadr
GO TO 1470

!     AT 604(?) SET UP SCRATCH FILE.

750 i    = 5
buff = ibuff3
op   = outrwd
ASSIGN 760 TO iretrn
GO TO 1350

!     PICK UP FIRST TEMPAX CARD = 4 WORDS.

760 last = 0
CALL READ (*1600,*1200,axic,z(1),4,noeor,iamt)
770 k = 0
setid = z(1)
ringid= z(2)

!     CHECK RING ID FOR PROPER RANGE OF VALUE

nnn = ringid
ASSIGN 780 TO ierrtn
GO TO 1440

780 iat = 3
790 k   = k + 1
iat = iat + 2
icrq= iat + 3 - icore
IF (icore < iat+3) GO TO 1310

!     ALL TEMPAX CARDS HAVING SAME SET AND RING ID MUST BE ABLE TO
!     HAVE 2 WORDS EACH FIT IN CORE.

z(iat  ) = z(3)
z(iat+1) = z(4)

CALL READ (*1600,*1130,axic,z(1),4,noeor,iamt)

!     DOES THIS CARD HAVE SAME SET AND RING ID AS LAST IN CURRENT SERIES

IF (z(1) /= setid ) GO TO 800
IF (z(2) /= ringid) GO TO 800
GO TO 790

!     WE HAVE A  K X 2  ARRAY OF  PHI-S  AND T-S.

!     CONVERT ALL  PHIS SUCH THAT (0.LE. PHI .LT.TWOPI)

800 iend   = iat + 1
ibegin = 5

DO  i = ibegin,iend,2
  angle = rz(i)
  IF (angle < 0.0) THEN
    GO TO   810
  ELSE IF (angle == 0.0) THEN
    GO TO   840
  ELSE
    GO TO   830
  END IF
  810 IF (angle < 0.0) THEN
    GO TO   820
  ELSE
    GO TO   840
  END IF
  820 angle = angle + 360.0
  GO TO 810
  
  830 IF (angle < 360.0) GO TO 840
  angle = angle - 360.0
  GO TO 830
  
  840 rz(i) = angle*raddeg
END DO

!     SIMPLE SORT FOR THE K X 2  MATRIX.
!     SORT IS PERFORMED ON COLUMN 1 ONLY

IF (k == 1) GO TO 950
istart = ibegin + 2
DO  i = istart,iend,2
  iat = i - 2
  IF (rz(i) >= rz(iat)) CYCLE
  
!     ROW NOT HIGH ENOUGH.  MOVE IT UP.
  
  850 iat = iat - 2
  IF (iat-ibegin > 0) THEN
    GO TO   860
  ELSE
    GO TO   870
  END IF
  860 IF (rz(i) < rz(iat)) GO TO 850
  iat = iat + 2
  GO TO 880
  870 iat = ibegin
  
!     THE ELEMENTS (I) AND (I+1) WILL BE MOVED UP TO POSITIONS (IAT) AND
!     (IAT+1) AND ELEMENTS (IAT) THRU (I-1) WILL BE  MOVED DOWN 1 ROW.
  
!     FIRST SAVE THE ROW BEING MOVED UP
  
  880 rz(iend+1) = rz(i)
  rz(iend+2) = rz(i+1)
  nmove = i - iat
  iat   = i + 2
  DO  j = 1,nmove
    iat = iat - 1
    rz(iat) = rz(iat-2)
  END DO
  
!     REPLACE SAVE ROW IN NEW SLOT
  
  rz(iat-2) = rz(iend+1)
  rz(iat-1) = rz(iend+2)
  
END DO

!     CHECK FOR ANY DUPLICATE ANGLES AND REMOVE THEM...

ibegin = ibegin + 2
910 DO  i = ibegin,iend,2
  IF (z(i) == z(i-2)) GO TO 930
END DO
GO TO 950

!     DUPLICATE, SHRINK LIST UP OVER IT.

930 iend = iend - 2
k=k-1
DO  j = i,iend,2
  z(j  ) = z(j+2)
  z(j+1) = z(j+3)
END DO
ibegin = i
IF (ibegin - iend < 0) THEN
  GO TO   910
END IF

!     SET UP K + 1  CARD

950 rz(iend+1) = rz(i5) + twopi
rz(iend+2) = rz(i6)

!     THERE ARE K CARDS NOW WITH SETID, AND RINGID, NOT INCLUDING THE
!     K + 1ST CARD

!     N+1 TEMP CARDS FOR S SET (PUT ON GEOM3)
!     N+1 TEMP CARDS FOR C SET (PUT ON SCRTCH FOR NOW)

!     NOTE FMMS-52  (10/04/67) PAGE -9- FOR FOLLOWING...

csset = 1
setid = setid + 100000000

!     CSSET = 0 FOR C-SET  AND NON-ZERO FOR S-SET.

ibegin = k + k + 7
icrq   = ibegin + 2 - icore
IF ((ibegin+2) > icore) GO TO 1310

960 nadd = 0
z(ibegin) = setid
DO  i = 1,nplus1
  nadd = nadd + 1000000
  
!     NI IS REAL
  
  ni   = i - 1
  nisq = (i-1)**2
  z(ibegin+1) = ringid + nadd
  iphi = 3
  it   = 4
  sum  = 0.0E0
  IF (ni    == 0) THEN
    GO TO   970
  ELSE
    GO TO  1010
  END IF
  970 IF (csset == 0.0) THEN
    GO TO   980
  ELSE
    GO TO  1000
  END IF
  980 DO  ik = 1,k
    iphi = iphi + 2
    it   = it   + 2
    sum  = sum  + (rz(it)+rz(it+2))*(rz(iphi+2)-rz(iphi))
  END DO
  1000 rz(ibegin+2) = 0.25*sum/pi
  GO TO 1060
  
!     NON-ZERO NI
  
  1010 IF (k == 1) GO TO 1050
  DO  ik = 1,k
    iphi = iphi + 2
    it   = it   + 2
    nphi = ni*rz(iphi  )
    nphi1= ni*rz(iphi+2)
    
    IF (csset == 0.0) THEN
      GO TO  1020
    ELSE
      GO TO  1030
    END IF
    
!     C-SET
    
    1020 a1 =  SIN(nphi1)
    a2 = -SIN(nphi )
    a3 =  COS(nphi1)
    a4 = -COS(nphi )
    GO TO 1040
    
!     S-SET
    
    1030 a1 = -COS(nphi1)
    a2 =  COS(nphi )
    a3 =  SIN(nphi1)
    a4 = -SIN(nphi )
    
    
    1040 sum = sum + (((rz(it)*rz(iphi+2) - rz(it+2)*rz(iphi))*  &
        (a1 + a2)/ni) + ((rz(it+2) - rz(it))*  &
        (a3 + a4 + nphi1*a1 + nphi*a2)/nisq))/ (rz(iphi+2) - rz(iphi))
  END DO
  
  1050 rz(ibegin+2) = sum/pi
  
  1060 IF (nogo ) EXIT
  IF (csset == 0.0) THEN
    GO TO  1080
  END IF
  1070 nfile = geom(3)
  GO TO 1090
  1080 nfile = scrtch
  1090 CALL WRITE (nfile,z(ibegin),3,noeor)
END DO
1105 IF (csset == 0.0) THEN
  GO TO  1120
END IF
1110 csset = 0
setid = setid + 100000000
GO TO 960

!     THIS SERIES OF TEMPAX CARDS COMPLETE GO FOR MORE IF LAST = 0

1120 IF (last == 0) THEN
  GO TO   770
ELSE
  GO TO  1140
END IF
1130 last = 1
GO TO 800

!     ALL TEMPAX CARDS COMPLETE. CLOSE SCRATCH, OPEN SCRATCH
!     AND COPY SCRATCH TO GEOM3.

1140 IF (nogo) GO TO 1210
CALL WRITE (scrtch,z(1),0,eor)
CALL CLOSE (scrtch,clorwd)
CALL OPEN  (*1640,scrtch,z(ibuff3),inprwd)

veor = 0
1150 CALL READ (*1610,*1170,scrtch,z(1),icore,noeor,iamt)
iamt = icore
1160 CALL WRITE (geom(3),z(1),iamt,veor)
IF (veor == 0.0) THEN
  GO TO  1150
ELSE
  GO TO  1180
END IF
1170 veor = 1
GO TO 1160

!     ALL  TEMPAX  CARDS  PROCESSED.

1180 CALL CLOSE (scrtch,clorwd)

!     PUT BITS IN TRAILER FOR TEMP CARDS WRITTEN

REC(1) = temp(1)
REC(2) = temp(2)
REC(3) = temp(3)
ASSIGN 1210 TO iretrn
GO TO 1290

!     RECORD LENGTH ERROR

1200 REC(1) = tempax(1)
REC(2) = tempax(2)
REC(3) = tempax(3)
ASSIGN 1210 TO ierrtn
GO TO 1490

!     CLOSE GEOM3

1210 i = 3
ASSIGN 1220 TO iretrn
GO TO 1380

!     CTRIAAX CARD

1700 REC(1) = ctriaa (1)
REC(2) = ctriaa (2)
REC(3) = ctriaa (3)
ncard  = 43
CALL locate (*1800,z(ibuff1),REC(1),flag)

!     RECORD HEADER FOR CTRIAAX

ASSIGN 1710 TO iheadr
iconb  = 2
iconso = 1
GO TO 1470
1710 CALL READ (*1600,*1770,axic,z(1),6,noeor,iamt)
z(1) = z(1)*1000
DO  i = 1,nplus1
  z(1) = z(1) + 1
  z(3) = z(3) + 1000000
  z(4) = z(4) + 1000000
  z(5) = z(5) + 1000000
  IF (nogo) CYCLE
  CALL WRITE (geom(2),z(1),6,noeor)
END DO
GO TO 1710

!     OUT OF CTRIAAX CARD

1770 IF (iamt == 0) THEN
  GO TO  1740
END IF
1730 ASSIGN 260 TO ierrtn
GO TO 1490

!     PUT BITS IN TRILER

1740 ASSIGN 260 TO iretrn
GO TO 1270
1800 IF (iconso == 1) GO TO 1740
ASSIGN 260 TO ierrtn

!     MISSING REQUIRED CCONEAX OR CTRIAAX OR CTRAPAX CARD

CALL page2 (3)
imsg = 362
WRITE (nout,185) ufm,imsg
WRITE (nout,1910) cdtype(3),cdtype(4),cdtype(43),cdtype(44),  &
    cdtype(45),cdtype(46)
1910 FORMAT (5X,'MINIMUM PROBLEM REQUIRES ',2A4,2H,  ,2A4,4H OR ,2A4,  &
    ' CARD.  NONE FOUND')
nogo = .true.
GO TO ierrtn, (260,240)

!     CTRAPAX CARD
!     ============

1750 REC(1) = ctrapa (1)
REC(2) = ctrapa (2)
REC(3) = ctrapa (3)
CALL locate (*1700,z(ibuff1),REC(1),flag)
iconb  = 1

!     RECORD HEADER FOR CTRAPAX

ASSIGN 1751 TO iheadr
iconso = 1
GO TO 1470
1751 CALL READ (*1600,*1753,axic,z(1),7,noeor,iamt)
z(1) = z(1)*1000
DO  i = 1,nplus1
  z(1) = z(1) + 1
  z(3) = z(3) + 1000000
  z(4) = z(4) + 1000000
  z(5) = z(5) + 1000000
  z(6) = z(6) + 1000000
  IF (nogo) CYCLE
  CALL WRITE (geom(2),z(1),7,noeor)
END DO
GO TO 1751

!     OUT OF CTRAPAX CARD

1753 IF (iamt == 0) THEN
  GO TO  1755
END IF
1754 ASSIGN 260 TO ierrtn
GO TO 1490

!     PUT BITS IN TRILER

1755 ASSIGN 260 TO iretrn
IF (nogo) GO TO 1300
CALL WRITE (ifile,z(1),iamt,eor)
i1 = (REC(2)-1)/16 + 2
i2 = REC(2) - (i1-2)*16 + 16
trail (i1) = orf(trail(i1),two(i2))
GO TO 1700

!     GEOM4 AND GEOM1 PROCESSING IS PERFORMED IN IFP3B ROUTINE
!                                                =====

1220 CALL ifp3b
GO TO 1570

!     UTILITY SECTION FOR IFP3
!     AXIS-SYMETRIC-CONICAL-SHELL DATA GENERATOR.
!     ==========================================

!     COMMON CODE FOR TRANSFER OF RECORD FROM AXIC FILE TO SOME
!     OTHER FILE

1230 CALL locate (*1300,z(ibuff1),REC(1),flag)
IF (nogo) GO TO 1300
CALL WRITE (ifile,REC(1),3,noeor)
1260 CALL READ (*1600,*1270,axic,z(1),icore,noeor,iamt)
iamt = icore
CALL WRITE (ifile,z(1),iamt,noeor)
GO TO 1260
1270 IF (nogo) GO TO 1300
IF (ifile == geom(3)) GO TO 1280
IF (ifile == geom(2) .AND. iconb == 1) GO TO 1300
1280 CALL WRITE (ifile,z(1),iamt,eor)

!     PUT BITS IN TRAILER

1290 i1 = (REC(2)-1)/16 + 2
i2 =  REC(2) - (i1-2)*16 + 16
trail(i1) = orf(trail(i1),two(i2))

1300 GO TO iretrn, (250,260,610,630,640,722,730,740,1210)

!     OUT OF CORE

1310 CALL page2 (4)
imsg = 363
WRITE  (nout, 185) imsg
WRITE  (nout,1320) icrq
1320 FORMAT (5X,'INSUFFICIENT CORE TO PROCESS AXIC DATA IN SUBROUTINE',  &
    'IFP3', /5X,'ADDITIONAL CORE NEEDED =',i8,' WORDS.')
nogo = .true.

!     GO TO FATAL ERROR RETURN

GO TO 1570

!     AXIC FILE NOT IN FIST

1330 CALL page2 (3)
imsg = 1061
WRITE (nout,55) sfm,imsg
WRITE (nout,60) axic,iname(11),iname(12),ifist
nogo = .true.

!     GO TO FATAL ERROR RETURN

GO TO 1570

!     OPEN A FILE AND GET THE TRAILER

1350 IF (nogo) GO TO 1360
CALL OPEN (*1370,FILE(i),z(buff),op)
openfl(i) = 1
IF (i > 4) GO TO 1360

!     WRITE THE HEADER RECORD

CALL WRITE (FILE(i),iname(2*i-1),2,eor)
trail(1) = FILE(i)
CALL rdtrl (trail(1))

1360 GO TO iretrn, (150,280,370,760,490)

1370 CALL page2 (3)
imsg = 1061
WRITE (nout,55) sfm,imsg
WRITE (nout,60) FILE(i),iname(2*i-1),iname(2*i),ifist
nogo = .true.
GO TO 1570

!     CLOSE A FILE

1380 IF (openfl(i) == 0.0) THEN
  GO TO  1410
END IF
1390 IF (i > 4) GO TO 1400
CALL WRITE (FILE(i),t65535(1),3,eor)
1400 CALL CLOSE (FILE(i),clorwd)
openfl(i) = 0
IF (i > 4) GO TO 1410
CALL wrttrl (trail(1))
1410 GO TO iretrn, (270,590,1220,480)

!     HARMONIC NUMBER ... ON CARD TYPE ...... IS OUT OF RANGE 0 TO 998

1420 IF (nnn < 999 .AND. nnn >= 0 .AND. nnn <= n) GO TO ierrtn, (140,400,397)
CALL page2 (3)
imsg = 364
WRITE  (nout,185 ) ufm,imsg
WRITE  (nout,1430) nnn,cdtype(2*ncard-1),cdtype(2*ncard),n
1430 FORMAT (5X,'HARMONIC NUMBER ',i6,4H on ,2A4,' CARD OUT OF 0 TO ',  &
    i4,' ALLOWABLE RANGE.')
nogo = .true.
GO TO ierrtn, (140,400,397)

!     RING ID OUT OF PERMISSABLE RANGE OF 1 TO 999999

1440 IF (nnn > 0 .AND. nnn <= 999999) GO TO ierrtn, (170,180,390,670,680,780)
CALL page2 (3)
imsg = 365
WRITE  (nout,185 ) ufm,imsg
WRITE  (nout,1450) nnn,cdtype(2*ncard-1),cdtype(2*ncard)
1450 FORMAT (5X,'RING ID',i10,4H on ,2A4,' CARD OUT OF 1 TO 999999 ',  &
    'ALLOWABLE RANGE')
nogo = .true.
GO TO ierrtn, (170,180,390,670,680,780)

!     CHECK BIT-IBIT IN TRAILER AND RETURN NON = ZERO OR NON-ZERO...

1460 i1 = (ibit-1)/16 + 2
i2 = ibit - (i1-2)*16 + 16
non = andf(axtrl(i1),two(i2))
GO TO ibitr, (300,310)

!     WRITE 3 WORD RECORD HEADER

1470 IF (nogo) GO TO 1480
CALL WRITE (ifile,REC(1),3,noeor)
1480 GO TO iheadr, (160,350,660,723,750,1710,1751)

!     END-OF-RECORD ON AXIC FILE

1490 CALL page2 (3)
imsg = 1063
WRITE  (nout,55) sfm,imsg
WRITE  (nout,1500) cdtype(2*ncard-1),cdtype(2*ncard)
1500 FORMAT (5X,'EOR ON AXIC FILE WHILE READING ',2A4,'CARD RECORDS.')
nogo = .true.
GO TO ierrtn, (260,460,610,722,730,1210)

!     MISSING REQUIRED CARD

1510 CALL page2 (3)
imsg = 362
WRITE  (nout,185 ) ufm,imsg
WRITE  (nout,1520) cdtype(2*ncard-1),cdtype(2*ncard)
1520 FORMAT (5X,'MINIMUM PROBLEM REQUIRES ',2A4,' CARD.  NONE FOUND.')
nogo = .true.
GO TO ierrtn, (260,140)

!     AXIC TRAILER BIT ON BUT CAN NOT LOCATE RECORD

1530 CALL page2 (3)
imsg = 1064
WRITE  (nout,55) sfm,imsg
WRITE  (nout,1540) cdtype(2*ncard-1),cdtype(2*ncard)
1540 FORMAT (5X,2A4,' CARD COULD NOT BE LOCATED ON AXIC FILE AS ', 'EXPECTED')
1550 nogo = .true.
GO TO 610
1560 CALL page2 (2)
WRITE (nout,1540) recid(1),recid(2),recid(3)
GO TO 1550

!     CLOSE ANY OPEN FILES AND RETURN

1570 DO  i = 1,6
  IF (openfl(i) == 0.0) THEN
    GO TO  1590
  END IF
  1580 CALL CLOSE (FILE(i),clorwd)
  openfl(i) = 0
END DO
IF (nogo) noflag = 32767
CALL conmsg (msg2,2,0)
RETURN

!     EOF ENCOUNTERED READING AXIC FILE.

1600 nfile = axic
in  = 11
in1 = 12
GO TO 1620
1610 nfile = scrtch
in  = 9
in1 = 10
1620 CALL page2 (3)
imsg = 3002
WRITE  (nout,55) sfm,imsg
WRITE  (nout,1630) iname(in),iname(in1),nfile
1630 FORMAT (5X,'EOF ENCOUNTERED WHILE READING DATA SET ',2A4,' (FILE',  &
    i4,') IN SUBROUTINE IFP3')
nogo = .true.
GO TO 1570

1640 i = 5
GO TO 1370
END SUBROUTINE ifp3
