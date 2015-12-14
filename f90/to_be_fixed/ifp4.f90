SUBROUTINE ifp4
     
!     HYDROELASTIC PREFACE ROUTINE
 
!     THIS PREFACE MODULE OPERATES ON FLUID RELATED INPUT DATA WHICH
!     EXISTS AT THIS POINT IN THE FORM OF CARD IMAGES ON THE AXIC DATA
!     BLOCK.
 
!     7/12/73 NO AXIAL SYMMETRY FIRST FIVE WORDS OF BNDFL NO WRITTEN
 
!     THE FOLLOWING LIST GIVES THE CARD IMAGES IFP4 WILL LOOK FOR ON THE
!     AXIC DATA BLOCK, THE CARD IMAGES IFP4 WILL GENERATE OR MODIFY, AND
!     THE DATA BLOCKS ONTO WHICH THE GENERATED OR MODIFIED CARD IMAGES
!     WILL BE PLACED.
 
!     IFP4 INPUT         IFP4 OUTPUT        DATA BLOCK
!     CARD IMAGE         CARD IMAGE         EFFECTED
!     -----------        -----------        ----------
!       AXIF               -NONE-             -NONE-
!       BDYLIST            -DATA-             MATPOOL
!       CFLUID2            CFLUID2            GEOM2
!       CFLUID3            CFLUID3            GEOM2
!       CFLUID4            CFLUID4            GEOM2
!       FLSYM              -DATA-             MATPOOL
!       FREEPT             SPOINT             GEOM2
!                          MPC                GEOM4
!       FSLIST             CFSMASS            GEOM2
!                          SPC                GEOM4
!       GRIDB              GRID               GEOM1
!       PRESPT             SPOINT             GEOM2
!                          MPC                GEOM4
!       RINGFL             GRID               GEOM1
!                          SEQGP              GEOM1
!       DMIAX              DMIG               MATPOOL
 
!     SOME OF THE ABOVE OUTPUT CARD IMAGES ARE A FUNCTION OF SEVERAL
!     INPUT CARD IMAGES
 
LOGICAL :: harms    ,anygb    ,END      ,any      ,g1eof    ,  &
    g2eof    ,g4eof    ,set102   ,press    ,bit      ,  &
    nogo     ,mateof   ,anygrd   ,bit2
INTEGER :: axif(2)  ,bdylst(2),cfluid(6),flsym(2) ,freept(2),  &
    fslst(2) ,gridb(2) ,prespt(2),ringfl(2),cfsmas(2),  &
    mpc(2)   ,mpcadd(2),TYPE(2)  ,spoint(2),cord(8)  ,  &
    spc(2)   ,spcadd(2),spc1(2)  ,ncord(4) ,geom1    ,  &
    subr(2)  ,buf(10)  ,last(10) ,axic     ,geom2    ,  &
    card(10) ,FILE     ,matpol   ,geom4    ,seqgp(2) ,  &
    scrt1    ,entrys   ,corsys   ,SPACE    ,core     ,  &
    scrt2    ,sysbuf   ,output   ,flag     ,eor      ,  &
    csf      ,rd       ,rdrew    ,wrt      ,wrtrew   ,  &
    cls      ,clsrew   ,saveid(5),mones(4) ,dmig(2)  ,  &
    buf1     ,buf2     ,buf3     ,buf4     ,buf5     ,  &
    words    ,bndfl(2) ,trail(7) ,z        ,point    ,  &
    dmiax(2) ,msg1(2)  ,msg2(2)  ,grid(2)
REAL :: rbuf(10) ,rcard(10),rz(4)
CHARACTER (LEN=23) :: ufm
COMMON /xmssg / ufm
COMMON /system/ sysbuf   ,output   ,nogo     ,dum34(34),iaxif
COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,clsrew   , cls
COMMON /zzzzzz/ z(1)
EQUIVALENCE     (z(1),rz(1)), (buf(1),rbuf(1)), (card(1),rcard(1))  &
    ,               (core,icore), (rhob,irhob), (bd,ibd)
DATA    axif  / 8815  ,88    /
DATA    bdylst/ 8915  ,89    /
DATA    cfluid/ 7815  ,78    , 7915  ,79   , 8015  ,80 /
DATA    flsym / 9115  ,91    /
DATA    freept/ 9015  ,90    /
DATA    fslst / 8215  ,82    /
DATA    gridb / 8115  ,81    /
DATA    prespt/ 8415  ,84    /
DATA    ringfl/ 8315  ,83    /
DATA    cfsmas/ 2508  ,25    /
DATA    bndfl / 9614  ,96    /
DATA    mpc   / 4901  ,49    /
DATA    spc   / 5501  ,55    /
DATA    spc1  / 5481  ,58    /
DATA    mpcadd/ 4891  ,60    /
DATA    spcadd/ 5491  ,59    /
DATA    spoint/ 5551  ,49    /
DATA    grid  / 4501  ,45    /
DATA    seqgp / 5301  ,53    /
DATA    dmiax / 214   ,2     /
DATA    dmig  / 114   ,1     /
DATA    cord   /1701  ,17    , 1901  ,19   , 2001 ,20  ,2201 ,22 /
DATA    ncord / 6     ,6     , 13    ,13   /
DATA    mones / -1    ,-1    , -1    ,-1   /
DATA    subr  / 4HIFP4,4H    /
DATA    degrad/ 1.7453292519943E-02/, minus1/ -1 /
DATA    geom1 , geom2,geom4  / 201,208,210  /
DATA    axic  , matpol,eor   / 215,214,1    /

!     NOTE  SCRATCH2 IN IFP4 IS EQUIVALENCED TO THE -FORCE- DATA BLOCK

DATA    scrt1 , scrt2 , noeor /  301,213,0  /
DATA    msg1  / 4HIFP4, 4HBEGN/, msg2       / 4HIFP4, 4HEND  /

!     DEFINE CORE AND BUFFER POINTERS

CALL conmsg (msg1,2,0)
icore = korsz(z)
buf1  = icore- sysbuf - 2
buf2  = buf1 - sysbuf - 2
buf3  = buf2 - sysbuf - 2
buf4  = buf3 - sysbuf - 2
buf5  = buf4 - sysbuf - 2
icore = buf5 - 1
icrq  = 100  - icore
IF (icore < 100) GO TO 2370

!     OPEN AXIC DATA BLOCK (IF NAME NOT IN FIST RETURN - NO MESSAGE)

CALL preloc (*2300,z(buf1),axic)

!     PICK UP AXIF CARD. (IF AXIF CARD NOT PRESENT - RETURN NO MESSAGE)

CALL locate (*2300,z(buf1),axif,flag)
CALL READ (*2320,*30,axic,z(1),icore,eor,words)
WRITE  (output,10) ufm
10 FORMAT (a23,' 4031, INSUFFICIENT CORE TO READ DATA ON AXIF CARD.')
WRITE  (output,20) icore
20 FORMAT (5X,'ADDITIONAL CORE NEEDED =',i8,' WORDS.')
GO TO 2310

!     DATA OF AXIF CARD IS NOW STORED

30 csf   = z(1)
g     = rz(2)
drho  = rz(3)
j     = 3
idrho = z(j)
bd    = rz(4)
nosym = z(j+2)
in    = 6
nn    = words - 1
ni    = nn
j     = nn - in + 1
harms = .false.
IF (j >= 1) harms =.true.
IF (.NOT.harms) GO TO 100

!     CONVERT USER INPUT LIST OF HARMONIC NUMBERS TO A LIST OF INDICES.

IF (j == 1) GO TO 40
CALL sort (0,0,1,1,z(in),j)
40 ii = nn + 1
ni = nn
DO  i = in,nn
  itemp = 2*z(i)
  IF (nosym == 0) THEN
    GO TO    60
  END IF
  50 IF (z(i) == 0) GO TO 60
  ni = ni + 1
  z(ni) = itemp + 1
  60 ni = ni + 1
  z(ni) = itemp + 2
END DO
n  = ni - ii + 1

!     SET MAXIMUM HARMONIC+1 FOR USE BY SDR2C AND VDRB

iaxif = z(nn) + 1

!     BEGIN GEOM1 PROCESSING
!     **********************

!     OPEN GEOM1 AND FIND CORD1C, CORD1S, CORD2C, OR CORD2S CARD IMAGE
!     WITH COORDINATE SYSTEM ID = CSF OF AXIF CARD. THEN NOTE TYPE
!     (CYLINDRICAL OR SPHERICAL, 2 OR 3 RESPECTIVELY)

100 FILE = geom1

!     BEFORE CALLING PRELOC ON GEOM1 CHECK FOR DATA ON GEOM1

trail(1) = geom1
CALL rdtrl (trail)
DO  i = 2,7
  IF (trail(i) == 0.0) THEN
    GO TO   110
  ELSE
    GO TO   120
  END IF
END DO
GO TO 150
120 CALL preloc (*2360,z(buf2),geom1)
DO  i = 1,4
  i2 = 2*i
  CALL locate (*140,z(buf2),cord(i2-1),flag)
  nsize = ncord(i)
  130 CALL READ (*2340,*140,geom1,z(ni+1),nsize,noeor,flag)
  IF (z(ni+1) == csf) GO TO 170
  GO TO 130
END DO

!     FALL THROUGH LOOP IMPLIES COORDINATE SYSTEM WAS NOT FOUND

150 nogo = .true.
WRITE  (output,160) ufm,csf
160 FORMAT (a23,' 4033, COORDINATE SYSTEM ID =',i20,' AS SPECIFIED ',  &
'ON AXIF CARD IS NOT PRESENT', /5X,' AMONG ANY OF CORD1C,',  &
    ' CORD1S, CORD2C, OR CORD2S CARD TYPES.', /5X,  &
    ' CYLINDRICAL TYPE ASSUMED FOR CONTINUING DATA CHECK.')
corsys = 2
GO TO 180
170 corsys = z(ni+2)
180 CALL CLOSE (geom1,clsrew)

!     READ INTO CORE FROM AXIC ALL GRIDB CARD IMAGES (5 WORDS / IMAGE)

anygb  = .false.
igridb = ni + 1
ngridb = ni
CALL locate (*210,z(buf1),gridb,flag)
anygb = .true.
SPACE = core- ni
CALL READ (*2320,*200,axic,z(igridb),SPACE,eor,nwords)
nogo = .true.
WRITE  (output,190) ufm
190 FORMAT (a23,' 4034, INSUFFICIENT CORE TO HOLD GRIDB CARD IMAGES.')
WRITE (output,20) SPACE
anygb = .false.
GO TO 210
200 ngridb = ni + nwords

!     IF ANY GRIDB IMAGES ARE PRESENT A BOUNDARY LIST IS FORMED IN CORE.

210 ibdyl = ngridb + 1
nbdyl = ngridb
IF (.NOT.anygb) GO TO 520
CALL locate (*520,z(buf1),bdylst,flag)
220 CALL READ (*2320,*330,axic,rhob,1,noeor,flag)
IF (irhob /= 1) GO TO 250
IF (idrho /= 1) GO TO 240
nogo = .true.
WRITE  (output,230) ufm
230 FORMAT (a23,' 4035, THE FLUID DENSITY HAS NOT BEEN SPECIFIED ON ',  &
    'A BDYLIST CARD AND', /5X,'THERE IS NO DEFAULT FLUID ',  &
'DENSITY SPECIFIED ON THE AXIF CARD.')
rhob = 1.0
GO TO 250
240 rhob = drho
250 END  = .false.
idfpre = 0
260 CALL READ (*2320,*2330,axic,idf,1,noeor,flag)
IF (idf /= 0) GO TO 270
idfpre = -1
GO TO 260
270 CALL READ (*2320,*2330,axic,idfaft,1,noeor,flag)

!     NOTE.......  ON INPUT   ID=0 IMPLIES AXIS
!                             ID=-1 IMPLIES END OF CARD


!     NOTE.......  ON OUTPUT  ID=0 IMPLIES UNDEFINED ID
!                             ID=-1 IMPLIES AXIS

IF (idfaft == -1) GO TO 280
IF (idfaft ==  0) idfaft = -1
GO TO 290
280 idfaft = 0
END = .true.

!     DO NOT PUT OUT ENTRY WHEN IDF = AXIS

290 IF (idf == -1) GO TO 320
IF (nbdyl+7 <= core) GO TO 310
WRITE  (output,300) ufm
300 FORMAT (a23,' 4036, INSUFFICIENT CORE TO BUILD BOUNDARY LIST ', 'TABLE.')
icrq = nbdyl + 7 - core
GO TO 2370
310 z(nbdyl +1) = idf
z(nbdyl +2) = 1
z(nbdyl +3) = 1
z(nbdyl +4) = 1
z(nbdyl +5) = idfpre
z(nbdyl +6) = idfaft
rz(nbdyl+7) = rhob
nbdyl = nbdyl + 7

!     ROTATE THE ID-S

320 idfpre = idf
idf = idfaft
IF (.NOT.END) GO TO 270
GO TO 220

!     SORT ENTRIES ON FIRST WORD OF EACH ENTRY.

330 CALL sort (0,0,7,1,z(ibdyl),nbdyl-ibdyl+1)
entrys = (nbdyl-ibdyl+1)/7

!     PASS THE RINGFL IMAGES INSERTING X1, X2, AND X3 IN THE APPROPRIATE
!     BDYLIST ENTRY.

CALL locate (*490,z(buf1),ringfl,flag)
340 CALL READ (*2320,*490,axic,buf,4,noeor,flag)
IF (corsys  /=  3) GO TO 360
IF (rbuf(3) /= 0.) GO TO 360
nogo = .true.
WRITE  (output,350) ufm,buf(1)
350 FORMAT (a23,' 5003, ZERO X2 VALUE ON RINGFL CARD WITH SPHERICAL ',  &
    'COORDINATES.  FLUID POINT ID =',i10)
360 IF (buf(corsys+1) == 0.0) THEN
  GO TO   410
END IF
370 nogo = .true.
IF (corsys == 3) GO TO 390
WRITE  (output,380) ufm,buf(1)
380 FORMAT (a23,'4042, COORDINATE SYSTEM IS CYLINDRICAL BUT RINGFL ',  &
    'CARD ID =',i20,' HAS A NON-ZERO X2 VALUE.')
GO TO 410
390 WRITE  (output,400) ufm,buf(1)
400 FORMAT (a23,' 4043, COORDINATE SYSTEM IS SPHERICAL BUT RINGFL ',  &
    'CARD ID =',i20,' HAS A NON-ZERO X3 VALUE.')
410 CALL bisloc(*340,buf(1),z(ibdyl),7,entrys,jpoint)
ntemp = ibdyl + jpoint - 1
IF (z(ntemp+1) == 1) GO TO 430
nogo = .true.
WRITE  (output,420) ufm,buf(1)
420 FORMAT (a23,' 4038, RINGFL CARD HAS ID =',i20,' WHICH HAS BEEN ', 'USED.')
GO TO 340

!     CHECK TO GET RANGE OF BDYLIST HAVING THIS SAME ID.
!     THEN FILL IN X1, X2, AND X3 IN THOSE ENTRIES.

430 nlist = ntemp
440 ntemp = ntemp - 7
IF (ntemp < ibdyl) GO TO 450
IF (z(ntemp) == z(ntemp+7)) GO TO 440
450 ilist = ntemp + 7
ntemp = nlist
460 ntemp = ntemp + 7
IF (ntemp > nbdyl) GO TO 470
IF (z(ntemp) == z(ntemp-7)) GO TO 460
470 nlist = ntemp - 1
DO  i = ilist,nlist,7
  z(i+1) = buf(2)
  z(i+2) = buf(3)
  z(i+3) = buf(4)
END DO
GO TO 340

!     CHECK TO SEE THAT X1, X2, AND X3 WERE FOUND FOR ALL ENTRIES.

490 DO  i = ibdyl,nbdyl,7
  IF (z(i+1) /= 1) CYCLE
  nogo = .true.
  WRITE  (output,500) ufm,z(i)
  500 FORMAT (a23,' 4040, ID =',i20,' APPEARS ON A BDYLIST CARD, BUT ',  &
      'NO RINGFL CARD IS PRESENT WITH THE SAME ID.')
END DO

!     OPEN GEOM1, OPEN SCRATCH1, COPY HEADER REC FROM GEOM1 TO SCRATCH1

520 CALL ifp4c (geom1,scrt1,z(buf2),z(buf3),g1eof)

!     COPY ALL DATA UP TO FIRST GRID CARD IMAGE.

CALL ifp4b (geom1,scrt1,any,z(nbdyl+1),core-nbdyl,grid,g1eof)
anygrd = any
IF (.NOT.anygb) GO TO 1040
IF (nbdyl < ibdyl) GO TO 1040

!     CREATE AND MERGE WITH GRIDS FROM GEOM1, GRIDS FROM GRIDB IMAGES.

FILE = geom1
IF (.NOT.any) GO TO 540
CALL READ (*2340,*530,geom1,last,8,noeor,flag)
CALL ifp4e (last(1))
GO TO 540
530 any = .false.
540 DO  i = igridb,ngridb,5
  card(1) = z(i)
  CALL ifp4e (card(1))
  card(2) = csf
  kid = z(i+4)
  CALL bisloc (*560,kid,z(ibdyl),7,entrys,point)
  ntemp   = ibdyl + point - 1
  card(3) = z(ntemp+1)
  card(4) = z(ntemp+2)
  card(5) = z(ntemp+3)
  card(corsys+2) = z(i+1)
  card(6) = z(i+2)
  card(7) = z(i+3)
  card(8) = 0
  
!     MERGE CARD IN
  
  IF (.NOT.any) GO TO 590
  550 IF (last(1) > card(1)) GO TO 590
  CALL WRITE (scrt1, last, 8, noeor)
  CALL READ (*2340,*580,geom1,last,8,noeor,flag)
  CALL ifp4e (last(1))
  GO TO 550
  560 nogo = .true.
  WRITE  (output,570) ufm,z(i),z(i+4)
  570 FORMAT (a23,' 4057, GRIDB CARD WITH ID =',i10,' HAS A REFERENCE ',  &
      'IDF =',i10,/5X,'WHICH DOES NOT APPEAR IN A BOUNDARY LIST')
  CYCLE
  580 any = .false.
  590 CALL WRITE (scrt1,card,8,noeor)
END DO

IF (.NOT.any) GO TO 620
610 CALL WRITE (scrt1,last,8,noeor)
CALL READ (*2340,*620,geom1,last,8,noeor,flag)
CALL ifp4e (last(1))
GO TO 610

!     FURTHER ALTERATIONS TO BOUNDARY LIST TABLE AT THIS TIME.
!     RADIAL LOCATION (RJ) AND VERTICAL LOCATION (ZJ)

620 nring = ngridb
IF (.NOT.harms) GO TO 1200
DO  i = ibdyl,nbdyl,7
  IF (corsys == 3) GO TO 630
  z(i+2) = z(i+3)
  CYCLE
  
  630 angle = rz(i+2)*degrad
  temp  = rz(i+1)
  rz(i+1) = temp*SIN(angle)
  rz(i+2) = temp*COS(angle)
END DO

!     LENGTH AND ASSOCIATED ANGLE COMPONENTS OF A CONICAL SECTION. L,C,S

IF (nogo) GO TO 780
DO  i = ibdyl,nbdyl,7
  rj = rz(i+1)
  zj = rz(i+2)
  
!     FIND R   , Z     AND  R   , Z     (RJL1,ZJL1,RJP1,ZJP1)
!           J-1   J-1        J+1   J+1
  
  IF (z(i+4) < 0.0) THEN
    GO TO   650
  ELSE IF (z(i+4) == 0.0) THEN
    GO TO   660
  ELSE
    GO TO   670
  END IF
  
!     SECONDARY ID IS AXIS
  
  650 rjl1 = 0
  zjl1 = zj
  GO TO 680
  
!     SECONDARY ID IS NOT AVAILABLE
  
  660 rjl1 = rj
  zjl1 = zj
  GO TO 680
  
!     FIND SECONDARY ID ENTRY
  
  670 kid = z(i+4)
  CALL bisloc (*2380,kid,z(ibdyl),7,entrys,point)
  ntemp = ibdyl + point - 1
  rjl1  = rz(ntemp+1)
  zjl1  = rz(ntemp+2)
  
!     SECONDARY ID ON PLUS SIDE
  
  680 IF (z(i+5) < 0.0) THEN
    GO TO   690
  ELSE IF (z(i+5) == 0.0) THEN
    GO TO   700
  ELSE
    GO TO   710
  END IF
  
!     SECONDARY ID IS AXIS
  
  690 rjp1 = 0
  zjp1 = zj
  GO TO 720
  
!     SECONDARY ID IS NOT AVAILABLE
  
  700 rjp1 = rj
  zjp1 = zj
  GO TO 720
  
!     FIND SECONDARY ID ENTRY
  
  710 kid = z(i+5)
  CALL bisloc (*2380,kid,z(ibdyl),7,entrys,point)
  ntemp = ibdyl + point - 1
  rjp1  = rz(ntemp+1)
  zjp1  = rz(ntemp+2)
  
!     COMPUTE AND INSERT L,C,S.
  
  720 IF (rj /= 0.0) GO TO 740
  nogo = .true.
  WRITE  (output,730) ufm,z(i)
  730 FORMAT (a23,' 4044, RINGFL CARD ID =',i20,' HAS SPECIFIED A ',  &
      'ZERO RADIAL LOCATION.')
  CYCLE
  
  740 temp1 = rjp1 - rj
  temp2 = 0.25/rj
  r = 0.5*(rjp1-rjl1+temp2*(temp1*temp1-(rjl1-rj)**2))
  zz= 0.5*(zjl1-zjp1+temp2*(temp1*(zj-zjp1)-(rj-rjl1)*(zjl1-zj)))
  rz(i+3) = SQRT(r*r + zz*zz)
  IF (rz(i+3) /= 0.0) GO TO 760
  nogo = .true.
  WRITE  (output,750) ufm,z(i)
  750 FORMAT (a23,' 4045, THE BOUNDARY LIST ENTRY FOR ID =',i9,  &
      ' HAS A ZERO CROSS-SECTION LENGTH.')
  CYCLE
  
  760 rz(i+4) = zz/rz(i+3)
  rz(i+5) =  r/rz(i+3)
END DO

!     SORT GRIDB IMAGES TO BE IN SORT ON RID AND PHI WITHIN EACH RID

780 ntemp = ngridb - igridb + 1
CALL sort (0,0,5,-2,z(igridb),ntemp)
CALL sort (0,0,5,-5,z(igridb),ntemp)

!     THE BOUNDARY FLUID DATA IS ADDED TO THE MATPOOL DATA BLOCK AS 1
!     LOCATE RECORD CONTAINING THE FOLLOWING.

!     1-3   LOCATE CODE  9614,96,0
!     4     CDF
!     5     G
!     6     DRHO
!     7     BD
!     8     NOSYM
!     9     M
!     10    S1
!     11    S2
!     12    N = NUMBER OF INDICES FOLLOWING
!     12+1  THRU  12+N  THE INDICES
!     13+N TO THE EOR IS THE BOUNDARY FLUID DATA


FILE  = matpol
iname = nbdyl + 1
nname = nbdyl
CALL ifp4c (matpol,scrt2,z(buf4),z(buf5),mateof)
IF (mateof) GO TO 930

!     IF ANY DMIAX CARDS ARE PRESENT THEN THEY ARE MERGED IN FRONT OF
!     DMIG CARDS IN THE DMIG RECORD.  FILE NAMES MAY NOT BE THE SAME ON
!     BOTH DMIG AND DMIAX CARDS.

CALL ifp4f (dmiax(2),matpol,bit)
CALL ifp4f (dmig(2) ,matpol,bit2)

!     LOCATE DMIAX CARDS, COPY THEM TO SCRT2 AS DMIG CARDS AND KEEP
!     LIST OF THEIR FILE NAMES.

IF (.NOT.bit .AND. .NOT.bit2) GO TO 900
CALL CLOSE (matpol,clsrew)
CALL preloc (*2360,z(buf4),matpol)

!     WRITE DMIG HEADER.

buf(1) = dmig(1)
buf(2) = dmig(2)
buf(3) = 120
CALL WRITE (scrt2,buf,3,noeor)
IF (.NOT.bit) GO TO 850
CALL locate (*850,z(buf4),dmiax,flag)
ASSIGN 800 TO iretrn

!     READ 9 WORD HEADER

790 GO TO iretrn(800,860)
800 CALL READ (*2340,*850,matpol,buf,9,noeor,flag)

!     SAVE NAME

z(iname  ) = buf(1)
z(iname+1) = buf(2)
nname = nname + 2
icrq  = nname + 2 - icore
IF (icrq > 0) GO TO 2370
810 CALL WRITE (scrt2,buf,9,noeor)

!     COPY THE COLUMN DATA.  FIRST THE COLUMN INDEX.

820 CALL READ (*2340,*2350,matpol,buf,2,noeor,flag)
CALL WRITE (scrt2,buf,2,noeor)
IF (buf(1) < 0.0) THEN
  GO TO   790
END IF

!     TERMS OF COLUMN

830 CALL READ (*2340,*2350,matpol,buf,2,noeor,flag)
CALL WRITE (scrt2,buf,2,noeor)
IF (buf(1) < 0.0) THEN
  GO TO   820
END IF
840 CALL READ (*2340,*2350,matpol,buf,1,noeor,flag)
CALL WRITE (scrt2,buf,1,noeor)
GO TO 830

!     DMIAX-S ALL COPIED.  NOW COPY ANY DMIG-S.

850 IF (.NOT.bit2) GO TO 890
CALL locate (*890,z(buf4),dmig,flag)
ASSIGN 860 TO iretrn

!     READ HEADER

860 CALL READ (*2320,*890,matpol,buf,9,noeor,flag)

!     CHECK THE NAME FOR BEING THE SAME AS ONE ON A DMIAX CARD

DO  i = iname,nname,2
  IF (buf(1) /= z(i  )) CYCLE
  IF (buf(2) /= z(i+1)) CYCLE
  
!     ERROR FOR NAME DOES MATCH THAT OF A DMIAX NAME
  
  nogo = .true.
  WRITE  (output,870) ufm,buf(1),buf(2)
  870 FORMAT (a23,' 4062, DMIG BULK DATA CARD SPECIFIES DATA BLOCK ',  &
      2A4,' WHICH ALSO APPEARS ON A DMIAX CARD.')
END DO

!     COPY THE COLUMN DATA

GO TO 810

!     WRITE THE END OF RECORD FOR DMIG CARDS

890 CALL WRITE (scrt2,0,0,eor)

!     TURN ON BIT FOR DMIG CARD TYPE

CALL ifp4g  (dmig(2),matpol)
CALL REWIND (matpol)
CALL fwdrec (*2340,matpol)

!     COPY EVERYTHING ON MATPOL TO SCRT2, EXCEPT FOR DMIG, DMIAX, AND
!     THE 2**31-1 RECORD.

900 CALL READ (*930,*2350,matpol,buf,3,noeor,flag)
!     2147483647  = 2**31-1
itwo31  =  2147483647
IF (buf(1) /= itwo31.AND.(buf(1) /= dmig(1).OR.buf(2) /= dmig(2))  &
    .AND.(buf(1) /= dmiax(1).OR.buf(2) /= dmiax(2))) GO TO 910
CALL fwdrec (*2340,matpol)
GO TO 900
910 CALL READ  (*2340,*920,matpol,z(nbdyl+1),core-nbdyl,noeor,flag)
CALL WRITE (scrt2,z(nbdyl+1),core-nbdyl,noeor)
GO TO 900
920 CALL WRITE (scrt2,z(nbdyl+1),flag,eor)
GO TO 900
930 mateof = .true.
CALL ifp4b (matpol,scrt2,any,z(nbdyl+1),core-nbdyl,bndfl,mateof)
card(1) = 0
card(2) = 0
card(3) = 0
card(4) = n
CALL locate (*940,z(buf1),flsym,flag)
CALL READ (*2320,*2330,axic,card,3,eor,flag)
940 CONTINUE
CALL WRITE (scrt2,z(1),5,noeor)
CALL WRITE (scrt2,card,4,noeor)
CALL WRITE (scrt2,z(ii),n,noeor)

!     OUTPUT ENTRIES TO MATPOOL DATA BLOCK.(TEMPORARILY ON SCRT2)

jgridb = igridb
jsave  = 0
DO  i = ibdyl,nbdyl,7
  
!     POSSIBILITY OF 2 FLUID ID-S HAVING SAME VALUE
  
  IF (jsave /= 0) jgridb = jsave
  jsave = 0
  IF (z(i) == z(i+7)) jsave = jgridb
  
!     IF RHO FOR A FLUID POINT IS ZERO WE DO NOT PUT OUT FLUID
!     DATA AND CONNECTED POINTS.
  
  IF (rz(i+6) == 0.0) THEN
    GO TO   960
  END IF
  950 CALL WRITE (scrt2,z(i),7,noeor)
  
!     APPEND GRIDB POINTS WITH THEIR ANGLES.
  
  960 IF (jgridb > ngridb) GO TO 1010
  IF (z(jgridb+4) - z(i) < 0.0) THEN
    GO TO   970
  ELSE IF (z(jgridb+4) - z(i) == 0.0) THEN
    GO TO   980
  ELSE
    GO TO  1010
  END IF
  970 jgridb = jgridb + 5
  GO TO 960
  
!     APPEND THE POINT
  
  980 IF (rz(i+6) == 0.0) THEN
    GO TO  1000
  END IF
  990 CALL WRITE (scrt2,z(jgridb),2,noeor)
  1000 jgridb = jgridb + 5
  GO TO 960
  
!     COMPLETE THE ENTRY
  
  1010 IF (rz(i+6) == 0.0) THEN
    GO TO  1030
  END IF
  1020 CALL WRITE (scrt2,mones,2,noeor)
END DO

!     COMPLETE RECORD.

CALL WRITE (scrt2,0,0,eor)
CALL ifp4b (matpol,scrt2,any,z(ngridb+1),core-ngridb,mones,mateof)

!  READ ALL RINGFL CARD IMAGES INTO CORE

1040 IF (anygb) GO TO 1060
IF (.NOT.anygrd) GO TO 1060

!     COPY GRID CARDS NOT COPIED AS A RESULT OF THE ABSENCE OF GRIDB
!     CARDS.

FILE = geom1
1050 CALL READ (*2340,*1060,geom1,card,8,noeor,flag)
CALL WRITE (scrt1,card,8,noeor)
GO TO 1050
1060 iring = ngridb + 1
nring = ngridb
CALL locate (*1090,z(buf1),ringfl,flag)
CALL READ (*2320,*1080,axic,z(iring),core-iring,noeor,flag)
WRITE  (output,1070) ufm
1070 FORMAT (a23,' 4047, INSUFFICIENT CORE TO HOLD RINGFL IMAGES.')
icrq = core - iring
WRITE (output,20) icrq
GO TO 2310
1080 nring = iring + flag - 1

!     OUTPUT HARMONIC GRID CARDS.

1090 IF (nring < iring) GO TO 1150

!     SORT RINGFL IDS

CALL sort (0,0,4,1,z(iring),flag)
card(2)  = 0
rcard(5) = 0.0

!     CARD(6) = -1 AS A FLAG TO TELL GP1 THIS IS A ONE DEGREE OF
!     FREEDOM POINT.

card(6) = -1
card(7) = 0
card(8) = 0
loop1140:  DO  i = ii,ni
  INDEX = z(i)*500000
  DO  k = iring,nring,4
    
!     CALL IFP4E TO CHECK ID RANGE 1 TO 99999
    
    CALL ifp4e (z(k))
    IF (k == iring) GO TO 1100
    IF (z(k) /= ztemp) GO TO 1100
    nogo = .true.
    WRITE (output,420) ufm,z(k)
    1100 ztemp   = z(k)
    card(1) = z(k) + INDEX
    IF (corsys == 3) GO TO 1110
    card(3) = z(k+1)
    card(4) = z(k+3)
    GO TO 1120
    1110 angle = rz(k+2)*degrad
    rcard(3) = rz(k+1)*SIN(angle)
    rcard(4) = rz(k+1)*COS(angle)
    IF (rcard(3) /= 0.0) GO TO 1120
    nogo = .true.
    WRITE (output,350) ufm,z(k)
    CYCLE loop1140
    1120 CALL WRITE (scrt1,card,8,noeor)
  END DO
END DO loop1140

!     COMPLETE GRID CARD RECORD.

1150 CALL WRITE (scrt1,0,0,eor)

!     CREATE AND OUTPUT SEQGP CARDS ONTO SCRT1.  COPY GEOM1 TO SCRT1 UP
!     TO AND INCLUDING SEQGP 3-WORD HEADER.

IF (nring < iring) GO TO 1210
CALL ifp4b (geom1,scrt1,any,z(nring+1),core-nring,seqgp,g1eof)

!     COPY ALL SEQGP CARDS OVER ALSO (ID-S MUST BE OF CORRECT VALUE).

FILE = geom1
IF (.NOT.any) GO TO 1170
1160 CALL READ (*2340,*1170,geom1,card,2,noeor,flag)
CALL ifp4e (card(1))
CALL WRITE (scrt1,card,2,noeor)
GO TO 1160

!     NOW OUTPUT SEQGP CARDS FOR HARMONICS OF EACH RINGFL.

1170 DO  i = ii,ni
  INDEX = z(i)*500000
  ntemp = z(i) - 1
  DO  k = iring,nring,4
    card(1) = z(k) + INDEX
    card(2) = z(k)*1000 + ntemp
    CALL WRITE (scrt1,card,2,noeor)
  END DO
END DO
1200 CALL WRITE (scrt1,0,0,eor)

!     COPY BALANCE OF GEOM1 TO SCRT1 (IF ANY MORE, WRAP UP, AND COPY
!     BACK)

1210 CALL ifp4b(geom1,scrt1,any,z(nring+1),core-nring,mones,g1eof)

!     IF THERE ARE NO HARMONICS THEN ONLY GRID CARDS ARE CREATED FROM
!     GRIDB CARDS.

!     IF (.NOT. HARMS) GO TO 2300
! === IF (.NOT. HARMS) SHOULD NOT GO TO 2300 HERE === G.CHAN/UNISYS 86

!     END OF GEOM1 PROCESSING

!     BEGIN GEOM2 PROCESSING
!     **********************

!     OPEN GEOM2, AND SCRT1. COPY HEADER FROM GEOM2 TO SCRT1.

CALL ifp4c (geom2,scrt1,z(buf2),z(buf3),g2eof)

!     PROCESS CFLUID2, CFLUID3, AND CFLUID4 CARDS.

DO  i = 1,3
  i2 = 2*i
  CALL locate (*1410,z(buf1),cfluid(i2-1),flag)
  
!     COPY DATA FROM GEOM2 TO SCRT1 UP TO POINT WHERE CFLUID CARDS GO
!     AND WRITE 3-WORD RECORD ID.
  
  CALL ifp4b (geom2,scrt1,any,z(ni+1),core-ni,cfluid(2*i-1),g2eof)
  1300 CALL READ (*2320,*1400,axic,card,i+4,noeor,flag)
  IF (card(i+3) /= 1) GO TO 1330
  IF (idrho     /= 1) GO TO 1320
  nogo = .true.
  WRITE  (output,1310) ufm,card(1)
  1310 FORMAT (a23,' 4058, THE FLUID DENSITY HAS NOT BEEN SPECIFIED ON ',  &
      'A CFLUID CARD WITH ID =',i10, /5X,  &
  'AND THERE IS NO DEFAULT ON THE AXIF CARD.')
  1320 rcard(i+3) = drho
  1330 IF (card(i+4) /= 1) GO TO 1360
  IF (ibd /= 1) GO TO 1350
  nogo = .true.
  WRITE  (output,1340) ufm,card(1)
  1340 FORMAT (a23,' 4059, THE FLUID BULK MODULUS HAS NOT BEEN SPECIFIED'  &
      ,      ' ON A CFLUID CARD WITH ID =',i10, /5X,'AND THERE IS NO ',  &
  'DEFAULT ON THE AXIF CARD.')
  1350 rcard(i+4) = bd
  
!     OUTPUT N IMAGES.
  
  1360 ntemp = i+2
  DO  k = 1,ntemp
    saveid(k) = card(k)
  END DO
  
  DO  k = ii,ni
    card(1) = saveid(1)*1000 + z(k)
    INDEX   = 500000*z(k)
    DO  l = 2,ntemp
      card(l) = saveid(l) + INDEX
    END DO
    card(ntemp+3) = (z(k)-1)/2
    CALL WRITE (scrt1,card,ntemp+3,noeor)
  END DO
  GO TO 1300
  
!     END OF CFLUID DATA
  
  1400 CALL WRITE (scrt1,0,0,eor)
END DO

!     CONSTRUCTION OF FSLIST TABLE IN CORE 3-WORDS/ENTRY

ifslst = ni + 1
nfslst = ni
CALL locate (*1600,z(buf1),fslst,flag)
1420 CALL READ (*2320,*1490,axic,rhob,1,noeor,flag)
IF (irhob /= 1) GO TO 1450
IF (idrho /= 1) GO TO 1440
nogo = .true.
WRITE  (output,1430) ufm
1430 FORMAT (a23,' 4048, THE FLUID DENSITY HAS NOT BEEN SPECIFIED ON ',  &
    'AN FSLIST CARD AND', /5X,'THERE IS NO DEFAULT FLUID ',  &
'DENSITY SPECIFIED ON THE AXIF CARD.')
rhob = 1.0
GO TO 1450
1440 rhob = drho
1450 CALL READ (*2320,*2330,axic,idf,1,noeor,flag)
IF (idf == 0) idf = -1
1460 CALL READ (*2320,*2330,axic,idfaft,1,noeor,flag)
IF (idfaft == -1) idfaft = -2
IF (idfaft ==  0) idfaft = -1
IF (nfslst+3 <= core) GO TO 1480
WRITE  (output,1470) ufm
1470 FORMAT (a23,' 4049, INSUFFICIENT CORE TO BUILD FREE SURFACE ',  &
    'LIST TABLE.')
icrq = nfslst + 3 - core
WRITE (output,20) icrq
GO TO 2310
1480 z(nfslst+1) = idf
z(nfslst+2) = idfaft
rz(nfslst+3)= rhob
nfslst = nfslst + 3
IF (idfaft == -2) GO TO 1420
idf = idfaft
GO TO 1460

!     TABLE IS COMPLETE. COPY GEOM2 DATA TO SCRT1 UP TO CFSMASS RECORD
!     SLOT

1490 IF (nfslst > ifslst) GO TO 1510
nogo = .true.
WRITE  (output,1500) ufm
1500 FORMAT (a23,' 4050, FSLIST CARD HAS INSUFFICIENT IDF DATA, OR ',  &
    'FSLIST DATA MISSING.')
GO TO 1600
1510 CALL ifp4b(geom2,scrt1,any,z(nfslst+1),core-nfslst,cfsmas,g2eof)
entrys =(nfslst-ifslst+1)/3
k = 0
DO  i = ifslst,nfslst,3
  IF (z(i+1) == -2) CYCLE
  k = k + 1000000
  rcard(4) = rz(i+2)*g
  DO  l = ii,ni
    INDEX = 500000*z(l)
    card(1) = k + z(l)
    card(2) = z(i) + INDEX
    IF (z(i) <= 0) card(2) = z(i+1) + INDEX
    card(3) = z(i+1) + INDEX
    IF (z(i+1) <= 0) card(3) = z(i) + INDEX
    card(5) = (z(l)-1)/2
    CALL WRITE (scrt1,card,5,noeor)
  END DO
END DO
CALL WRITE (scrt1,0,0,eor)

!     BEGIN GEOM4 PROCESSING
!     **********************

!     OPEN GEOM4 AND SCRT2 AND COPY HEADER RECORD FROM GEOM4 TO SCRT2.

1600 CALL ifp4c (geom4,scrt2,z(buf4),z(buf5),g4eof)

!     COPY ALL DATA ON GEOM4 TO SCRT2 UP TO AND INCLUDING 3-WORD RECORD
!     HEADER OF MPC-RECORD.

CALL ifp4b (geom4,scrt2,any,z(nfslst+1),core-nfslst,mpc,g4eof)

!     COPY ANY MPC IMAGES HAVING A SET ID .LT. 103 TO SCRT2. ERROR
!     MESSAGE IF ANY HAVE ID = 102.  MAINTAIN A LIST OF SETID-S LESS
!     THAN 102.

impc = nfslst + 1
nmpc = nfslst
idlast = 0
FILE = geom4
set 102 = .false.
IF (.NOT.any) GO TO 1650

!     PICK UP SET ID

1610 CALL READ (*2340,*1650,geom4,id,1,noeor,flag)
IF (id > 102) GO TO 1660
IF (id /= 102) GO TO 1630
nogo = .true.
WRITE  (output,1620) ufm
1620 FORMAT (a23,' 4051, AN MPC CARD HAS A SET ID SPECIFIED = 102. ',  &
    ' SET 102 IS ILLEGAL WHEN FLUID DATA IS PRESENT.')
1630 CALL WRITE (scrt2,id,1,noeor)

!     ADD ID TO LIST IF NOT IN LIST

IF (id == idlast) GO TO 1640
nmpc = nmpc + 1
z(nmpc) = id

!     3 WORD GROUPS

1640 CALL READ (*2340,*2350,geom4,card,3,noeor,flag)
CALL WRITE (scrt2,card,3,noeor)
IF (card(1) == -1) GO TO 1610
GO TO 1640

!     NOW POSITIONED TO OUTPUT MPC CARDS FOR SET 102

1650 id = 0

!     IF G FROM AXIF CARD IS NON-ZERO FREEPT DATA IS NOW PROCESSED.

1660 ispnt = nmpc + 1
nspnt = nmpc
press = .false.
IF (g == 0.0) GO TO 1780

!     IF THERE IS NO FREE SURFACE LIST, FREEPT CARDS ARE NOT USED.

IF (nfslst < ifslst) GO TO 1780
CALL sort (0,0,3,1,z(ifslst),nfslst-ifslst+1)
CALL locate (*1780,z(buf1),freept,flag)

!     PICK UP A 3-WORD FREEPT OR PRESPT IMAGE (IDF,IDP,PHI)

1670 CALL READ (*2320,*1770,axic,card,3,noeor,flag)

!     START MPC CARD

angle = rcard(3)*degrad
idf   = card(1)
card(1) = 102
card(3) = 0
IF (press) GO TO 1700

!     LOOK UP RHOB IN FSLIST TABLE

CALL bisloc (*1680,idf,z(ifslst),3,entrys,point)
ntemp = ifslst + point + 1
rcard(4) = -ABS(rz(ntemp)*g)
GO TO 1710
1680 nogo = .true.
WRITE  (output,1690) ufm,idf
1690 FORMAT (a23,' 4052, IDF =',i10,' ON A FREEPT CARD DOES NOT ',  &
    'APPEAR ON ANY FSLIST CARD.')
GO TO 1710
1700 rcard(4) = -1.0
1710 CALL WRITE (scrt2,card,4,noeor)
set102 = .true.

!     ADD SPOINT TO CORE LIST

IF (nspnt+1 <= core) GO TO 1730
WRITE  (output,1720) ufm
1720 FORMAT (a23,' 4053, INSUFFICIENT CORE TO PERFORM OPERATIONS ',  &
    'REQUIRED AS A RESULT OF FREEPT OR PRESPT DATA CARDS')
icrq = nspnt + 1 - core
WRITE (output,20) icrq
GO TO 2310
1730 nspnt = nspnt + 1
z(nspnt) = card(2)
card(2)  = 0

!     HARMONIC COEFFICIENT DATA

DO  i = ii,ni
  card(1) = 500000*z(i) + idf
  nn = (z(i)-1)/2
  IF (MOD(z(i),2) == 0) GO TO 1740
  rcard(3) = SIN(FLOAT(nn)*angle)
  GO TO 1750
  1740 rcard(3) = COS(FLOAT(nn)*angle)
  1750 CALL WRITE (scrt2,card,3,noeor)
END DO
CALL WRITE (scrt2,mones,3,noeor)
GO TO 1670

!     CREATE MPC CARDS AND SPOINTS AS A RESULT OF PRESPT DATA.

1770 IF (press) GO TO 1790
1780 CALL locate (*1790,z(buf1),prespt,flag)
press = .true.
GO TO 1670

!     ANY SPOINTS IN CORE ARE AT THIS TIME OUTPUT TO GEOM2.

1790 IF (nspnt < ispnt) GO TO 1830

!     COPY DATA FROM GEOM2 TO SCRT1 UP TO AND INCLUDING THE 3-WORD
!     RECORD HEADER FOR SPOINTS

FILE = geom2
CALL ifp4b (geom2,scrt1,any,z(nspnt+1),core-nspnt,spoint,g2eof)
IF (.NOT.any) GO TO 1820
1800 CALL READ (*2340,*1810,geom2,z(nspnt+1),core-nspnt,noeor,flag)
CALL WRITE (scrt1,z(nspnt+1),core-nspnt,noeor)
GO TO 1800
1810 CALL WRITE (scrt1,z(nspnt+1),flag,noeor)
1820 CALL WRITE (scrt1,z(ispnt),nspnt-ispnt+1,eor)

!     COPY BALANCE OF GEOM2 TO SCRT1,CLOSE THEM, AND SWITCH DESIGNATIONS

1830 CALL ifp4b (geom2,scrt1,any,z(nmpc+1),core-nmpc,-1,g2eof)

!     END OF GEOM2 PROCESSING
!     ***********************

!     COPY BALANCE OF MPC IMAGES ON GEOM4 TO SCRT2, COMPLETE LIST OF MPC
!     SETS.

FILE = geom4
IF (id == 0) GO TO 1930
GO TO 1910

!     3-WORD GROUPS

1900 CALL READ (*2340,*2350,geom4,card,3,noeor,flag)
CALL WRITE (scrt2,card,3,noeor)
IF (card(1) /= -1) GO TO 1900
CALL READ (*2340,*1930,geom4,id,1,noeor,flag)
1910 IF (id == idlast) GO TO 1920

!     ADD ID TO LIST

idlast = id
nmpc   = nmpc + 1
z(nmpc)= id
1920 CALL WRITE (scrt2,id,1,noeor)
GO TO 1900
1930 CALL WRITE (scrt2,0,0,eor)
TYPE(1) = mpcadd(1)
TYPE(2) = mpcadd(2)

!     GENERATION OF MPCADD OR SPCADD CARDS FROM USER ID-S.  FIRST
!     OUTPUT MANDATORY MPCADD OR SPCADD.

1940 CALL ifp4f (TYPE(2),geom4,bit)
IF (.NOT.set102 .AND. nmpc < impc .AND. .NOT.bit) GO TO 2020
CALL ifp4b (geom4,scrt2,any,z(nmpc+1),core-nmpc,TYPE,g4eof)
IF (.NOT. set102) GO TO 1950
card(1) = 200000000
card(2) = 102
card(3) = -1
CALL WRITE (scrt2,card,3,noeor)

!     NOW FROM USER ID-S

1950 IF (nmpc < impc) GO TO 1980
DO  i = impc,nmpc
  card(1) = z(i) + 200000000
  card(2) = z(i)
  nn = 3
  IF (.NOT.set102) GO TO 1960
  card(3) = 102
  nn = 4
  1960 card(nn) = -1
  CALL WRITE (scrt2,card,nn,noeor)
END DO

!     IF USER MPCADD OR SPCADD CARDS ARE PRESENT, NOW CHANGE THEIR ID-S
!     AND ADD THE 102 SET IF IT EXISTS.

1980 IF (.NOT.any) GO TO 2010
1990 CALL READ (*2340,*2010,geom4,id,1,noeor,flag)
id = id + 200000000
CALL WRITE (scrt2,id,1,noeor)
IF (set102) CALL WRITE (scrt2,102,1,noeor)
2000 CALL READ  (*2340,*2350,geom4,id,1,noeor,flag)
CALL WRITE (scrt2,id,1,noeor)
IF (id == -1) GO TO 1990
GO TO 2000

2010 CALL WRITE (scrt2,0,0,eor)
2020 IF (TYPE(1) == spcadd(1)) GO TO 2270

!     START LIST OF SPC AND SPC1 ID-S

ispc = nfslst + 1
nspc = nfslst
set102 = .false.
idlast = 0

!     CHECK BIT FOR SPC CARDS

CALL ifp4f (spc(2),geom4,bit)
IF (.NOT.bit) GO TO 2080

!     COPY GEOM4 TO SCRT2 UP TO SPC CARDS

CALL ifp4b (geom4,scrt2,any,z(ispc),core-ispc,spc,g4eof)

!     COPY SPC IMAGES KEEPING LIST OF ID-S.

2030 CALL READ (*2340,*2070,geom4,id,1,noeor,flag)
IF (id == idlast) GO TO 2060
IF (id /=    102) GO TO 2050
nogo = .true.
WRITE  (output,2040) ufm
2040 FORMAT (a23,' 4055, SET ID = 102 MAY NOT BE USED FOR SPC CARDS ',  &
    'WHEN USING THE HYDROELASTIC-FLUID ELEMENTS.')
GO TO 2060
2050 nspc = nspc + 1
z(nspc) = id
idlast  = id
2060 CALL WRITE (scrt2,id,1,noeor)
CALL READ  (*2340,*2350,geom4,card,3,noeor,flag)
CALL WRITE (scrt2,card,3,noeor)
GO TO 2030
2070 CALL WRITE (scrt2,0,0,eor)

!     CHECK FOR ANY SPC1 IMAGES

2080 CALL ifp4f (spc1(2),geom4,bit)
IF (.NOT.bit .AND. g /= 0.0) GO TO 2260

!     COPY FROM GEOM4 TO SCRT2 UP TO SPC1 DATA.

CALL ifp4b(geom4,scrt2,any,z(nspc+1),core-nspc-2,spc1,g4eof)

!     COPY SPC1-S UP TO SETID .GE. 103.  SET 102 IS ILLEGAL FOR USER.

IF (.NOT.bit) GO TO 2150
2090 CALL READ (*2340,*2150,geom4,id,1,noeor,flag)
IF (id. LT. 102) GO TO 2100
IF (id /= 102) GO TO 2160
nogo = .true.
WRITE (output,2040) ufm

!     ADD ID TO LIST IF NOT YET IN LIST

2100 IF (nspc < ispc) GO TO 2120
DO  i = ispc,nspc
  IF (id == z(i)) GO TO 2130
END DO

!     ADD ID TO LIST

2120 nspc = nspc + 1
z(nspc) = id
2130 CALL WRITE (scrt2,id,1,noeor)
CALL READ  (*2340,*2350,geom4,id,1,noeor,flag)
CALL WRITE (scrt2,id,1,noeor)
2140 CALL READ  (*2340,*2350,geom4,id,1,noeor,flag)
CALL WRITE (scrt2,id,1,noeor)
IF (id == -1) GO TO 2090
GO TO 2140

!     IF G IS ZERO AND THERE ARE FSLST ENTRIES, GENERATE SPC1-S NOW.

2150 id = 0
2160 IF (g /= 0.0 .OR. nfslst < ifslst) GO TO 2190

!     GENERATION OF HARMONIC SPC1-S

DO  i = ifslst,nfslst,3
  IF (z(i) == -1) CYCLE
  card(1) = 102
  card(2) = 0
  CALL WRITE (scrt2,card,2,noeor)
  DO  j = ii,ni
    CALL WRITE (scrt2,z(i)+500000*z(j),1,noeor)
  END DO
  CALL WRITE (scrt2,minus1,1,noeor)
END DO
set102 = .true.

!     COMPLETE COPYING OF SPC1 CARDS TO SCRT2 WITH SETID-S .GE. 103

2190 IF (id == 0) GO TO 2250

!     ADD ID TO LIST IF NOT YET IN

IF (nspc < ispc) GO TO 2220
2200 DO  i = ispc,nspc
  IF (id == z(i)) GO TO 2230
END DO

!     ID NOT IN LIST, THUS ADD IT.

2220 nspc = nspc + 1
z(nspc) = id

!     CONTINUE COPYING DATA TO NEXT ID

2230 CALL WRITE (scrt2,id,1,noeor)
2240 CALL READ  (*2340,*2350,geom4,id,1,noeor,flag)
CALL WRITE (scrt2,id,1,noeor)
IF (id /= -1) GO TO 2240
CALL READ  (*2340,*2250,geom4,id,1,noeor,flag)
GO TO 2200

!     END OF SPC1 CARD IMAGES.

2250 CALL WRITE (scrt2,0,0,eor)

!     SORT LIST OF SPC AND SPC1 ID-S

CALL sort (0,0,1,1,z(ispc),nspc-ispc+1)

!     SPCADD WORK (USE MPCADD LOGIC)

2260 TYPE(1) = spcadd(1)
TYPE(2) = spcadd(2)
impc = ispc
nmpc = nspc
GO TO 1940

!     ALL PROCESSING COMPLETE ON GEOM4

2270 CALL ifp4b (geom4,scrt2,any,z(1),core,mones,g4eof)

!     END OF GEOM4 PROCESSING
!     ***********************

!     AXIC FILE NOT IN FIST OR AXIF CARD IS MISSING, THUS DO NOTHING.

2300 CALL CLOSE  (axic,clsrew)
CALL conmsg (msg2,2,0)
RETURN

!     FATAL ERROR NO MORE PROCESSING POSSIBLE

2310 nogo = .true.
GO TO 2300

!     END OF FILE ON AXIC

2320 FILE = axic
GO TO 2340

!     END OF RECORD ON AXIC

2330 FILE = axic
GO TO 2350

!     END OF FILE OR END OF RECORD ON -FILE-, OR FILE NOT IN FIST.

2340 ier = -2
GO TO 2390
2350 ier = -3
GO TO 2390
2360 ier = -1
GO TO 2390
2370 ier = -8
FILE = icrq
GO TO 2390
2380 ier = -37
2390 CALL mesage (ier,FILE,subr)
RETURN
END SUBROUTINE ifp4
