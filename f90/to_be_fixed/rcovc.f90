SUBROUTINE rcovc
     
!     RCOVC COMPUTES REACTION FORCES AND GENERATES OUTPUT DATA BLOCKS
!     FOR DISPLACEMENTS, APPLIED LOADS, AND REACTION FORCES.
 
 LOGICAL :: incore     ,uflag      ,pflag      ,non0       ,  &
qflag      ,END        ,once       ,complx     ,  &
    supres     ,KEEP
INTEGER :: dry        ,step       ,fss        ,rfno       ,  &
    uinms      ,ua         ,rss        ,sysbuf     ,  &
    utypo      ,sof2       ,sof3       ,buf1       ,  &
    buf2       ,casess     ,sof1       ,ougv1      ,  &
    opg1       ,oqg1       ,scr1       ,eqss       ,  &
    soln       ,energy     ,pvec       ,uvec       ,  &
    NAME(2)    ,casecc(2)  ,srd        ,pg         ,  &
    substr(4)  ,iz(2)      ,rc         ,FILE       ,  &
    qa         ,pa         ,namef(2)   ,idbuf(146) ,  &
    disp(3)    ,oload(3)   ,spcf(3)    ,dofs(32)   ,  &
    buf3       ,nfwd(3)    ,comps(3)   ,buf(1)     ,  &
    acce(3)    ,velo(3)    ,buf4       ,scr2       ,  &
    scr6       ,scr7       ,scr8       ,scr3
DIMENSION       mcba(7)    ,rbuf(1)    ,rdbuf(7)   ,DATA(12)
CHARACTER (LEN=1) ::
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm swm*27
COMMON /xmssg / ufm        ,uwm        ,uim        ,sfm        , swm
COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
    rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
    nosort     ,uthres     ,pthres     ,qthres
COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
    buf3       ,buf4       ,sof1       ,sof2       , sof3
COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
    iopt       ,rss(2)     ,energy     ,uimpro     ,  &
    range(2)   ,ireq       ,lreq       ,lbasic
COMMON /system/ sysbuf     ,nout
COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
    rew        ,norew
COMMON /unpakx/ utypo      ,iru        ,nru        ,incu
COMMON /condas/ phi        ,twophi     ,raddeg
COMMON /zzzzzz/ z(1)
EQUIVALENCE     (buf(1)    ,z(1))
EQUIVALENCE     (z(1)      ,iz(1)),    (buf(1)     ,rbuf(1))   ,  &
    (idbuf(1)  ,rdbuf(1))
DATA    casess, ougv1      ,opg1       ,oqg1       ,scr1       ,  &
    pg    , scr3       ,scr6       ,scr7       ,scr8       , scr2  /  &
    101   , 201        ,202        ,203        ,301        ,  &
    105   , 303        ,306        ,307        ,308        , 302   /
DATA    srd   / 1          /
DATA    eqss  , soln       ,uvec       ,pvec       /  &
    4HEQSS, 4HSOLN     ,4HUVEC     ,4HPVEC     /
DATA    NAME               ,casecc                 ,substr     /  &
    4HRCOV, 4HC        ,4HCASE     ,4HCC       ,4HSUBS     ,  &
    4HTRUC, 4HTURE     ,4H         /
DATA    comps / 4HCOMP, 4HONEN     ,4HT        /

!     INITIALIZE

IF (dry < 0) RETURN
sof1 = korsz(z) - lreq - sysbuf + 1
sof2 = sof1 - sysbuf - 1
sof3 = sof2 - sysbuf
buf1 = sof3 - sysbuf
buf2 = buf1 - sysbuf
buf3 = buf2 - sysbuf
buf4 = buf3 - sysbuf
lcore= buf4 - 1
IF (lcore <= 0) GO TO 6313
CALL sofopn (z(sof1),z(sof2),z(sof3))

!     ================================================
!     THIS CARD SHOULD BE ADDED WHEN SDR3 IS FIXED

!     IF (RFNO .EQ. 9) NOSORT = 1

!     ================================================
pa = 0
qa = 0
uflag = .false.
pflag = .false.
qflag = .false.

!     CHECK OUTPUT REQUESTS ON CASESS

CALL gopen (casess,z(buf1),rdrew)
nccrec = 1
FILE   = casess
110 CALL fread (casess,z,2,1)
nccrec = nccrec + 1
IF (iz(1) /= casecc(1) .OR. iz(2) /= casecc(2)) GO TO 110
120 CALL READ (*130,*9003,casess,idbuf,35,1,i)
IF (idbuf(17) /= 0) pflag = .true.
IF (idbuf(20) /= 0) uflag = .true.
IF (idbuf(29) /= 0 .AND. rfno >= 8) uflag = .true.
IF (idbuf(32) /= 0 .AND. rfno >= 8) uflag = .true.
IF (idbuf(35) /= 0) qflag = .true.
IF (pflag .AND. uflag .AND. qflag) GO TO 130
GO TO 120
130 CALL CLOSE (casess,rew)

IF (buf(ireq  ) == 1) uflag = .true.
IF (buf(ireq+1) == 1) pflag = .true.
IF (buf(ireq+2) == 1) qflag = .true.
IF (energy == 0) GO TO 135
uflag = .true.
IF (rfno >= 3 .AND. rfno <= 8) pflag = .true.
IF (rfno >= 3 .AND. rfno <= 8) qflag = .true.
135 CONTINUE

IF (.NOT.(uflag .OR. pflag .OR. qflag)) GO TO 900

!     COMPUTE THE APPLIED STATIC LOADS FOR THE REQUESTED SUBSTRUCTURE
!     IF WE ARE PRINTING THE SOLUTION SUBSTRUCTURE CHECK IF THE LOADS
!     ARE ON A GINO FILE.

IF (rfno == 3) pflag = .false.
IF (.NOT.pflag .AND. (.NOT.qflag .OR. rfno == 3)) GO TO 150
IF (rss(1) /= fss(1) .OR. rss(2) /= fss(2)) GO TO 140
pa = pg
mcba(1) = pg
CALL rdtrl (mcba)
IF (mcba(1) > 0) GO TO 150
140 pa = scr3
CALL rcovsl (rss,pvec,0,scr6,scr7,scr8,pa,z(1),z(1),sof3-1, .false.,rfno)
IF (pa <= 0) pflag = .false.

!     GET THE DISPLACEMENT VECTOR AND IF RIGID FORMAT 8 THEN
!     CALCULATE THE VELOCITIES AND ACCELERATIONS.

150 IF (.NOT.uflag .AND. .NOT.qflag) GO TO 170
mcba(1) = ua
CALL rdtrl (mcba)
IF (mcba(1) > 0) GO TO 160
ua = scr2
CALL mtrxi (ua,rss,uvec,0,rc)
IF (rc == 1) GO TO 160
155 ua = 0
WRITE (nout,63190) swm,rss
uflag = .false.
qflag = .false.
energy = 0

160 IF (rfno /= 8 .OR. .NOT.(uflag.OR.qflag)) GO TO 170
CALL rcovva (ua,0,scr1,0,0,0,rss,z(1),z(1),z(1))
IF (ua <= 0) GO TO 155
ua = scr1

!     COMPUTE THE SPCF REACTIONS IF OUTPUT REQUESTS WERE SPECIFIED

170 IF (qflag) CALL rcovqv
IF (qa <= 0) qflag = .false.

!     OUTPUT PROCESSING


!     IF IOPT IS EQUAL TO ONE THEN THE OUTPUT WILL BE SORTED BY SUBCASE
!     IF EQUAL TO TWO IT WILL BE SORTED BY SUBSTRUCTURE

np = buf(ireq+3)
ns = buf(ireq+4)

!     FIND THE LENGTH AND TYPE OF THE VECTORS TO BE OUTPUT

CALL softrl (rss,uvec,mcba)
nmodes = mcba(2)
nsteps = mcba(2)
IF (rfno == 9) nsteps = nsteps/3
complx = .false.
IF (mcba(5) >= 3) complx = .true.
nword = 1
IF (complx) nword = 2

!     PERFORM GENERAL INITIALIZATION OF OFP ID RECORD

idbuf( 3) = 0
idbuf( 6) = 0
idbuf( 7) = 0
idbuf( 8) = 0
idbuf(10) = 8
IF (complx) idbuf(10) = 14
DO  i = 11,50
  idbuf(i) = 0
END DO

!     INITALIZE THE UNPACK COMMON BLOCK

utypo = 1
IF (complx) utypo = 3
iru   = 1
nru   = mcba(3)
incu  = 1

!     ALLOCATE OPEN CORE

isets = 1
lsets = 100
ivect = isets + lsets
isil  = ivect + (nru*nword)
ieqss = isil  + np
IF (ieqss+2 > lcore) GO TO 6313


!                          OPEN CORE DIAGRAM FOR /RCOVCX/

!                       +----------------------------------+
!          Z(ISETS)     I                                  I
!                       I     CASECC SET INFORMATION       I
!                       I                                  I
!                       +----------------------------------+
!          Z(IVECT)     I                                  I
!                       I     VECTOR TO BE PRINTED         I
!                       I                                  I
!                       +----------------------------------+
!          Z(ISIL )     I                                  I
!                       I     SCALAR INDEX LIST FROM EQSS  I
!                       I                                  I
!                       +----------------------------------+
!          Z(IEQSS)     I                                  I
!                       I     EQSS DATA IN TRIPLES OF      I
!                       I        (1) EXTERNAL GRID ID      I
!                       I        (2) INTERNAL POINT INDEX  I
!                       I        (3) COMPONENT CODE        I
!                       I     DATA FOR EACH BASIC SUB-     I
!                       I     STRUCTURE TERMINATED BY      I
!                       I     THREE (-1)S                  I
!                       I                                  I
!                       I     NOTE  EQSS DATA MAY NOT BE   I
!                       I     IN CORE IF SPILL LOGIC       I
!                       I     INVOKED.                     I
!                       I                                  I
!                       +----------------------------------+
!          Z(ISEQ)      I                                  I
!                       I     SYMMETRY SEQUENCE            I
!                       I                                  I
!                       +----------------------------------+
!          Z(ICOMB)     I                                  I
!                       I     VECTOR CONTRIBUTING TO THE   I
!                       I     LINEAR COMBINATION FOR THE   I
!                       I     SYMMETRY SEQUENCE            I
!                       I                                  I
!                       +----------------------------------+

!     READ SIL FROM EQSS INTO OPEN CORE AT ISIL

CALL sfetch (rss,eqss,srd,rc)
n = ns + 1
CALL sjump (n)
DO  i = 1,np
  CALL suread (z(isil+i-1),1,nwds,rc)
  CALL suread (j,1,nwds,rc)
END DO

!     READ EQSS DATA INTO OPEN CORE AT IEQSS IF IT WILL FIT.  IF IOPT
!     EQUALS 2, READ ONLY ONE GROUP AND PRCESS ONE BASIC SUBSTRUCTURE
!     A TIME.

incore = .false.
neqss  = ieqss + 2
CALL sfetch (rss,eqss,srd,rc)
n = 1
CALL sjump (n)
nss = ns
IF (iopt == 2) nss = 1
iss = 0

!     TOP OF LOOP OVER BASIC SUBSTRUCTURES WHEN PROCESSING ONE AT A TIME

475 iss = iss + 1
k   = lcore - ieqss + 1
j   = ieqss
item= eqss
DO  i = 1,nss
  CALL suread (z(j),k,nwds,rc)
  IF (rc == 3) GO TO 6107
  IF (rc /= 2) GO TO 490
  j = j + nwds
  IF (j+3 > lcore) GO TO 490
  iz(j  ) = -1
  iz(j+1) = -1
  iz(j+2) = -1
  j       = j + 3
  neqss   = j - 1
  k       = k - nwds - 3
  IF (k <= 0) GO TO 490
END DO
incore  = .true.
GO TO 491

!     EQSS WILL NOT FIT IN CORE

490 neqss = ieqss + 2
491 iseq  = neqss + 1

!     WRITE HEADER RECORDS ON OUTPUT DATA BLOCKS AND POSITION BOTH
!     INPUT AND OUTPUT DATA BLOCKS AFTER THE HEADER RECORD

DO  i = 1,3
  SELECT CASE ( i )
    CASE (    1)
      GO TO 495
    CASE (    2)
      GO TO 492
    CASE (    3)
      GO TO 493
  END SELECT
  
!     CHECK DISPLACEMENT VECTOR
  
  495 IF (.NOT.uflag) CYCLE
  in   = ua
  iout = ougv1
  GO TO 494
  
!     CHECK LOAD VECTOR
  
  492 IF (.NOT.pflag) CYCLE
  in   = pa
  iout = opg1
  GO TO 494
  
!     CHECK READTIONS VECTOR
  
  493 IF (.NOT.qflag) CYCLE
  in   = qa
  iout = oqg1
  
!     POSITION FILES
  
  494 CALL gopen (in,z(buf1),rdrew)
  CALL CLOSE (in,norew)
  IF (iss > 1) CYCLE
  CALL OPEN (*496,iout,z(buf2),wrtrew)
  CALL fname (iout,namef)
  CALL WRITE (iout,namef,2,1)
  CALL CLOSE (iout,norew)
  CYCLE
  
!     OUTPUT FILE PURGED - TURN OFF REQUEST FLAG
  
  496 WRITE (nout,63140) swm,iout
  IF (iout == ougv1) uflag = .false.
  IF (iout == opg1 ) pflag = .false.
  IF (iout == oqg1 ) qflag = .false.
END DO

!     SETUP FOR LOOP OVER SUBCASES

isc = 0
DO  i = 1,3
  nfwd(i) = 0
END DO

!     POSITION CASESS TO FIRST CASECC SUBCASE

FILE = casess
CALL OPEN (*9001,casess,z(buf3),rdrew)
DO  i = 1,nccrec
  CALL fwdrec (*9002,casess)
END DO
END = .false.

!     TOP OF LOOP OVER SUBCASES

540 isc   = isc + 1
itype = 1
IF (END) GO TO 596

!     READ OUTPUT REQUESTS FROM CASECC RECORD

CALL READ  (*545,*9003,casess,0,-3,0,nwds)
CALL fread (casess,lid  ,1  ,0)
CALL fread (casess,0    ,-12,0)
CALL fread (casess,oload,3  ,0)
CALL fread (casess,disp ,3  ,0)
CALL fread (casess,0    ,-6 ,0)
CALL fread (casess,acce ,3  ,0)
CALL fread (casess,velo ,3  ,0)
CALL fread (casess,spcf ,3  ,0)
CALL fread (casess,0    ,-1 ,0)

!     SET OUTPUT TYPE AND MEDIA - IF NO REQUEST IN CASE CONTROL
!     THE DEFAULT VALUES ARE REAL AND PRINTER

iform = 1
IF (complx)  iform = 2
IF (disp(2)  == 0) disp(2) = 1
IF (disp(3)  == 0) disp(3) = iform
IF (disp(3)  < 0) nosort  = 1
IF (oload(2) == 0) oload(2)= 1
IF (oload(3) == 0) oload(3)= iform
IF (oload(3) < 0) nosort  = 1
IF (spcf(2)  == 0) spcf(2) = 1
IF (spcf(3)  == 0) spcf(3) = iform
IF (spcf(3)  < 0) nosort  = 1
IF (velo(2)  == 0) velo(2) = 1
IF (velo(3)  == 0) velo(3) = iform
IF (velo(3)  < 0) nosort  = 1
IF (acce(2)  == 0) acce(2) = 1
IF (acce(3)  == 0) acce(3) = iform
IF (acce(3)  < 0) nosort  = 1
GO TO 548

!     END OF CASE CONTROL RECORDS - CHECK IF THIS IS REALLY THE END

545 END = .true.
IF (rfno <= 2) GO TO 860
IF (rfno == 3 .AND. isc > nmodes) GO TO 860
IF (rfno >= 8 .AND. isc > nsteps) GO TO 860
GO TO 596

!     READ TITLE, SUBTITLE, AND LABEL.  WILL REPLACE RIGHTMOST WORDS OF
!     SUBTITLE WITH BASIC SUBSTRUCTURE NAME

548 CALL fread (casess,idbuf(51),96,0)
DO  i = 1,3
  idbuf(i+101) = substr(i)
  idbuf(i+133) = comps(i)
END DO
idbuf(  105) = substr(4)
idbuf(  106) = rss(1)
idbuf(  107) = rss(2)

!     READ SYMMETRY SEQUENCE AND SET INFORMATION

nwds =-1
iz(isets  ) = 0
iz(isets+1) = 0
CALL fread (casess,0,-31,0)
CALL fread (casess,lcc,1,0)
lskip = 167 - lcc
CALL fread (casess,0,lskip,0)
CALL READ (*9002,*590,casess,lseq,1,0,n)
IF (neqss+lseq > lcore) GO TO 6313
IF (lseq > 0)CALL READ (*9002,*590,casess,z(iseq),lseq,0,n)
icomb = iseq + lseq
IF (icomb+nru > lcore) GO TO 6313
CALL READ (*9002,*590,casess,z(isets),lsets,0,nwds)
k = lsets

!     MUST EXPAND SETS PORTION OF OPEN CORE

560 n = lcore - neqss
IF (n > 0) GO TO 570
IF (.NOT. incore) GO TO 6313
incore = .false.
neqss  = ieqss + 2
GO TO 560
570 DO  i = isil,neqss
  iz(lcore-i+1) = iz(neqss-i+1)
END DO
ivect = ivect + n
isil  = isil  + n
ieqss = ieqss + n
neqss = neqss + n
CALL READ (*9002,*580,casess,z(isets+lsets),n,0,nwds)
k     = k + n
GO TO 560
580 nwds  = k + nwds
590 nsets = isets + nwds

!     PROCESS OUTPUT ITYPE

596 once  = .false.
jeqss = ieqss - 3
iskip = 0
IF (itype == 1 .AND. .NOT.uflag) GO TO 855
IF (itype == 2 .AND. .NOT.pflag) GO TO 855
IF (itype == 3 .AND. .NOT.qflag) GO TO 855
IF (itype == 4 .AND. .NOT.uflag) GO TO 855
IF (itype == 5 .AND. .NOT.uflag) GO TO 855

!     FOR EACH BASIC SUBSTRUCTURE CURRENTLY BEING PROCESSED, CONSTRUCT
!     ONE OFP ID AND DATA RECORD PAIR.  THE BASIC LOOP IS ABOVE THE
!     VECTOR PROCESSING BECAUSE OUTPUT REQUESTS CAN CHANGE FOR EACH
!     BASIC

DO  js = 1,nss
  jss  = iss  + js - 1
  nreq = ireq + (jss-1)*lbasic + 5
  kpoint = buf(nreq+12)
  
!     STATICS
  
  IF (rfno > 2) GO TO 603
  IF (js   > 1) GO TO 598
  iappro   = 1
  idbuf(4) = isc
  idbuf(5) = lid
  GO TO 598
  
!     FOR NORMAL MODES GET MODE NUMBER, EIGENVALUE AND FREQUENCY
  
  603 IF (rfno /= 3) GO TO 612
  IF (js   > 1) GO TO 598
  CALL sfetch (fss,soln,srd,rc)
  n = 1
  CALL sjump (n)
  j = isc - 1
  IF (j == 0) GO TO 611
  DO  i = 1,j
    CALL suread (mcba(1),7,nwds,rc)
  END DO
  611 CONTINUE
  CALL suread (mode,1,nwds,rc)
  CALL suread (i,1,nwds,rc)
  CALL suread (eigen ,1,nwds,rc)
  CALL suread (eigeni,1,nwds,rc)
  CALL suread (value,1,nwds,rc)
  
  iappro = 2
  IF (complx) iappro = 9
  idbuf(4) = isc
  idbuf(5) = mode
  rdbuf(6) = eigen
  rdbuf(7) = 0.0
  IF (complx) rdbuf(7) = eigeni
  GO TO 598
  
!     FOR DYNAMICS GET THE TIME OR FREQUENCY
  
  612 IF (rfno /= 8 .AND. rfno /= 9) GO TO 598
  IF (js > 1) GO TO 598
  CALL sfetch (fss,soln,srd,rc)
  n = 1
  CALL sjump (n)
  j = isc - 1
  IF (j == 0) GO TO 614
  DO  i = 1,j
    CALL suread (mcba(1),1,nwds,rc)
  END DO
  614 CONTINUE
  CALL suread (value,1,nwds,rc)
  
  iappro = 5
  IF (rfno == 9) iappro = 6
  idbuf(4) = isc
  rdbuf(5) = value
  idbuf(8) = lid
  
!     GET SUBCASE OR MODE REQUEST
  
  598 IF (rfno > 2) GO TO 599
  isub = isc
  iloc = 5
  GO TO 600
  599 IF (rfno /= 3) GO TO 607
  isub = mode
  iloc = 6
  GO TO 600
  607 isub = isc
  iloc = 11
  600 iset = buf(nreq+iloc)
  IF (iset < 0) GO TO 608
  IF (iset == 0) GO TO 835
  
!     FIND THE REQUESTED SET
  
  jset = isets
  601 CONTINUE
  IF (iset == iz(jset)) GO TO 602
  jset = jset + iz(jset+1) + 2
  IF (jset < nsets) GO TO 601
  
!     SET NOT FOUND, ISSUE WARNING AND PRINT ALL INSTEAD.
  
  WRITE (nout,63650) uwm,iset
  buf(nreq+iloc) = -1
  GO TO 608
  
!     FIND IF CURRENT SUBCASE OR MODE IS IN REQUESTED SET
  
  602 next = 1
  kset = iz(jset+1)
  CALL setfnd (*835,iz(jset+2),kset,isub,next)
  
!     SO FAR SO GOOD - IF NORMAL MODES OR DYNAMICS PROBLEM CHECK IF
!     EIGEN VALUE, TIME OR FREQUENCY IS IN REQUESTED RANGE
  
  608 CONTINUE
  IF (rfno < 3) GO TO 609
  IF (value < rbuf(nreq+7)) GO TO 835
  IF (value > rbuf(nreq+8)) GO TO 835
  
  609 SELECT CASE ( itype )
    CASE (    1)
      GO TO 615
    CASE (    2)
      GO TO 640
    CASE (    3)
      GO TO 650
    CASE (    4)
      GO TO 652
    CASE (    5)
      GO TO 654
  END SELECT
  
!     PROCESS DISPLACEMENT REQUESTS
  
  615 iopst = disp(1)
  IF (buf(nreq+2) > -2) iopst = buf(nreq+2)
  IF (iopst == 0 .AND. lseq == 0) GO TO 835
  IF (once) GO TO 705
  once  = .true.
  
  idc   = disp(2)
  iform = IABS(disp(3))
  idbuf(2) = 1
  IF (rfno == 3) idbuf(2) = 7
  thresh = uthres
  supres = .false.
  in   = ua
  iout = ougv1
  GO TO 660
  
!     PROCESS OLOAD REQUESTS
  
  640 iopst = oload(1)
  IF (buf(nreq+3) > -2) iopst = buf(nreq+3)
  IF (iopst == 0 .AND. lseq == 0) GO TO 835
  IF (once) GO TO 705
  once = .true.
  
  idc    = oload(2)
  iform  = IABS(oload(3))
  thresh = pthres
  supres = .true.
  idbuf(2) = 2
  in   = pa
  iout = opg1
  GO TO 660
  
!     PROCESS SPCFORCE (ACTUALLY, ALL REACTIONS) REQUESTS
  
  650 iopst = spcf(1)
  IF (buf(nreq+4) > -2 ) iopst = buf(nreq+4)
  IF (iopst == 0 .AND. lseq == 0) GO TO 835
  IF (once) GO TO 705
  once = .true.
  
  idc    = spcf(2)
  iform  = IABS(spcf(3))
  thresh = qthres
  supres = .true.
  idbuf(2) = 3
  in   = qa
  iout = oqg1
  GO TO 660
  
!     PROCESS VELOCITY REQUESTS
  
  652 iopst = velo(1)
  IF (buf(nreq+9) > -2) iopst = buf(nreq+9)
  IF (iopst == 0 .AND. lseq == 0) GO TO 835
  IF (once) GO TO 705
  once = .true.
  
  idc   = velo(2)
  iform = IABS(velo(3))
  idbuf(2) = 10
  thresh = uthres
  supres = .false.
  in   = ua
  iout = ougv1
  GO TO 660
  
!     PROCESS ACCELERATION REQUESTS
  
  654 iopst = acce(1)
  IF (buf(nreq+10) > -2) iopst = buf(nreq+10)
  IF (iopst == 0 .AND. lseq == 0) GO TO 835
  IF (once) GO TO 705
  once = .true.
  
  idc   = acce(2)
  iform = IABS(acce(3))
  idbuf(2) = 11
  thresh = uthres
  supres = .false.
  in   = ua
  iout = ougv1
  
!     OPEN FILES AND UNPACK VECTOR TO BE PRINTED
  
  660 FILE = in
  CALL gopen (in,z(buf1),rd)
  CALL gopen (iout,z(buf2),wrt)
  it  = itype
  IF (itype > 3) it = 1
  IF (lseq > 0) GO TO 664
  n   = nfwd(it)
  IF (n <= 0) GO TO 663
  DO  i = 1,n
    CALL fwdrec (*9002,in)
  END DO
  nfwd(it) = 0
  663 CALL unpack (*673,in,z(ivect))
  GO TO 675
  
!     FORM LINEAR COMBINATION FOR SYMMETRY SEQUENCE
  
  664 n = nfwd(it) - lseq
  IF (n < 0) THEN
    GO TO   665
  ELSE IF (n == 0) THEN
    GO TO   669
  ELSE
    GO TO   667
  END IF
  665 n = -n
  DO  i = 1,n
    CALL bckrec(in)
  END DO
  GO TO 669
  667 DO  i = 1,n
    CALL fwdrec (*9002,in)
  END DO
  669 DO  i = 1,nru
    z(ivect+i-1) = 0.0E0
  END DO
  DO  i = 1,lseq
    CALL unpack (*672,in,z(icomb))
    DO  j = 1,nru
      z(ivect+j-1) = z(ivect+j-1) + z(iseq+i-1)*z(icomb+j-1)
    END DO
  END DO
  nfwd(it) = 0
  GO TO 675
  673 n = nru*nword
  DO  i = 1,n
    z(ivect+i-1) = 0.0
  END DO
  
!     IF EQSS DATA NOT IN CORE, POSITION THE SOF
  
  675 IF (incore) GO TO 705
  CALL sfetch (rss,eqss,srd,rc)
  nskip = iss + iskip
  CALL sjump (nskip)
  jeqss = ieqss
  
!     INSERT SUBSTRUCTURE NAME IN IDREC WRITE IT OUT
  
  705 idbuf(1) = idc + 10*iappro
  IF (complx .AND. js == 1) idbuf(2) = idbuf(2) + 1000
  idbuf(  9) = iform
  idbuf(138) = buf(nreq)
  idbuf(139) = buf(nreq+1)
  KEEP = .false.
  
!     FIND THE REQUESTED OUTPUT SET
  
  next  = 1
  jset  = isets
  njset = jset + 1
  IF (iopst < 0) GO TO 730
  710 IF (iopst == iz(jset)) GO TO 730
  jset  = jset + iz(jset+1) + 2
  IF (jset < nsets) GO TO 710
  
!     SET NOT FOUND. ISSUE A WARNING AND PRINT ALL INSTEAD
  
  WRITE (nout,63650) uwm,iopst
  i = itype + 1
  IF (itype > 3) i = i + 4
  buf(nreq+i) = -1
  iopst = -1
  
!     FOR EACH GRID POINT ID IN EQSS FOR THE CURRENT SUBSTRUCTURE WHICH
!     IS A MEMBER OF THE REQUESTED OUTPUT SET, WRITE A LINE OF OUTPUT
  
  730 IF (incore) GO TO 780
  CALL suread (z(jeqss),3,nwds,rc)
  IF (rc /= 1) GO TO 830
  GO TO 790
  780 jeqss = jeqss + 3
  IF (iz(jeqss) > 0) GO TO 790
  GO TO 830
  
  790 IF (iopst < 0) GO TO 800
  IF (next > iz(jset+1)) GO TO 830
  kset = iz(jset+1)
  kid  = iz(jeqss )
  CALL setfnd (*730,iz(jset+2),kset,kid,next)
  
!     WRITE A LINE OF OUTPUT
  
  800 icode = iz(jeqss+2)
  CALL decode (icode,dofs(1),n)
  dofs(n+1) = -1
  jsil = iz(jeqss+1) + isil - 1
  k    = 0
  non0 = .false.
  DO  i = 1,6
    IF (dofs(k+1)+1 /= i) GO TO 815
    j = ivect + (iz(jsil)-1)*nword + k*nword
    k = k + 1
    DATA(i) = z(j)
    IF (complx) GO TO 805
    IF (supres .AND. DATA(i) == 0.0) GO TO 815
    IF (ABS(DATA(i)) < thresh) GO TO 815
    non0 = .true.
    CYCLE
    805 DATA(6+i) = z(j+1)
    IF (iform /= 3 .OR. DATA(i)+DATA(6+i) == 0.0) GO TO 810
    DATA(i)   = SQRT(z(j)**2 + z(j+1)**2)
    DATA(6+i) = ATAN2(z(j+1),z(j))*raddeg
    IF (DATA(6+i) < -.000005) DATA(6+i) = DATA(6+i) + 360.0
    810 IF (supres .AND. DATA(i)+DATA(6+i) == 0.0) GO TO 815
    IF (ABS(DATA(i)) < thresh .AND. ABS(DATA(6+i)) < thresh) GO TO 815
    non0 = .true.
    CYCLE
    815 DATA(  i) = 0.0
    DATA(6+i) = 0.0
  END DO
  IF (.NOT.non0) GO TO 825
  IF (.NOT. KEEP) CALL WRITE (iout,idbuf,146,1)
  CALL WRITE (iout,10*iz(jeqss)+idc,1,0)
  CALL WRITE (iout,kpoint,1,0)
  CALL WRITE (iout,DATA,6*nword,0)
  KEEP = .true.
  825 CONTINUE
  IF (next <= iz(jset+1) .OR. iopst < 0) GO TO 730
  
!     IF NO DATA WAS WRITTEN FOR THIS BASIC BACKREC THE OFP FILE
!     OVER THE PREVIOUSLY WRITTEN ID RECORD
  
  830   IF (KEEP) CALL WRITE (iout,0,0,1)
  IF (iz(jeqss) < 0 .OR. (.NOT.incore .AND. rc /= 1)) CYCLE
  
!     NO MORE OUTPUT FOR THIS BASIC - SKIP EQSS DATA
  
  835 CONTINUE
  IF (incore) GO TO 836
  IF (once  ) GO TO 837
  iskip = iskip + 1
  CYCLE
  837 n = 1
  CALL sjump (n)
  CYCLE
  836 jeqss = jeqss + 3
  IF (iz(jeqss) > 0) GO TO 836
END DO

!     GO BACK AND DO ANOTHER OUTPUT TYPE

CALL CLOSE (in,norew)
CALL CLOSE (iout,norew)
855 IF (once) GO TO 856
it = itype
IF (itype > 3) it = 1
nfwd(it) = nfwd(it) + 1
856 itype = itype + 1
IF (itype <= 3) GO TO 596
IF (itype <= 5 .AND. rfno >= 8) GO TO 596
IF (.NOT.END) GO TO 540
IF (rfno == 3 .AND. isc < nmodes) GO TO 540
IF (rfno >= 8 .AND. isc < nsteps) GO TO 540

!     ALL SUBCASES PROCESSED,  IF IOPT EQ 2, GO BACK AND PROCESS
!     NEXT BASIC SUBSTRUCTURE

860 CALL CLOSE (casess,rew)
IF (iopt == 1 .OR. iss == ns) GO TO 870
CALL sfetch (rss,eqss,srd,rc)
n = iss + 1
CALL sjump (n)
GO TO 475

!     WRITE TRAILERS AND EOF ON OUTPUT DATA BLOCKS

870 DO  i = 2,7
  mcba(i) = 1
END DO
IF (.NOT.uflag) GO TO 885
CALL gopen (ougv1,z(buf1),wrt)
CALL CLOSE (ougv1,rew)
mcba(1) = ougv1
CALL wrttrl (mcba)
885 IF (.NOT.pflag) GO TO 890
CALL gopen (opg1,z(buf1),wrt)
CALL CLOSE (opg1,rew)
mcba(1) = opg1
CALL wrttrl (mcba)
890 IF (.NOT.qflag) GO TO 900
CALL gopen (oqg1,z(buf1),wrt)
CALL CLOSE (oqg1,rew)
mcba(1) = oqg1
CALL wrttrl (mcba)

!     NORMAL MODULE TERMINATION

900 CALL sofcls
RETURN

!     ERROR PROCESSING

6107 n = 7
CALL smsg (n,item,rss)
GO TO 9200
6313 WRITE (nout,63130) swm,rss
GO TO 9200
9001 n = 1
GO TO 9100
9002 n = 2
GO TO 9100
9003 n = 3
GO TO 9100
9100 CALL mesage (n,FILE,NAME)
9200 CALL sofcls
DO  i = 101,111
  CALL CLOSE (i,rew)
END DO
DO  i = 201,203
  CALL CLOSE (i,rew)
END DO
DO  i = 301,308
  CALL CLOSE(i,rew)
END DO
RETURN

!     DIAGNOSTICS FORMAT STATEMENTS

63130 FORMAT (a27,' 6313, INSUFFICIENT CORE FOR RCOVR MODULE WHILE ',  &
    'TRYING TO PROCESS', /34X,'PRINTOUT DATA BLOCKS FOR ', 'SUBSTRUCTURE',2A4)
63140 FORMAT (a27,' 6314, OUTPUT REQUEST CANNOT BE HONORED.', /34X,  &
    'RCOVR MODULE OUTPUT DATA BLOCK',i4,' IS PURGED.')
63190 FORMAT (a27,' 6319, DISPLACEMENT MATRIX FOR SUBSTRUCTURE ',2A4,  &
    ' MISSING.' /5X,'DISPLACEMENT OUTPUT REQUESTS CANNOT BE ',  &
    'HONORED.  SPCFORCE OUTPUT REQUESTS CANNOT BE HONORED UN',  &
    'LESS THE', /5X,'REACTIONS HAVE BEEN PREVIOUSLY COMPUTED.')
63650 FORMAT (a25,' 6365, REQUESTED OUTPUT SET ID',i6,' IS NOT DECLARED'  &
    ,      ' IN CASE CONTROL, ALL OUTPUT WILL BE PRODUCED.')
END SUBROUTINE rcovc
