SUBROUTINE sdr2d
     
!     SDR2D PERFORMS THE FINAL STRESS AND FORCE RECOVERY COMPUTATIONS.
!     CASE CONTROL AND THE DISPLACEMENT VECTOR FILE ARE PROCESSED IN
!     PARALLEL.  THE ESTA IS PASSED ONCE FOR EACH VECTOR IN UGV FOR
!     WHICH A STRESS OR FORCE OUTPUT REQUEST EXISTS.  THE ESTA IS HELD
!     COMPLETELY IN CORE IF POSSIBLE.  STRESS OUTPUT IS WRITTEN ON OES1.
!     FORCE OUTPUT IS WRITTEN ON OEF1.
 
 LOGICAL :: axic     ,axsine   ,axcosi   ,eofcc
!    1,               IDSTRS   ,IDFORC   ,ILOGIC(2)
 INTEGER :: stresx   ,forcex   ,ugvvec   ,estawd   ,eltype  ,  &
     tload    ,eldef    ,elemid   ,gptt     ,oes1    ,  &
     oef1     ,oeigr    ,esta     ,edt      ,elesta  ,  &
     kdefrm(2),app      ,sort2    ,spcf     ,displ   ,  &
     vel      ,acc      ,stress   ,force    ,cstm    ,  &
     casecc   ,eqexin   ,sil      ,bgpdt    ,pg      ,  &
     qg       ,ugv      ,phig     ,eigr     ,opg1    ,  &
     oqg1     ,ougv1    ,ocb      ,buf(50)  ,dtype   ,  &
     FILE     ,buf1     ,buf2     ,buf3     ,buf4    ,  &
     buf5     ,buf6     ,buf7     ,symflg   ,outfl   ,  &
     sta      ,rei      ,ds0      ,ds1      ,frq     ,  &
     trn      ,bk0      ,branch   ,sysbuf   ,  &
     date     ,plots    ,qtype2   ,eol      ,bk1     ,  &
     time     ,setno    ,fsetno   ,z        ,retx    ,  &
     formt    ,flag     ,eof      ,cei      ,pla     ,  &
     bufa     ,bufb     ,ofile    ,device   ,oef1l   ,  &
     pugv1    ,xsetns   ,sdest    ,buf8     ,oes1l   ,  &
     opte     ,xset0    ,xsetnf   ,fdest    ,pcomps  ,  &
     sorc     ,tloads   ,tmprec   ,itr(7)   ,comps
 INTEGER :: pcomp(2) ,pcomp1(2),pcomp2(2),buf0     ,bufm1   ,  &
     nmes1l(2),nmef1l(2)
 REAL :: zz(1)    ,bufr(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 COMMON /BLANK / app(2)   ,sort2    ,idum(2)  ,comps
 COMMON /sdr2c1/ ipcmp    ,npcmp    ,ipcmp1   ,npcmp1   ,ipcmp2  ,  &
     npcmp2   ,nstrop
 COMMON /sdr2x1/ ieigen   ,ieldef   ,itload   ,isymfl   ,iloads  ,  &
     idispl   ,istr     ,ielf     ,iacc     ,ivel    ,  &
     ispcf    ,ittl     ,ilsym    ,ifrout   ,isload  , idload   ,isorc
 COMMON /sdr2x2/ casecc   ,cstm     ,mpt      ,dit      ,eqexin  ,  &
     sil      ,gptt     ,edt      ,bgpdt    ,pg      ,  &
     qg       ,ugv      ,est      ,phig     ,eigr    ,  &
     opg1     ,oqg1     ,ougv1    ,oes1     ,oef1    ,  &
     pugv1    ,oeigr    ,ophig    ,pphig    ,esta    ,  &
     gptta    ,harms    ,xycdb    ,scr3     ,pcomps  , oes1l    ,oef1l
COMMON /sdr2x4/ nam(2)   ,END      ,mset     ,icb(7)   ,ocb(7)  ,  &
    mcb(7)   ,dtype(8) ,icstm    ,ncstm    ,ivec    ,  &
    ivecn    ,temp     ,deform   ,FILE     ,buf1    ,  &
    buf2     ,buf3     ,buf4     ,buf5     ,any     ,  &
    all      ,tloads   ,eldef    ,symflg   ,branch  ,  &
    ktype    ,loads    ,spcf     ,displ    ,vel     ,  &
    acc      ,stress   ,force    ,kwdest   ,kwdedt  ,  &
    kwdgpt   ,kwdcc    ,nrigds   ,sta(2)   ,rei(2)  ,  &
    ds0(2)   ,ds1(2)   ,frq(2)   ,trn(2)   ,bk0(2)  ,  &
    bk1(2)   ,cei(2)   ,pla(22)  ,nrings   ,nharms  ,  &
    axic     ,knset    ,isopl    ,strspt   ,ddrmm   , isopl8
COMMON /sdr2x7/ elesta(100)        ,bufa(100),bufb(4076)
COMMON /sdr2x8/ elwork(300)
COMMON /zzzzzz/ z(1)
COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,clsrew
COMMON /system/ sysbuf   ,opte     ,nogo     ,intap    ,mpcn    ,  &
    spcn     ,method   ,loadnn   ,symm     ,stftmp  ,  &
    page     ,line     ,tline    ,maxlin   ,date(3) ,  &
    time     ,echo     ,plots    ,dum23(35),iheat
COMMON /unpakx/ qtype2   ,i2       ,j2       ,incr2
COMMON /zntpkx/ xx(4)    ,ixx      ,eol      ,eor
COMMON /zblpkx/ y(4)     ,iy
COMMON /sdr2de/ buf6     ,coef     ,deftmp   ,diff     ,diff1   ,  &
    device   ,estawd   ,elemid   ,eltype   ,eof     ,  &
    eofcc    ,ireqx    ,flag     ,fn       ,forcex  ,  &
    fsetno   ,formt    ,icc      ,i        ,iedt    ,  &
    isetno   ,isetf    ,isets    ,idef     ,isymn   ,  &
    sdest    ,ix       ,isetnf   ,iseq     ,iretrn  ,  &
    irecx    ,isave    ,fdest    ,ipart    ,ilist   ,  &
    igptta   ,icore    ,ielem    ,iesta    ,buf8    ,  &
    jforc    ,jstrs    ,jany     ,jlist    ,j       ,  &
    ktype1   ,khi      ,kx       ,k        ,klo     ,  &
    kn       ,ktypex   ,kfrq     ,kcount   ,lsym    ,  &
    m        ,midvec   ,nwdsa    ,nwdstr   ,nlogic  ,  &
    nwds     ,ndef     ,n        ,n1       ,n2      ,  &
    notset   ,nsets    ,nsetf    ,nwords   ,nx      ,  &
    tgrid(4) ,nwdfor   ,ngptt    ,nesta    ,nvects  ,  &
    nlist    ,ofile    ,outfl    ,retx     ,setno   ,  &
    stresx   ,SAVE     ,tload    ,ugvvec   ,ixsets  ,  &
    nxsets   ,ixsetf   ,nxsetf   ,xsetns   ,xsetnf  ,  &
    sorc     ,tmprec   ,buf7     ,tgrd(33)
EQUIVALENCE     (buf(1),bufr(1))   ,(z(1),zz(1))
!    1,               (IDSTRS,ILOGIC(1)) ,(IDFORC,ILOGIC(2))
DATA    buf   / 50*0/, kdefrm/104,1/, xset0/100000000/
DATA    nmes1l/ 4HOES1, 4HL   /    ,nmef1l  / 4HOEF1, 4HL   /
DATA    pcomp /  5502,  55    /, pcomp1/  5602,  56    /,  &
    pcomp2/  5702,  57    /

!     PERFORM GENERAL INITIALIZATION

bufm1  = korsz(z) - sysbuf + 1
buf0   = bufm1 - sysbuf - 1
buf1   = buf0  - sysbuf - 1
IF (comps /= -1) buf1 = bufm1
buf2   = buf1  - sysbuf - 1
i2     = 1
incr2  = 1
icc    = 0
ilist  = 1
nlist  = 0
jlist  = 1
kfrq   = 0
axsine = .false.
axcosi = .false.
sorc   = 0

!     READ TRAILER ON INPUT FILE. SET PARAMETERS.

icb(1) = ugv
CALL rdtrl (icb)
IF (icb(1) /= ugv) GO TO 770
nvects = icb(2)
IF (icb(5) > 2) GO TO 10

!     REAL VECTOR.

ktype  = 1
qtype2 = 1
ktype1 = 2
nwds   = 8
ktypex = 0
GO TO 20

!     COMPLEX VECTOR.

10 ktype  = 2
qtype2 = 3
ktype1 = 3
nwds   = 14
ktypex = 1000

!     OPEN CASE CONTROL AND SKIP HEADER. THEN BRANCH ON APPROACH.

20 FILE  = casecc
CALL OPEN (*740,casecc,z(buf1),rdrew)
CALL fwdrec (*750,casecc)
eofcc = .false.

SELECT CASE ( branch )
  CASE (    1)
    GO TO 100
  CASE (    2)
    GO TO  30
  CASE (    3)
    GO TO 100
  CASE (    4)
    GO TO  60
  CASE (    5)
    GO TO  70
  CASE (    6)
    GO TO  70
  CASE (    7)
    GO TO 100
  CASE (    8)
    GO TO  60
  CASE (    9)
    GO TO  30
  CASE (   10)
    GO TO 100
END SELECT
!            STA,REI,DS0,DS1,FRQ,TRN,BK0,BK1,CEI,PLA

!     EIGENVALUES - READ LIST OF MODE NOS. AND EIGENVALUES INTO CORE.
!     BUCKLING POSSIBLE HERE TOO

30 FILE = eigr
CALL OPEN (*740,eigr,z(buf2),rdrew)
CALL fwdrec (*750,eigr)
CALL fwdrec (*750,eigr)
i = ilist
m = 8 - ktype
iskip = 0
INDEX = 2
IF (app(1) /= rei(1)) GO TO 40

!     CHECK TO SEE IF ALL GENERALIZED MASS VALUES ARE ZERO

35 CALL READ (*750,*37,eigr,buf,m,0,flag)
IF (buf(6) == 0.0) GO TO 35
INDEX = 0
37 CALL skprec (eigr,-1)
40 CALL READ (*750,*50,eigr,buf(1),m,0,flag)
IF (app(1) /= rei(1)) GO TO 45
IF (INDEX == 2) GO TO 45

!     MATCH CORRECT MODE NOS. AND EIGENVALUES WITH PROPER
!     FORCES AND STRESSES WHEN USING GIVENS METHOD WITH F1.GT.0.0

IF (INDEX == 1) GO TO 45
IF (buf(6) /= 0.0) GO TO 43
iskip  = iskip + 1
GO TO 40
43 INDEX  = 1
45 z(i  ) = buf(1) - iskip
z(i+1) = buf(3)
z(i+2) = buf(4)
i      = i + ktype1
GO TO 40
50 CALL CLOSE (eigr,clsrew)
nlist  = i - ktype1
icc    = i
GO TO 100

!     DIFF. STIFF. PHASE 1 OR BUCKLING PHASE 1 - SKIP 1ST DATA RECORD ON
!     CC.

60 CALL fwdrec (*750,casecc)
IF (app(1) == bk1(1)) GO TO 30
GO TO 100

!     FREQUENCY OR TRANSIENT RESPONSE - READ LIST INTO CORE.

70 FILE = pg
CALL OPEN (*740,FILE,z(buf2),rdrew)
i  = ilist
m  = 3
ix = 1
IF (app(1) == frq(1) .OR. app(1) == trn(1)) ix = 2
80 CALL READ (*750,*90,FILE,buf(1),m,0,flag)
z(i  ) = buf(m)
z(i+1) = 0
i = i + ix
m = 1
GO TO 80
90 CALL CLOSE (FILE,clsrew)
nlist = i - ix
icc   = i

!     ALLOCATE CORE FOR CASE CONTROL, EDT, GPTT, ESTA, VECTOR
!     BALANCE OF REQUIRED BUFFERS
!       BUF1 = CASECC     BUF5 = GPTT
!       BUF2 = VECTOR     BUF6 = EDT
!       BUF3 = OES1       BUF7 = EQEXIN
!       BUF4 = OEF1       BUF8 = ESTA
!     SOME OF THE ABOVE MAY NOT BE REQUIRED AND THUS WILL NOT BE
!     ALLOCATED..

100 buf3 = buf2 - sysbuf - 1
IF (stress == 0) buf3 = buf2
buf4 = buf3 - sysbuf - 1
IF (force  == 0) buf4 = buf3
buf5 = buf4 - sysbuf - 1
IF (tloads == 0) buf5 = buf4
buf6 = buf5 - sysbuf - 3
IF (kwdedt == 0) buf6 = buf5
buf7 = buf6 - sysbuf - 1
IF (isopl == 0) buf7 = buf6
buf8 = buf7 - sysbuf - 1

!     IF COMPOSITE ELEMENTS ARE PRESENT, READ PCOMPS INTO CORE

IF (comps /= -1) GO TO 109
FILE   = pcomps
n      = -1
CALL preloc (*760,z(buf2),pcomps)
ipcmp  = icc + 1
ipcmp1 = ipcmp
ipcmp2 = ipcmp
npcmp  = 0
npcmp1 = 0
npcmp2 = 0
n      = -2

CALL locate (*106,z(buf2),pcomp,idx)
CALL READ (*760,*106,pcomps,z(ipcmp),buf2-ipcmp,1,npcmp)
CALL mesage (-8,0,nam)
106 ipcmp1 = ipcmp1 + npcmp
ipcmp2 = ipcmp1

CALL locate (*107,z(buf2),pcomp1,idx)
CALL READ (*760,*107,pcomps,z(ipcmp1),buf2-ipcmp1,1,npcmp1)
CALL mesage (-8,0,nam)
107 ipcmp2 = ipcmp2 + npcmp1

CALL locate (*108,z(buf2),pcomp2,idx)
CALL READ (*760,*108,pcomps,z(ipcmp2),buf2-ipcmp2,1,npcmp2)
CALL mesage (-8,0,nam)
108 icc = ipcmp2 + npcmp2 - 1

CALL CLOSE (pcomps,clsrew)

!     IF ESTA FITS IN CORE BUF8 MAY BE BUF7 SINCE IT WILL ONLY BE USED
!     TO READ ESTA IN ONCE..

109 iedt   = icc + kwdcc + 1
igptta = iedt + kwdedt
itr(1) = eqexin
CALL rdtrl (itr)
neqex  = 2*itr(2)
IF (isopl8 /= 8) neqex = 0
ieqex  = igptta + kwdgpt
ivec   = ieqex + neqex
ivecn  = ivec + ktype*icb(3) - 1

!     IF CONICAL SHELL DOUBLE VECTOR SPACE

IF (axic .AND. ktype == 1) ivecn = ivecn + icb(3)*ktype
iesta  = ivecn + 1
midvec = (ivec+ivecn)/2 + 1
IF (axic .AND. ktype == 1) midvec = 0
IF (axic .AND. ktype == 1) ivecn  = ivecn - icb(3)*ktype
IF (kwdest <= (buf7-iesta)) buf8 = buf7

!     OPEN ESTA

FILE = esta
CALL OPEN (*740,esta,z(buf8),rdrew)

!     REMAINING CORE

icore = buf8 - iesta
nesta = 0

!     WILL ESTA FIT IN CORE

IF (icore <= 0) CALL mesage (-8,0,nam)
IF (kwdest > icore) GO TO 140

!     ESTA WILL FIT. READ IT IN PLACING A ZERO WORD AT END OF EACH
!     RECORD.

i = iesta
110 CALL READ (*130,*120,esta,z(i),icore,1,nwords)
CALL REWIND (esta)
icore = buf8 - iesta
GO TO 140
120 i = i + nwords + 1
z(i-1) = 0
icore = icore - nwords - 1
GO TO 110

!     ALL ESTA NOW IN CORE

130 nesta = i - 1
CALL CLOSE (esta,clsrew)
IF (nesta > iesta) GO TO 140
WRITE  (opte,135) uwm
135 FORMAT (a25,' 3303, STRESSES OR FORCES REQUESTED FOR SET(S) ',  &
    'WHICH CONTAIN NO VALID ELEMENTS.')
GO TO 640

!     OPEN INPUT FILE. SKIP HEADER RECORD.

140 FILE  = ugv
CALL OPEN (*740,ugv,z(buf2),rdrew)
CALL fwdrec (*750,ugv)

!     IF ANY ISOPARAMETRIC ELEMENTS PRESENT, GET SECOND RECORD OF EQEXIN

IF (isopl == 0) GO TO 148
FILE = eqexin
CALL OPEN (*740,eqexin,z(buf7),rdrew)
CALL fwdrec (*750,eqexin)
CALL fwdrec (*750,eqexin)
isopl = eqexin
IF (isopl8 /= 8) GO TO 145
CALL fread (eqexin,z(ieqex),neqex,0)
CALL bckrec (eqexin)
145 CONTINUE

!     IF ANY STRESS OUTPUT IS REQUESTED,
!     OPEN OES1 AND WRITE HEADER RECORD

148 IF (stress == 0) GO TO 155
FILE = oes1
CALL OPEN (*151,oes1,z(buf3),wrtrew)
CALL fname (oes1,ocb)
DO  i = 1,3
  ocb(i+2) = date(i)
END DO
ocb(6) = time
ocb(7) = 1
CALL WRITE (oes1,ocb,7,1)
GO TO 155
151 CALL mesage (1,oes1,nam)
stress = 0

!     IF ANY STRESS OR FORCE OUTPUT IS REQUESTED AND COMPOSITE ELEMENTS
!     ARE PRESENT, OPEN OES1L AND OEF1L AND WRITE HEADER RECORDS

155 IF (comps /= -1 .OR. (stress == 0 .AND. force == 0)) GO TO 160
ilayer = 0
FILE   = oes1l
CALL OPEN (*158,oes1l,z(bufm1),wrtrew)
CALL WRITE (oes1l,nmes1l,2,1)
FILE = oef1l
CALL OPEN (*158,oef1l,z(buf0),wrtrew)
CALL WRITE (oef1l,nmef1l,2,1)
GO TO 160
158 CALL mesage (1,FILE,nam)
stress = 0
force  = 0

!     IF ANY FORCE OUTPUT IS REQUESTED,
!     OPEN OEF1 AND WRITE HEADER RECORD

160 IF (force == 0) GO TO 180
FILE = oef1
CALL OPEN (*171,oef1,z(buf4),wrtrew)
CALL fname (oef1,ocb)
DO  i = 1,3
  ocb(i+2) = date(i)
END DO
ocb(6) = time
ocb(7) = 1
CALL WRITE (oef1,ocb,7,1)
GO TO 180
171 CALL mesage (1,oef1,nam)
force = 0
180 IF (stress == 0 .AND. force == 0) GO TO 640

!     INITIALIZE UGV VEC, WHICH WILL BE THE NUMBER OF THE VECTOR WE
!     ARE NOW POSITIONED TO READ.

ugvvec = 1
isvvec = ivec
isvvcn = ivecn
iflag  = 0

!     READ A RECORD IN CASE CONTROL. SET SYMMETRY FLAG.

190 CALL READ (*610,*200,casecc,z(icc+1),kwdcc+1,1,flag)
CALL mesage (8,0,nam)
GO TO 640
200 ix  = icc + isymfl
symflg = z(ix)
ncc = icc + flag

!     FOR CONICAL SHELL SET SORC FLAG

ix = icc + isorc
IF (iflag  == 1) sorc   = isvsrc
IF (symflg == 0) sorc   = z(ix)
IF (sorc   == 1) axsine = .true.
IF (sorc   == 2) axcosi = .true.
IF (axic .AND. symflg == 0) isvsrc = sorc
ivec  = isvvec
ivecn = isvvcn
iflag = 0
IF (axic .AND. axsine .AND. axcosi .AND. ugvvec == 3) iflag = 1
IF (axic .AND. sorc == 0) GO TO 620

!     DETERMINE IF OUTPUT REQUEST IS PRESENT.
!     IF NOT, TEST FOR RECORD SKIP ON UGV THEN GO TO END OF THIS
!     REQUEST. IF SO, SET POINTERS TO SET DEFINING REQUEST.

210 ix = icc + istr
stresx = z(ix  )
sdest  = z(ix+1)
xsetns = -1
ix     = icc + ielf
forcex = z(ix  )
fdest  = z(ix+1)
xsetnf = -1
nstrop = z(icc+183)

!     DEBUG PRINTOUT


IF (comps == -1 .AND. nstrop > 1) ilayer = ilayer + 1
IF (stresx > 0.0) THEN
  GO TO   220
ELSE
  GO TO   240
END IF
220 ix = icc + ilsym
isetno = ix + z(ix) + 1
230 isets  = isetno + 2
nsets  = z(isetno+1) + isets - 1
IF (z(isetno) == stresx) GO TO 235
isetno = nsets + 1
IF (isetno <= ncc) GO TO 230
stresx = -1
GO TO 240

!     IF REQUIRED, LOCATE PRINT/PUNCH SUBSET FOR STRESSES

235 IF (stresx < xset0) GO TO 240
xsetns = sdest/10
sdest  = sdest - 10*xsetns
IF (xsetns == 0) GO TO 240
ixstns = ix + z(ix) + 1
236 ixsets = ixstns + 2
nxsets = z(ixstns+1) + ixsets - 1
IF (z(ixstns) == stresx) GO TO 240
ixstns = nxsets + 1
IF (ixstns < ncc) GO TO 236
stresx = -1
240 IF (forcex > 0.0) THEN
  GO TO   250
ELSE
  GO TO   270
END IF
250 ix = icc + ilsym
isetno = ix + z(ix) + 1
260 isetf  = isetno + 2
nsetf  = z(isetno+1) + isetf - 1
IF (z(isetno) == forcex) GO TO 265
isetno = nsetf + 1
IF (isetno <= ncc) GO TO 260
forcex = -1
GO TO 290

!     IF REQUIRED, LOCATE PRINT/PUNCH SUBSET FOR FORCES

265 IF (forcex < xset0) GO TO 290
xsetnf = fdest/10
fdest  = fdest - 10*xsetnf
IF (xsetnf == 0) GO TO 290
ixstnf = ix + z(ix) + 1
266 ixsetf = ixstnf + 2
nxsetf = z(ixstnf+1) + ixsetf - 1
IF (z(ixstnf) == forcex) GO TO 290
ixstnf = nxsetf + 1
IF (ixstnf < ncc) GO TO 266
forcex = -1
270 IF (stresx /= 0 .OR. forcex /= 0 .OR. axic) GO TO 290

!     NO REQUESTS THIS CC RECORD FOR STRESSES OR FORCES.
!     THUS SKIP CORRESPONDING UGV RECORD UNLESS SYMFLG IS ON, IN WHICH
!     CASE WE SKIP NO UGV RECORD SINCE THE SYMMETRY CASE HAS NO UGV
!     VECTOR, BUT IN FACT WOULD HAVE USED A SUMMATION OF THE IMMEDIATELY
!     PRECEEDING LSYM VECTORS.

!     IF END OF CC AND NO STRESS OR FORCE OUTPUT REQUEST WE ARE DONE

IF (eofcc ) GO TO 620
IF (symflg == 0.0) THEN
  GO TO   280
ELSE
  GO TO   190
END IF
280 CALL fwdrec (*750,ugv)
ugv vec = ugv vec + 1
GO TO 570

!     THERE IS A REQUEST FOR STRESSES AND OR FORCES
!     FIRST DETERMINE APPROPRIATE GPTT AND EDT RECORDS IF REQUIRED

290 ix     = icc + itload
tloads = z(ix)
ngptt  = 0
IF (tloads == 0) GO TO 370
FILE   = gptt
CALL CLOSE (gptt,clsrew)
CALL OPEN (*740,gptt,z(buf5),rdrew)

!     SKIP NAME

CALL READ (*750,*751,gptt,buf,2,0,n)

!     PICK UP 3 WORDS OF SET INFORMATION

295 CALL READ (*750,*751,gptt,buf,3,0,n)
IF (buf(1) /= tloads) GO TO 295
deftmp = bufr(2)
tmprec = buf(3)

370 ix    = icc + ieldef
eldef = z(ix)
IF (eldef == 0 .OR. kwdedt == 0) GO TO 430
FILE  = edt
CALL preloc (*740,z(buf6),edt)
CALL locate (*390,z(buf6),kdefrm,flag)
idef  = iedt
i     = idef
380 CALL READ (*750,*390,edt,buf(1),3,0,flag)
IF (buf(1) == eldef) GO TO 410
GO TO 380
390 buf(1) = eldef
buf(2) = 0
CALL mesage (-30,46,buf)
400 CALL READ (*750,*420,edt,buf(1),3,0,flag)
IF (buf(1) /= eldef) GO TO 420
410 z(i  ) = buf(2)
z(i+1) = buf(3)
i = i + 2
IF (i < igptta) GO TO 400
CALL mesage (-8,0,nam)
420 ndef = i - 2
CALL CLOSE (edt,clsrew)

!     UNPACK VECTOR INTO CORE

430 coef1 = 1.0
IF (symflg == 0) GO TO 490

!     SYMMETRY SEQUENCE-- BUILD VECTOR IN CORE.

ix   = icc + ilsym
lsym = z(ix)

!     IF SYMFLG IS NEGATIVE, THIS IS A REPEAT SUBCASE.  USE PRESENT
!     VECTOR IN CORE.

IF (symflg < 0 .AND. app(1) == sta(1)) GO TO 530
IF (symflg < 0) GO TO 190
DO  i = ivec,ivecn
  zz(i) = 0.0
END DO
IF (lsym > ugv vec-1) GO TO 780
limit = lsym
IF (iflag == 1) limit = 1
DO  i = 1,limit
  CALL bckrec (ugv)
END DO
isymn = ix + lsym
i = ix + 1
IF (iflag == 1) i = i + 1
j2 = icb(3)
460 coef = zz(i)
CALL intpk (*480,ugv,0,qtype2,0)
470 CALL zntpki
ix = ivec + ixx - 1
IF (ktype == 1) GO TO 471
zz(ix+j2) = zz(ix+j2) + coef*xx(1)
zz(ix   ) = zz(ix)    + coef*xx(2)
GO TO 472
471 CONTINUE
zz(ix)= zz(ix) + coef*xx(1)
472 CONTINUE
IF (eol   == 0) GO TO 470
480 IF (iflag == 1) GO TO 485
i = i + 1
IF (i <= isymn) GO TO 460
GO TO 530

!     CONICAL SHELL BOTH CASE
!     2 VECTORS IN CORE -
!     2-ND VECTOR IS NOW IN CORE AT Z(IVEC) THRU Z(IVECN)...
!     GET 1-ST VECTOR AND PUT IT AT Z(IVECN+1) THRU Z(2*IVECN-MIDVEC+1)


485 midvec = ivec
ivec   = ivecn + 1
ivecn  = ivecn + (ivecn-midvec+1)
coef1  = zz(icc + ilsym+1)

!     IF FALL HERE AND SORC=1 THE VECTOR IN CORE IS THE SINE VECTOR AND
!     IF SORC=2 THE VECTOR IN CORE IS THE COSINE VECTOR.  THUS THE FIRST
!     VECTOR WAS THE OTHER VECTOR RESPECTIVELY
!     BY THE WAY THE VECTOR IN CORE IS THE SECOND VECTOR.

CALL bckrec (ugv)
CALL bckrec (ugv)

!     NOT SYMMETRY-- UNPACK VECTOR.

490 j2 = icb(3)
IF (iflag == 1) GO TO 515
IF (ugvvec > nvects) GO TO 620
515 DO  i = ivec,ivecn
  zz(i) = 0.0
END DO
CALL intpk (*500,ugv,0,qtype2,0)
491 CALL zntpki
ix = ivec + ixx-1
IF (ktype == 1) GO TO 492
zz(ix   ) = coef1*xx(2)
zz(ix+j2) = coef1*xx(1)
GO TO 493
492 CONTINUE
zz(ix) = coef1*xx(1)
493 CONTINUE
IF (eol == 0) GO TO 491
495 IF (app(1) /= trn(1)) GO TO 520
CALL fwdrec (*520,ugv)
ugvvec = ugvvec + 1
CALL fwdrec (*520,ugv)
ugvvec = ugvvec + 1
GO TO 520
500 CONTINUE
GO TO 495
520 IF (iflag /= 1) ugvvec = ugvvec + 1
IF (iflag == 1) CALL skprec (ugv,1)

!     READY NOW TO SWEEP THROUGH THE ESTA ONCE.
!     SDR2E DOES ALL THE PROCESSING OF PHASE II ELEMENT COMPUTATIONS.
!     THE ESTA FILE, BE IT IN CORE OR NOT, IS SWEPT THRU ONCE FOR THE
!     FOLLOWING CALL.

530 IF (iflag == 1) sorc = sorc + 1
IF (sorc  == 3) sorc = 1
CALL sdr2e (*640,ieqex,neqex)

!     CONCLUDE PROCESSING OF THIS VECTOR
!     INITIALIZE FOR NEXT VECTOR
!     CANCEL THIS INITIALIZATION IN SOME CASES IF A REPEAT CASE.

570 SELECT CASE ( branch )
  CASE (    1)
    GO TO 580
  CASE (    2)
    GO TO 581
  CASE (    3)
    GO TO 620
  CASE (    4)
    GO TO 581
  CASE (    5)
    GO TO 590
  CASE (    6)
    GO TO 582
  CASE (    7)
    GO TO 620
  CASE (    8)
    GO TO 581
  CASE (    9)
    GO TO 581
  CASE (   10)
    GO TO 580
END SELECT

580 IF (.NOT.eofcc) GO TO 190
GO TO 589
581 jlist = jlist + ktype1
IF (.NOT.eofcc) GO TO 190
GO TO 589

!     TRANSIENT RESPONSE

582 jlist = jlist + 2
IF (jlist <= nlist .AND. .NOT.eofcc) GO TO 190
IF (jlist > nlist  .OR. ugvvec > nvects) GO TO 620
GO TO 490

!     PROCESS ANY REMAINING VECTORS WITH LAST CC RECORD

589 IF (ugv vec <= nvects .AND. symflg == 0) GO TO 210
GO TO 620

!     FREQUENCY RESPONSE, PICK UP NEXT VECTOR UNLESS ALL FREQUENCIES
!     COMPLETED

590 jlist = jlist + 2
IF (jlist <= nlist .AND. ugv vec <= nvects) GO TO 210
kfrq  = 0
jlist = ilist
DO  i = ilist,nlist,2
  z(i+1) = 0
END DO
IF (ugv vec <= nvects) GO TO 190
GO TO 620

!     EOF HIT ON CASECC FILE
!     PROCESS ANY MORE VECTORS USING LAST CASECC RECORD

610 eofcc = .true.
IF (nvects >= ugv vec) GO TO 210

!     WRITE TRAILERS AND CLOSE ANY OPEN FILES

620 ocb(2) = 63
IF (stress == 0) GO TO 630
ocb(1) = oes1
CALL wrttrl (ocb(1))
IF (comps /= -1 .OR. ilayer == 0) GO TO 630
ocb(1) = oes1l
CALL wrttrl (ocb(1))
630 IF (force == 0) GO TO 640
ocb(1) = oef1
CALL wrttrl (ocb(1))
IF (comps /= -1 .OR. ilayer == 0) GO TO 640
ocb(1) = oef1l
CALL wrttrl (ocb(1))
640 DO  i = 1,12
  SELECT CASE ( i )
    CASE (    1)
      GO TO 650
    CASE (    2)
      GO TO 660
    CASE (    3)
      GO TO 670
    CASE (    4)
      GO TO 680
    CASE (    5)
      GO TO 690
    CASE (    6)
      GO TO 700
    CASE (    7)
      GO TO 710
    CASE (    8)
      GO TO 720
    CASE (    9)
      GO TO 721
    CASE (   10)
      GO TO 725
    CASE (   11)
      GO TO 726
    CASE (   12)
      GO TO 728
  END SELECT
  650 FILE = oes1
  GO TO 730
  660 FILE = oef1
  GO TO 730
  670 FILE = ugv
  GO TO 730
  680 FILE = casecc
  GO TO 730
  690 FILE = edt
  GO TO 730
  700 FILE = gptt
  GO TO 730
  710 FILE = pg
  GO TO 730
  720 FILE = eigr
  GO TO 730
  721 FILE = esta
  GO TO 730
  725 FILE=eqexin
  GO TO 730
  726 FILE = oes1l
  GO TO 730
  728 FILE = oef1l
  730 CALL CLOSE (FILE,clsrew)
END DO
RETURN

740 n = 1
GO TO 760
750 n = 2
GO TO 760
751 n = 3
GO TO 760
760 CALL mesage (n,FILE,nam)
GO TO 640

!     UGV FILE PURGED, CAN NOT PROCESS STRESSES OR FORCES

770 CALL mesage (30,76,0)
GO TO 640
780 ocb(1) = lsym
ocb(2) = ugv vec - 1
CALL mesage (30,92,ocb(1))
GO TO 620
END SUBROUTINE sdr2d
