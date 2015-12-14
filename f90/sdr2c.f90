SUBROUTINE sdr2c
     
!     SDR2C PROCESSES OUTPUT REQUESTS FOR SINGLE-POINT FORCES OF
!     CONSTRAINT, LOADS, DISPLACEMENTS, VELOCITIES, ACCELERATIONS AND
!     EIGENVECTORS.
 
 LOGICAL :: anyout,axic  ,ddrmm ,axsine,axcosi
 INTEGER :: app   ,sort2 ,spcf  ,displ ,vel   ,acc   ,stress,  &
     force ,cstm  ,casecc,eqexin,sil   ,bgpdt ,pg    ,  &
     qg    ,ugv   ,phig  ,eigr  ,opg1  ,oqg1  ,ougv1 ,  &
     pugv1 ,ocb   ,sorc  ,dtype ,FILE  ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,buf5  ,symflg,outfl ,sta   ,rei   ,  &
     ds0   ,ds1   ,frq   ,trn   ,bk0   ,date  ,sysbuf,  &
     branch,plots ,qtype2,eol   ,bk1   ,time  ,setno ,  &
     fsetno,z     ,retx  ,formt ,flag  ,eof   ,cei   ,  &
     pla   ,oharms,blanks,harms ,xsetno,xset0 ,dest  ,  &
     pbuff(4)     ,extra ,axif  ,edt   ,platit(12)   , buf(50)
 REAL :: zz(1) ,bufr(11)     ,pbufr(4)
 DIMENSION       date(3)
 COMMON /BLANK / app(2),sort2
 COMMON /sdr2x1/ ieigen,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym ,ifrout, isload,idload,isorc
 COMMON /sdr2x2/ casecc,cstm  ,mpt   ,dit   ,eqexin,sil   ,gptt  ,  &
     edt   ,bgpdt ,pg    ,qg    ,ugv   ,est   ,phig  ,  &
     eigr  ,opg1  ,oqg1  ,ougv1 ,oes1  ,oef1  ,pugv1 ,  &
     oeigr ,ophig ,pphig ,esta  ,gptta ,harms
COMMON /sdr2x4/ nam(2),END   ,mset  ,icb(7),ocb(7),mcb(7),dtype(8)  &
    ,               icstm ,ncstm ,ivec  ,ivecn ,temp  ,deform,FILE  ,  &
    buf1  ,buf2  ,buf3  ,buf4  ,buf5  ,any   ,all   ,  &
    tloads,eldef ,symflg,branch,ktype ,loads ,spcf  ,  &
    displ ,vel   ,acc   ,stress,force ,kwdest,kwdedt,  &
    kwdgpt,kwdcc ,nrigds,sta(2),rei(2),ds0(2),ds1(2),  &
    frq(2),trn(2),bk0(2),bk1(2),cei(2),pla(22)      ,  &
    nrings,nharms,axic  ,knset ,isopl ,strspt,ddrmm
COMMON /zzzzzz/ z(1)
COMMON /condas/ pi    ,twopi ,raddeg,degra ,s4pisq
COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
COMMON /system/ ksystm(65)
COMMON /unpakx/ qtype2,i2    ,j2    ,incr2
COMMON /zntpkx/ xx(4),ixx    ,eol   ,eor
COMMON /zblpkx/ y(4) ,iy
EQUIVALENCE     (ksystm( 1),sysbuf) ,(ksystm(15),date(1)) ,  &
    (ksystm(18),time  ) ,(ksystm(20),plots  ) ,  &
    (ksystm(38),axif  ) ,(ksystm(56),iheat  ) ,  &
    (buf(1),bufr(1)),(z(1),zz(1)),(pbuff(1),pbufr(1))
DATA    buf   / 50*0     /
DATA    blanks/ 4H       /
DATA    xset0 / 100000000/
DATA    platit/ 4HLOAD,4H fac,4HTOR ,9*0/
DATA    mmreig/ 4HMMRE   /

!     IF THIS IS A DYNAMIC-DATA-RECOVERY-MATRIX-METHOD REIG PROBLEM
!     THEN ALL EIGENVECTORS ARE TO BE OUTPUT FOR THE DDRMM MODULE.

setno = 0
IF (ddrmm .AND. ireq == idispl) setno = -1

!     PERFORM GENERAL INITIALIZATION

buf1  = korsz(z) - sysbuf + 1
buf2  = buf1 - sysbuf
buf3  = buf2 - sysbuf
buf4  = buf3 - sysbuf
buf5  = buf4 - sysbuf
iseq  = 1
m8    =-8
i2    = 1
incr2 = 1
kplot = 0
extra = 0
axsine = .false.
axcosi = .false.

!     READ SECOND RECORD OF EQEXIN OR EQDYN INTO CORE.

FILE = eqexin
CALL gopen (eqexin,z(buf1),0)
CALL skprec (eqexin,1)
CALL READ (*1320,*30,eqexin,z,buf5,1,neqex)
CALL mesage (m8,0,nam)
30 CALL CLOSE (eqexin,clsrew)
itabl= 1
kn   = neqex/2
icc  = neqex
ilist= neqex + 1

!     INITIALIZE FOR PROCESSING SPECIFIC REQUEST.

40 IF (iseq-2 < 0) THEN
  GO TO    50
ELSE IF (iseq-2 == 0) THEN
  GO TO    60
ELSE
  GO TO    70
END IF

!     LOAD VECTOR.

50 IF (loads == 0 .OR. app(1) == rei(1) .OR. app(1) == cei(1) .OR.  &
    app(1) == bk1(1)) GO TO 1180
infil = 115
outfl = opg1
ireq  = iloads
GO TO 90

!     SINGLE-POINT FORCES OF CONSTRAINT.

60 IF (spcf == 0) GO TO 1180
infil = qg
outfl = oqg1
ireq  = ispcf
GO TO 90

!     DISPLACEMENT VECTOR OR EIGENVECTOR

70 IF (displ /= 0 .OR. vel /= 0 .OR. acc /= 0 .OR. plots /= 0) GO TO 80
GO TO 1180
80 infil = ugv
outfl = ougv1
jtj   = vel + acc
IF (.NOT.(app(1) == mmreig .AND. displ == 0 .AND. jtj /= 0)) GO TO 88
IF (vel == 0) GO TO 84
ireq = ivel
GO TO 90
84 ireq = iacc
GO TO 90
88 ireq = idispl

!     READ TRAILER ON INPUT FILE. SET PARAMETERS.

90 icb(1) = infil
CALL rdtrl (icb)
IF (icb(1) /= infil) GO TO 1200
nvects = icb(2)
IF (icb(5) > 2) GO TO 100

!     REAL VECTOR.

ktype  = 1
qtype2 = 1
ktype1 = 2
nwds   = 8
ktypex = 0
GO TO 110

!     COMPLEX VECTOR.

100 ktype  = 2
qtype2 = 3
ktype1 = 3
nwds   = 14
ktypex = 1000

!     OPEN CASE CONTROL AND SKIP HEADER. THEN BRANCH ON APPROACH.

110 CALL gopen (casecc,z(buf1),0)
pbuff(2) = 1
SELECT CASE ( branch )
  CASE (    1)
    GO TO 190
  CASE (    2)
    GO TO 120
  CASE (    3)
    GO TO 190
  CASE (    4)
    GO TO 150
  CASE (    5)
    GO TO 160
  CASE (    6)
    GO TO 160
  CASE (    7)
    GO TO 190
  CASE (    8)
    GO TO 150
  CASE (    9)
    GO TO 120
  CASE (   10)
    GO TO 190
END SELECT

!     EIGENVALUES - READ LIST OF MODE NOS. AND EIGENVALUES INTO CORE.

120 FILE = eigr
CALL gopen (eigr,z(buf2),0)
CALL skprec (eigr,1)
IF (app(1) == cei(1)) pbuff(2) = 5
IF (app(1) == rei(1)) pbuff(2) = 4
i = ilist
m = 8 - ktype
iskip = 0
INDEX = 2
IF (app(1) /= rei(1)) GO TO 130

!     CHECK TO SEE IF ALL GENERALIZED MASS VALUES ARE ZERO

125 CALL READ (*1320,*127,eigr,buf,m,0,flag)
IF (buf(6) == 0.0) GO TO 125
INDEX = 0
127 CALL skprec (eigr,-1)
130 CALL READ (*1320,*140,eigr,buf,m,0,flag)
IF (app(1) /= rei(1)) GO TO 135
IF (INDEX == 2) GO TO 135

!     MATCH CORRECT MODE NOS. AND EIGENVALUES WITH PROPER
!     EIGENVECTORS WHEN USING GIVENS METHOD WITH F1.GT.0.0

IF (INDEX  ==  1) GO TO 135
IF (buf(6) /= 0.) GO TO 133
iskip = iskip + 1
GO TO 130
133 INDEX = 1
135 z(i  ) = buf(1) - iskip
z(i+1) = buf(3)
z(i+2) = buf(4)
i = i + ktype1
GO TO 130
140 CALL CLOSE (eigr,clsrew)
nlist = i - ktype1
icc   = i
GO TO 190

!     DIFF. STIFF. PHASE 1 OR BUCKLING PHASE 1 - SKIP 1ST DATA RECORD ON
!     CC.

150 CALL skprec (casecc,1)
pbuff(2) = 4
IF (app(1) == bk1(1)) GO TO 120
pbuff(2) = 1
GO TO 190

!     FREQUENCY OR TRANSIENT RESPONSE - READ LIST INTO CORE.

160 FILE = pg
CALL OPEN (*1310,FILE,z(buf2),rdrew)
i  = ilist
m  = 3
ix = 1
pbuff(2) = 3
IF (app(1) == frq(1)) pbuff(2) = 2
IF (app(1) == frq(1) .OR. app(1) == trn(1)) ix = 2
170 CALL READ (*1320,*180,FILE,buf,m,0,flag)
z(i  ) = buf(m)
z(i+1) = 0
i = i + ix
m = 1
GO TO 170
180 CALL CLOSE (FILE,clsrew)
nlist = i - ix
icc   = i

!     OPEN OUTPUT FILE. WRITE HEADER RECORD.

190 FILE   = outfl
anyout = .false.
CALL OPEN (*1200,outfl,z(buf2),wrtrew)
ocb(1) = outfl
CALL fname (outfl,buf)
DO  i = 1,3
  buf(i+2) =  date(i)
END DO
buf(6) = time
buf(7) = 1
CALL WRITE (outfl,buf,7,1)

!     OPEN INPUT FILE. SKIP HEADER RECORD.

FILE   = infil
CALL OPEN (*1190,infil,z(buf3),rdrew)
CALL fwdrec (*1320,infil)

!     SET PARAMETERS TO KEEP CASE CONTROL AND VECTORS IN SYNCH.

eof    = 0
jcount = 0
kcount = 1
jlist  = ilist
kfrq   = 0
incore = 0
kwds   = 0

!     READ A RECORD IN CASE CONTROL. SET SYMMETRY FLAG.

230 CALL READ (*1160,*240,casecc,z(icc+1),buf5-icc,1,ncc)
CALL mesage (m8,0,nam)
240 ix     = icc + isymfl
itemp  = icc + harms

!     OHARMS WILL BE 1 GREATER THAN THE MAXIMUM OUTPUT HARMONIC

oharms = z(itemp)
IF (oharms < 0 .AND. axif /= 0) oharms = axif
IF (oharms < 0) oharms = nharms

!     IF A FLUID PROBLEM CONVERT USER HARMONIC TO INTERNAL HARMONIC MAX.

IF (oharms == 0) GO TO 243
IF (axif   == 0) GO TO 243
oharms = oharms - 1
oharms = 2*oharms + 3
243 symflg = z(ix)
IF (symflg == 0) sorc = z(icc+isorc)
IF (sorc == 1) axsine = .true.
IF (sorc == 2) axcosi = .true.
iflag  = 0
IF (axic .AND. axsine .AND. axcosi .AND. jcount == 2) iflag = 1
ivec   = icc + ncc + 1

!     DETERMINE IF OUTPUT REQUEST IS PRESENT.
!     IF NOT, TEST FOR RECORD SKIP ON INFIL  THEN GO TO END OF THIS
!     REQUEST.
!     IF SO, SET POINTERS TO SET DEFINING REQUEST.

250 ireqx  = icc + ireq
setno  = z(ireqx  )
dest   = z(ireqx+1)
formt  = IABS(z(ireqx+2))
xsetno = -1
IF (setno < 0.0) THEN
  GO TO   300
ELSE IF (setno == 0.0) THEN
  GO TO   260
ELSE
  GO TO   280
END IF
260 IF (symflg /= 0) GO TO 1000
IF (app(1) /= frq(1)) GO TO 270
IF (iseq  == 3) GO TO 300
270 IF (plots /= 0) GO TO 300
CALL fwdrec (*1320,infil)
jcount = jcount + 1
GO TO 1000
280 ix     = icc + ilsym
isetno = ix + z(ix) + 1
290 iset   = isetno + 2
nset   = z(isetno+1) + iset - 1
IF (z(isetno) == setno) GO TO 295
isetno = nset + 1
IF (isetno < ivec) GO TO 290
setno  = -1
GO TO 300

!     IF REQUIRED, LOCATE PRINT/PUNCH SUBSET.

295 IF (setno < xset0) GO TO 300
xsetno = dest/10
dest   = dest - 10*xsetno
IF (xsetno == 0) GO TO 300
ixsetn = ix + z(ix) + 1
296 ixset  = ixsetn + 2
nxset  = z(ixsetn+1) + ixset - 1
IF (z(ixsetn) == xsetno) GO TO 300
ixsetn = nxset + 1
IF (ixsetn < ivec) GO TO 296
xsetno = -1

!     UNPACK VECTOR INTO CORE (UNLESS VECTOR IS ALREADY IN CORE).

300 IF (incore /= 0) GO TO 400
ivecn = ivec + ktype*icb(3) - 1
IF (ivecn  >= buf5) CALL mesage (m8,0,nam)
IF (symflg == 0) GO TO 360

!     SYMMETRY SEQUENCE - BUILD VECTOR IN CORE.

ix   = icc + ilsym
lsym = z(ix)

!     IF SYMFLG IS NEGATIVE THIS IS A REPEAT SUBCASE. BCKREC VECTOR
!     AND READ IT INTO CORE.

IF (symflg < 0 .AND. app(1) == sta(1)) GO TO 358
IF (symflg < 0) GO TO 230
DO  i = ivec,ivecn
  zz(i) = 0.
END DO
DO  i = 1,lsym
  CALL bckrec (infil)
END DO
isymn = ix + lsym
i = ix + 1
330 coef = zz(i)
CALL intpk (*350,infil,0,qtype2,0)
340 CALL zntpki
ix = ivec + ixx - 1
zz(ix) = zz(ix) + coef*xx(1)
IF (ktype == 2) zz(ix+1) = zz(ix+1) + coef*xx(2)
IF (eol == 0) GO TO 340
350 i = i + 1
IF (i <= isymn) GO TO 330
GO TO 400

!     REPEAT SUBCASE

358 jcount = jcount - 1
CALL bckrec (infil)

!     NOT SYMMETRY - UNPACK VECTOR.

360 j2= icb(3)
IF (jcount >= nvects) GO TO 1170
CALL unpack (*370,infil,z(ivec))
GO TO 390
370 DO  i = ivec,ivecn
  zz(i)  = 0.
END DO
390 jcount = jcount + 1

!     TEST FOR CONTINUATION FROM HERE.

400 IF (setno /= 0) GO TO 410
IF (app(1) == frq(1)) GO TO 1040

!     PREPARE TO WRITE ID RECORD ON OUTPUT FILE.

410 SELECT CASE ( branch )
  CASE (    1)
    GO TO 420
  CASE (    2)
    GO TO 430
  CASE (    3)
    GO TO 420
  CASE (    4)
    GO TO 420
  CASE (    5)
    GO TO 440
  CASE (    6)
    GO TO 560
  CASE (    7)
    GO TO 420
  CASE (    8)
    GO TO 430
  CASE (    9)
    GO TO 430
  CASE (   10)
    GO TO 420
END SELECT

!     NORMAL STATICS OR DIFF.STIFF. PHASE O OR 1 OR BUCKLING PHASE 0.

420 buf(2) = dtype(iseq)
ix = icc + isload
buf(5) = z(icc+1)
buf(6) = 0
buf(7) = 0
buf(8) = z(ix)
pbuff(2) = 1
pbuff(3) = z(ix)
pbuff(4) = 0
IF (branch /= 10) GO TO 610
ix = icc + ittl + 84
z(ix  ) = platit(1)
z(ix+1) = platit(2)
z(ix+2) = platit(3)
CALL int2al (jcount,z(ix+3),platit(4))
GO TO 610

!     EIGENVALUES OR BUCKLING PHASE 1.

430 IF (iseq == 2) buf(2) = ktypex + 3
IF (iseq == 3) buf(2) = ktypex + 7
buf(5) = z(jlist  )
buf(6) = z(jlist+1)
buf(7) = z(jlist+2)
buf(8) = 0
!     PBUFF(2) = 2  THIS CARD WAS REMOVED SINCE LEVEL 16. NO LONGER NEED
pbuff(3) = buf(5)
IF (app(1) == bk1(1)) pbuff(3) = -buf(5)
pbuff(4) = buf(6)
IF (app(1) /= bk1(1) .AND. app(1) /= cei(1))  &
    pbufr(4) = SQRT(ABS(bufr(6)))/twopi
IF (app(1) == cei(1)) pbufr(4) = ABS(bufr(7))/twopi
GO TO 610

!     FREQUENCY RESPONSE.

440 ix = icc + idload
buf(8) = z(ix)
buf(6) = 0
buf(7) = 0
pbuff(2) = 2
pbuff(3) = buf(8)
IF (iseq == 3) GO TO 520
buf(2) = dtype(iseq) + ktypex
GO TO 441
520 IF (kcount-2 < 0) THEN
  GO TO   530
ELSE IF (kcount-2 == 0) THEN
  GO TO   540
ELSE
  GO TO   550
END IF
530 buf(2) = 1001
GO TO 441
540 buf(2) = 1010
GO TO 441
550 buf(2) = 1011
GO TO 441
441 CONTINUE
IF (kfrq /= 0) GO TO 510

!     FIRST TIME FOR THIS LOAD VECTOR ONLY - MATCH LIST OF USER
!     REQUESTED FREQS WITH ACTUAL FREQS. MARK FOR OUTPUT EACH ACTUAL
!     FREQ WHICH IS CLOSEST TO USER REQUEST.

kfrq   = 1
ix     = icc + ifrout
fsetno = z(ix)
IF (fsetno <= 0) GO TO 460
ix     = icc + ilsym
isetnf = ix  + z(ix) + 1
450 isetf  = isetnf + 2
nsetf  = z(isetnf+1) + isetf - 1
IF(z(isetnf) == fsetno) GO TO 480
isetnf = nsetf + 1
IF (isetnf < ivec) GO TO 450
fsetno = -1
460 DO  j = ilist,nlist,2
  z(j+1) = 1
END DO
GO TO 510
480 DO  i = isetf,nsetf
  k    = 0
  diff = 1.e25
  bufr(1) = zz(i)
  DO  j = ilist,nlist,2
    IF (z(j+1) /= 0) CYCLE
    diff1 = ABS(zz(j)-bufr(1))
    IF (diff1 >= diff) CYCLE
    diff = diff1
    k = j
  END DO
  IF (k /= 0) z(k+1) = 1
END DO

!     DETERMINE IF CURRENT FREQ IS MARKED FOR OUTPUT.

510 IF (z(jlist+1) == 0) GO TO 1000
buf(5)   = z(jlist)
pbuff(4) = buf(5)
GO TO 610

!     TRANSIENT RESPONSE.

560 buf(5) = z(jlist)
IF (kcount - 2 < 0) THEN
  GO TO   570
ELSE IF (kcount - 2 == 0) THEN
  GO TO   580
ELSE
  GO TO   590
END IF
570 buf(2) = 1
GO TO 600
580 buf(2) = 10
GO TO 600
590 buf(2) = 11
600 IF (ireq == iloads) buf(2) = 2
IF (ireq == ispcf ) buf(2) = 3
ix = icc + idload
buf(8) = z(ix)
buf(6) = 0
buf(7) = 0
pbuff(2) = 3 + 10*(kcount-1)
pbuff(3) = buf(8)
pbuff(4) = buf(5)
GO TO 441

!     WRITE ID RECORD ON OUTPUT FILE.

610 IF (setno == 0 .AND. plots /= 0) GO TO 880
buf(1) = dest + 10*branch
buf(3) = 0

!     IF CONICAL SHELL PROBLEM, SET MINOR ID = 1000 FOR USE BY OFP

IF (axic) buf(3) = 1000
buf(4) = z(icc+1)
IF (ddrmm) buf(4) = 9999
buf(9) = IABS(z(ireqx+2))
IF (buf(9) == 1 .AND. ktype == 2) buf(9) = 2
formt  = buf(9)
buf(10)= nwds
CALL WRITE (outfl,buf,50,0)
ix = icc + ittl
CALL WRITE (outfl,z(ix),96,1)

!     BUILD DATA RECORD ON OUTPUT FILE.

IF (setno /= -1) GO TO 650

!     SET .EQ. ALL  -  OUTPUT ALL POINTS DEFINED IN EQEXIN.

kx = 1
n  = neqex - 1
ASSIGN 640 TO retx
GO TO 700
640 kx = kx + 2
IF (kx <= n) GO TO 700
GO TO 880

!     SET .NE. ALL  -  OUTPUT ONLY POINTS DEFINED IN SET.

650 jharm = 0
651 i = iset
660 IF (i   == nset) GO TO 680
IF (z(i+1) > 0) GO TO 680
n = -z(i+1)
buf(1) = z(i)
ibufsv = buf(1)
i = i + 1
ASSIGN 670 TO retx
GO TO 1210
670 buf(1) = ibufsv + 1
ibufsv = buf(1)
IF (buf(1) <= n) GO TO 1210
GO TO 690
680 buf(1) = z(i)
ASSIGN 690 TO retx
GO TO 1210
690 i = i + 1
IF (i <= nset) GO TO 660
jharm = jharm + 1
IF (.NOT.axic .AND. axif == 0) GO TO 880
IF (jharm <= oharms) GO TO 651
GO TO 880

!     PICK UP POINTER TO GRID POINT DATA AND GRID POINT TYPE.

700 buf(1) = z(kx)
IF (iflag == 1 .AND. buf(1) >= 1000000) GO TO retx, (640,670,690)
j = z(kx+1)/10
buf(2) = z(kx+1) - 10*j
j = ivec + ktype*(j-1)
IF (buf(2) == 1) GO TO 770

!     SCALAR OR EXTRA POINT.

buf(3) = z(j)
IF (ktype == 2) GO TO 720
IF (iseq <= 2 .AND. bufr(3) == 0.0 .AND. sort2 < 0) GO TO retx, (640,670,690)
DO  k = 4,8
  buf(k) = 0
END DO
GO TO 840

!     COMPLEX SCALAR OR EXTRA POINT.

720 buf(4) = z(j+1)
IF (iseq <= 2 .AND. bufr(3) == 0.0 .AND. bufr(4) == 0.0 .AND.  &
    sort2 < 0) GO TO retx, (640,670,690)
DO  k = 5,14
  buf(k) = 0
END DO
IF (formt /= 3) GO TO 840
redner = SQRT(bufr(3)**2 + bufr(4)**2)
IF (redner == 0.0) THEN
  GO TO   740
ELSE
  GO TO   750
END IF
740 bufr(4) = 0.0
GO TO 760
750 bufr(4) = ATAN2(bufr(4),bufr(3))*raddeg
IF (bufr(4) < -0.00005) bufr(4) = bufr(4) + 360.0
760 bufr(3) = redner
GO TO 840

!     GRID POINT.

770 flag = 0
IF (ktype == 2) GO TO 790
DO  k = 1,6
  bufr(k+2) = zz(j)
  IF (bufr(k+2) /= 0.0 .OR. sort2 >= 0) flag = 1
  j = j + 1
END DO
IF (iseq <= 2 .AND. flag == 0) GO TO retx, (640,670,690)
GO TO 840

!     COMPLEX GRID POINT.

790 DO  k = 1,11,2
  bufr(k+2) = zz(j  )
  bufr(k+3) = zz(j+1)
  IF (bufr(k+2) /= 0. .OR. bufr(k+3) /= 0. .OR. sort2 >= 0) flag = 1
  IF (formt /= 3) GO TO 830
  redner = SQRT(bufr(k+2)**2 + bufr(k+3)**2)
  IF (redner == 0.0) THEN
    GO TO   800
  ELSE
    GO TO   810
  END IF
  800 bufr(k+3) = 0.0
  GO TO 820
  810 bufr(k+3) = ATAN2( bufr(k+3),bufr(k+2) )*raddeg
  IF (bufr(k+3) < -0.00005) bufr(k+3) = bufr(k+3) + 360.0
  820 bufr(k+2) = redner
  830 j = j + 2
END DO
IF (iseq <= 2 .AND. flag == 0) GO TO retx, (640,670,690)

!     WRITE ENTRY ON OUTPUT FILE.

!     IF COMPLEX  TRANSPOSE DATA FOR OFP (REAL TOP, IMAG BOTTOM)

840 IF (nwds /= 14) GO TO 850
itemp   = buf( 4)
buf( 4) = buf( 5)
buf( 5) = buf( 7)
buf( 7) = buf(11)
buf(11) = buf( 8)
buf( 8) = buf(13)
buf(13) = buf(12)
buf(12) = buf(10)
buf(10) = buf( 6)
buf( 6) = buf( 9)
buf( 9) = itemp

850 anyout = .true.

!     IF CONICAL SHELL DECODE GRID POINT NUMBER IF GREATER THAN 1000000.

IF (.NOT.axic) GO TO 870
IF (buf(1) >= 1000000) GO TO 860
buf(2) = blanks
GO TO 870
860 itemp = buf(1)/1000000

!     STOP OUTPUT WHEN PRESENT HARMONIC EXCEEDS OUTPUT HARMONIC SIZE REQ

IF (itemp > oharms) GO TO 880
buf(1) = buf(1) - itemp*1000000
buf(2) = itemp - 1
GO TO 876

!     IF A FLUID PROBLEM THEN A CHECK IS MADE ON THE HARMONIC ID

870 IF (axif > 0.0) THEN
  GO TO   861
ELSE
  GO TO   876
END IF
861 IF (buf(1) < 500000) GO TO 876
itemp = buf(1) - MOD(buf(1),500000)
itemp = itemp/500000

!     STOP THE OUTPUT IF THE HARMONIC IS GREATER THAN THE OUTPUT
!     REQUEST FOR HARMONICS

IF (itemp >= oharms) GO TO 880

!     DETERMINE DESTINATION FOR ENTRY.

876 id = buf(1)
buf(1) = 10*id + dest
IF (xsetno < 0.0) THEN
  GO TO   878
ELSE IF (xsetno == 0.0) THEN
  GO TO   871
ELSE
  GO TO   872
END IF
871 buf(1) = 10*id
GO TO 878
872 ix = ixset
873 IF (ix ==  nxset) GO TO 874
IF (z(ix+1) > 0) GO TO 874
IF (id >= z(ix) .AND. id <= -z(ix+1)) GO TO 878
ix = ix + 2
GO TO 875
874 IF (id == z(ix)) GO TO 878
ix = ix + 1
875 IF (ix <= nxset) GO TO 873
GO TO 871

!     NOW WRITE ENTRY.

878 CALL WRITE (outfl,buf(1),nwds,0)
buf(1) = id
kwds = kwds + nwds
GO TO retx, (640,670,690)

!     IF PLOTS ARE REQUESTED, READ THE CSTM INTO CORE.
!     IF FIRST VECTOR, OPEN PUGV1 AND WRITE HEADER RECORD.

880 CONTINUE
extra = 0
IF (iseq /= 3 .OR. plots == 0 .OR. (kcount /= 1 .AND.  &
    app(1) /= trn(1))) GO TO 990
IF (symflg < 0) GO TO 990
FILE = cstm
CALL OPEN (*900,cstm,z(buf5),rdrew)
CALL fwdrec (*1320,cstm)
icstm = ivecn + 1
CALL READ (*1320,*890,cstm,z(icstm),buf5-icstm,1,ncstm)
CALL mesage (m8,0,nam)
890 CALL CLOSE (cstm,clsrew)
CALL pretrs (z(icstm),ncstm)
900 IF (jcount /= 1) GO TO 902
CALL makmcb (mcb,pugv1,j2,2,qtype2)
FILE = pugv1
CALL OPEN (*902,pugv1,z(buf4),wrtrew)
kplot = 1
CALL fname (pugv1,buf)
CALL WRITE (pugv1,buf,2,1)

!     IF PLOT FILE IS PURGED, NO PLOT FILE CAN BE PREPARED.
!     IF TRANSIENT PROBLEM, REMOVE EXTRA POINTS FROM VECTOR
!     NOW IN CORE THUS CREATING A G-SET VECTOR.

902 extra = 0
IF (kplot == 0) GO TO 990
IF (app(1) /= trn(1) .AND. app(1) /= frq(1) .AND. app(1) /= cei(1)  &
    ) GO TO 910
DO  i = 1,neqex,2
  j = z(i+1)/10
  k = z(i+1) - 10*j
  IF (k /= 3) CYCLE
  extra = 1
  j = ktype*j + ivec - ktype
  z(j) = 1
  IF (ktype == 2) z(j+1) = 1
END DO
IF (extra == 0) GO TO 910
j = ivec
DO  i = ivec,ivecn
  IF (z(i) == 1) CYCLE
  z(j) = z(i)
  j = j + 1
END DO
ivecn = j - 1

!     PASS THE BGPDT. FOR EACH ENTRY, ROTATE THE TRANSLATION COMPONENTS
!     OF UGV TO BASIC (IF REQUIRED). WRITE THESE COMPONENTS ON PUGV1.

910 FILE = bgpdt
CALL OPEN (*990,bgpdt,z(buf5),rdrew)
CALL fwdrec (*1320,bgpdt)
k = 0
i = ivec
pbuff(1) = z(icc+1)
CALL WRITE (pugv1,pbuff,4,1)
l = 3*ktype
CALL bldpk (qtype2,qtype2,pugv1,0,0)
920 CALL READ (*1320,*980,bgpdt,buf(7),4,0,flag)
itemp = 0
DO  j = 1,l
  ll = i + j - 1
  bufr(j) = zz(ll)
END DO
IF (buf(7) < 0.0) THEN
  GO TO   950
ELSE IF (buf(7) == 0.0) THEN
  GO TO   940
END IF

!     TRANSFORM TO BASIC

930 IF (qtype2 == 1) GO TO 935
j = buf(2)
buf(2) = buf(3)
buf(3) = buf(5)
buf(5) = buf(4)
buf(4) = j
935 itemp  = 19
CALL transs (bufr(7),bufr(11))
CALL gmmats (bufr(11),3,3,0,bufr(1),3,1,0,buf(itemp+1))
IF (qtype2 == 1) GO TO 940
CALL gmmats (bufr(11),3,3,0,bufr(4),3,1,0,buf(itemp+4))
j       = buf(21)
buf(21) = buf(23)
buf(23) = buf(24)
buf(24) = buf(22)
buf(22) = j
940 iy = (i-ivec+k)/ktype
DO  j = 1,l,ktype
  iy = iy + 1
  ll = itemp + j
  y(1) = bufr(ll)
  IF (ktype == 2) y(2) = bufr(ll+1)
  CALL zblpki
END DO
i = i + 6*ktype
GO TO 920

!     CHECK FOR FLUID POINTS

950 i = i + ktype
IF (buf(7) /= -2) GO TO 920
iy = (i-ivec+k)/ktype + 2
y(1) = bufr(1)
IF (qtype2 == 3) y(2) = bufr(2)
CALL zblpki
k = k + 5*ktype
GO TO 920
980 CALL bldpkn (pugv1,0,mcb)
CALL CLOSE (bgpdt,clsrew)

!     CONCLUDE PROCESSING OF THIS VECTOR.

990 IF (setno /= 0) CALL WRITE (outfl,0,0,1)
1000 SELECT CASE ( branch )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1020
  CASE (    3)
    GO TO 1170
  CASE (    4)
    GO TO 1020
  CASE (    5)
    GO TO 1040
  CASE (    6)
    GO TO 1110
  CASE (    7)
    GO TO 1170
  CASE (    8)
    GO TO 1020
  CASE (    9)
    GO TO 1020
  CASE (   10)
    GO TO 1020
END SELECT

!     NORMAL STATICS.

1010 IF (jcount < nvects) GO TO 230
IF (eof == 0) GO TO 230
GO TO 1170

!     EIGENVALUES OR DIFF. STIFF PHASE1 OR BUCKLING PHASE 1.

1020 jlist = jlist + ktype1
1030 IF (jcount >= nvects) GO TO 1170
IF (eof == 0) GO TO 230
GO TO 250

!     FREQUENCY RESPONSE.

1040 IF (iseq   <= 2) GO TO 1090
IF (kcount == 3) GO TO 1080
n = ivecn - 1
IF (extra == 0) GO TO 1045
CALL bckrec (infil)
CALL unpack (*1041,infil,z(ivec))
GO TO 1045
1041 DO  i = ivec,n
  zz(i) = 0.0
END DO
GO TO 1055
1045 CONTINUE
omega = twopi*zz(jlist)
DO  i = ivec,n,2
  bufr(1) = -omega*zz(i+1)
  zz(i+1) =  omega*zz(i  )
  zz(i  ) =  bufr(1)
END DO
1055 IF (kcount == 2) GO TO 1060
ireq   = ivel
GO TO 1070
1060 ireq   = iacc
1070 kcount = kcount + 1
incore = 1
GO TO 250
1080 kcount = 1
ireq   = idispl
1090 incore = 0
jlist  = jlist + 2
IF (jlist <= nlist .AND. jcount < nvects) GO TO 250
kfrq   = 0
jlist  = ilist
DO  i = ilist,nlist,2
  z(i+1) = 0
END DO
IF (jcount < nvects) GO TO 230
GO TO 1170

!     TRANSIENT RESPONSE.

1110 IF (iseq <= 2) GO TO 1150
IF (kcount - 2 < 0) THEN
  GO TO  1120
ELSE IF (kcount - 2 == 0) THEN
  GO TO  1130
ELSE
  GO TO  1140
END IF
1120 ireq   = ivel
kcount = 2
GO TO 250
1130 ireq   = iacc
kcount = 3
GO TO 250
1140 ireq   = idispl
kcount = 1
1150 jlist  = jlist + 2
IF (jlist <= nlist .AND. jcount < nvects) GO TO 250
GO TO 1170

!     HERE WHEN END-OF-FILE ENCOUNTERED ON CASE CONTROL.

1160 eof = 1
SELECT CASE ( branch )
  CASE (    1)
    GO TO 1170
  CASE (    2)
    GO TO 1030
  CASE (    3)
    GO TO 1170
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1170
  CASE (    6)
    GO TO 1170
  CASE (    7)
    GO TO 1170
  CASE (    8)
    GO TO 1030
  CASE (    9)
    GO TO 1030
  CASE (   10)
    GO TO 1030
END SELECT

!     CONCLUDE PROCESSING OF CURRENT INPUT FILE.

1170 CALL CLOSE (casecc,clsrew)
CALL CLOSE (infil ,clsrew)
CALL CLOSE (outfl ,clsrew)
IF (kplot /= 0) CALL CLOSE (pugv1,clsrew)
IF (kplot /= 0) CALL wrttrl (mcb)
ocb(2) = kwds/65536
ocb(3) = kwds - 65536*ocb(2)
IF (.NOT.anyout) GO TO 1180
CALL wrttrl (ocb)

!     TEST FOR ALL INPUT FILES PROCESSED.

1180 iseq = iseq + 1
IF (iseq <= 3) GO TO 40
CALL CLOSE (casecc,clsrew)
RETURN

!     HERE IF ABNORMAL CONDITION.
!     CLOSE ALL FILES, JUST TO BE SURE

1190 CALL CLOSE (outfl ,clsrew)
1200 CALL CLOSE (infil ,clsrew)
CALL CLOSE (casecc,clsrew)
ix = iseq + 75
CALL mesage (30,ix,0)
GO TO 1180

!     BINARY SEARCH ROUTINE.
!     =====================

1210 klo = 1
khi = kn
IF (axic) buf(1) = jharm*1000000 + buf(1)
IF (axif > 0.0) THEN
  GO TO  1213
ELSE
  GO TO  1220
END IF
1213 buf(1) = jharm*500000 + buf(1)
1220 k  = (klo+khi+1)/2
1230 kx = 2*k - 1
IF (buf(1)-z(kx) < 0.0) THEN
  GO TO  1240
ELSE IF (buf(1)-z(kx) == 0.0) THEN
  GO TO   700
ELSE
  GO TO  1250
END IF
1240 khi = k
GO TO 1260
1250 klo = k
1260 IF (khi-klo-1 < 0) THEN
  GO TO  1300
ELSE IF (khi-klo-1 == 0) THEN
  GO TO  1270
ELSE
  GO TO  1220
END IF
1270 IF (k == klo) GO TO 1280
k = klo
GO TO 1290
1280 k = khi
1290 klo = khi
GO TO 1230
1300 GO TO retx, (640,670,690)

!     FATAL FILE ERRORS

1310 n = -1
GO TO 1330
1320 n = -2
GO TO 1330
1330 CALL mesage (n,FILE,nam)
GO TO 1330
END SUBROUTINE sdr2c
