SUBROUTINE sdcompx (*,zi,zr,zd)
     
!     SDCOMP PERFORMS THE TRIANGULAR DECOMPOSITION OF A SYMMETRIC
!     MATRIX. THE MATRIX MAY BE REAL OR COMPLEX AND ITS PRECISION MAY
!     BE SNGL OR DBL
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(OUT)                     :: zi(1)
 REAL, INTENT(OUT)                        :: zr(2)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(2)
 EXTERNAL lshift  ,andf    ,orf
 LOGICAL :: GO      ,spill   ,splout  ,splin   ,rowone
 INTEGER :: prc     ,words   ,rlcmpx  ,clos    ,buf1    ,buf2    ,  &
     buf3    ,buf4    ,buf5    ,rc      ,prec    ,typea   ,  &
      config  ,power   ,dba     ,dbl     ,dbc     ,  &
     scr1    ,scr2    ,sysbuf  ,forma   ,sym     ,sqr     ,  &
     scra    ,scrb    ,c5max   ,blk     ,pcmax   ,savg    ,  &
     null(20),col     ,c       ,s       ,sprow   ,sturm   ,  &
     groups  ,cavg    ,cmax    ,sc      ,prevc   ,row     ,  &
     frstpc  ,pcavg   ,pcrow   ,pcsqr   ,sx      ,ci      ,  &
     scr3    ,wb      ,scrc    ,scrd    ,spflg   ,start   ,  &
wa      ,chlsky  ,begn    ,END     ,dbname(2)        ,  &
    pcgrou  ,ablk    ,bblk    ,subnam(5)                 ,  &
    key(1)  ,orf     ,statfl  ,andf    ,two24   ,two25   ,  &
    mtype(2),ireal(2),icmplx(2)
REAL :: SAVE(6) ,minds
DOUBLE PRECISION :: mindd   ,xdns(1) ,ddr     ,ddc     , rd      ,dsave3
CHARACTER (LEN=10) :: unuse   ,addi    ,unadd
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON  /xmssg /  ufm     ,uwm     ,uim     ,sfm
COMMON  /sfact /  dba(7)  ,dbl(7)  ,dbc(7)  ,scr1    ,scr2    ,  &
    lcore   ,ddr     ,ddc     ,power   ,scr3    , mindd   ,chlsky
COMMON  /ntime /  nitems  ,tmio    ,tmbpak  ,tmipak  ,tmpak   ,  &
    tmupak  ,tmgstr  ,tmpstr  ,tmt(4)  ,tml(4)
COMMON  /sturmx/  sturm   ,shftpt  ,KEEP    ,ptshft  ,nr
COMMON  /system/  ksystm(63)
COMMON  /names /  rdnrw   ,rdrew   ,wrt     ,wrtrew  ,rew     ,  &
    norew   ,eofnrw  ,rsp     ,rdp     ,csp     ,  &
    cdp     ,sqr     ,rect    ,diag    ,lowtri  , uprtri  ,sym
COMMON  /TYPE  /  prc(2)  ,words(4),rlcmpx(4)
COMMON  /zzzzzz/  xns(1)
COMMON  /sdcomx/  row     ,c       ,spflg   ,start   ,frstpc  ,  &
    lastpl  ,lasti   ,sc      ,iac     ,nzzadr  ,  &
    wa      ,wb      ,prevc   ,nzzz    ,sprow   ,  &
    s       ,blk(15) ,ablk(15),bblk(20)
COMMON  /packx /  itype1  ,itype2  ,i1      ,j1      ,incr1
COMMON  /unpakx/  itype3  ,i2      ,j2      ,incr2
EQUIVALENCE       (nrow,dba(3)) ,(forma,dba(4)) ,(typea,dba(5) ) ,  &
    (jstr,blk(5)) ,(col  ,blk(4)) ,(nterms,blk(6)) ,  &
    (xdns(1),xns(1)),(row,key(1)) ,(dsr  ,ddr    ) ,  &
    (rs  ,rd    ) ,(dsc  ,ddc   ) ,(minds,mindd  )
EQUIVALENCE       (ksystm( 1),sysbuf)  ,(ksystm( 2),nout       ) ,  &
    (ksystm(28),config)  ,(ksystm(40),nbpw       ) ,  &
    (ksystm(57),statfl)  ,(dbname( 1),subnam(4)  )
DATA     subnam/  4HSDCO,2HMP,3*1H   / ,  &
nkey  /  6 / ,  begn/  4HBEGN/ , END   / 4HEND  /      ,  &
    two24 /  16777216   /, two25 /   33554432       /
DATA     ireal ,  icmplx     /  4HREAL,   4H    , 4HCOMP, 4HLEX /
DATA     unuse ,  addi       / '    UNUSED',    'ADDITIONAL'    /

!     STATEMENT FUNCTIONS

nbrwds(i) = i + nwds*(i*(i+1))/2
sx(x)     = x - SQRT(AMAX1(x*(x+2.) + cmax*4. - cons, 1.)) - 1.0
maxc(j)   = SQRT(FLOAT(2*j)/fnwds - FLOAT(4*cmax)) - 1.0

!     BUFFER ALLOCATION

subnam(3) = begn
CALL conmsg (subnam,5,0)
buf1 = lcore- sysbuf
buf2 = buf1 - sysbuf
buf3 = buf2 - sysbuf
buf4 = buf3 - sysbuf
buf5 = buf4 - sysbuf
x    = 1.0
rkhr = 1.0E-10

!     INITIALIZATION AS A FUNCTION OF TYPE OF A MATRIX
!     RC   = 1 IF A IS REAL, 2 IF A IS COMPLEX
!     PREC = 1 IF A IS SINGLE, 2 IF A IS DOUBLE
!     NOTE - PRC(1) = 1, PRC(2) = 2, AND
!            PRC(3) = WORDS(1) = 1, PRC(4) = WORDS(2) = 2

rc = rlcmpx(typea)
mtype(1) = ireal(1)
mtype(2) = ireal(2)
IF (rc == 1) GO TO 10
mtype(1) = icmplx(1)
mtype(2) = icmplx(2)
10 prec  = prc(typea)
nwds  = words(typea)
fnwds = nwds
sturm = 0

!     CHECK INPUT PARAMETERS

IF (dba(2) /= dba(3)) GO TO 2300
icrq = 100 - buf5
IF (buf5 < 100) GO TO 2310
IF (nrow ==   1) GO TO 1900

!     GENERAL INITIALIZATION

loop   = 1
ispill = buf5 - MAX0(100,nrow/100)
fcmax  = 0.
20 ispill = ispill - (loop-1)*nrow/100
nspill = ispill
krow   = nrow + 1
icrq   =-ispill
IF (ispill <= 0) GO TO 2310
zi(ispill) = 0
pcgrou = 0
pcavg  = 0
pcsqr  = 0
pcmax  = 0
csqr   = 0.0
savg   = 0
clos   = ALOG(FLOAT(nrow)) + 5.0
clos   = 999999
pcrow  = -clos
zi(1)  = -nrow
icrq   = nrow - buf5
IF (nrow >= buf5) GO TO 2310
DO  i = 2,nrow
  zi(i)  = 0
END DO
CALL fname (dba,dbname)
power  = 0
scra   = scr3
scrb   = IABS(dbc(1))
GO     =.true.
spill  =.false.
time   = 0.
groups = 0
cmax   = 0
cons   = 2*ispill/nwds
c5max  = maxc(ispill)
dsr    = 1.0
dsc    = 0.
minds  = 1.e+25
IF (prec == 1) GO TO 40
ddr    = 1.0
ddc    = 0.d0
mindd  = 1.d+25
40 CONTINUE
cavg   = 0
cspill = 0.

!     THE FOLLOWING CODE GENERATES THE ACTIVE COLUMN VECTOR FOR EACH
!     ROW, SPILL GROUPS AND TIMING AND USER INFORMATION ABOUT THE
!     DECOMPOSITION

blk(1)  = dba(1)
ablk(1) = scra
ablk(2) = typea
ablk(3) = 0
CALL gopen ( dba,zi(buf1),rdrew)
CALL gopen (scra,zi(buf2),wrtrew)
jlist = 1
row   = 1
jj    = 0
kk    = 0
nlist = 0

!     BEGIN A ROW BY LOCATING THE DIAGONAL ELEMENT

50 blk(8) = -1
kr = krow
60 CALL getstr (*70,blk)
IF (prec == 2) jstr = 2*(jstr-1) + 1
IF (col > row) GO TO 70
IF (col+nterms-1 >= row) GO TO 90
CALL endget (blk)
GO TO 60
70 kk = kk + 1
zi(kk) = row
GO = .false.
80 IF (blk(8) /= 1) CALL skprec (blk,1)
row = row + 1
IF (row <= nrow) GO TO 50
GO TO 600

!     DIAGONAL TERM IS LOCATED - COMPLETE ENTRIES IN THE FULL COLUMN
!     VECTOR AND SAVE THE TERMS FROM EACH STRING IN CORE

90 IF (.NOT. GO) GO TO 80
jstr   = jstr + (row-col)*nwds
nterms = nterms - (row-col)
col  = row
100 zi(kr  ) = col
zi(kr+1) = nterms
kr   = kr + 2
nstr = jstr + nterms*nwds - 1
DO  jj = jstr,nstr
  zr(kr) = xns(jj)
  kr = kr + 1
END DO
n = col + nterms - 1
DO  j = col,n
  IF (zi(j) < 0.0) THEN
    GO TO   120
  ELSE IF (zi(j) == 0.0) THEN
    GO TO   130
  ELSE
    GO TO   160
  END IF
  120 m = IABS(zi(j))
  zi(j) = row
  IF (m /= 1) zi(j+1) = -(m-1)
  CYCLE
  130 i = j
  140 i = i - 1
  IF (i <= 0) GO TO 2000
  IF (zi(i) < 0.0) THEN
    GO TO   150
  ELSE IF (zi(i) == 0.0) THEN
    GO TO   140
  ELSE
    GO TO  2010
  END IF
  150 m = IABS(zi(i))
  zi(i) = -(j-i)
  zi(j) = row
  left  = m - (j-i+1)
  IF (left > 0) zi(j+1) = -left
  CYCLE
  160 IF (zi(j) > row .AND. zi(j) < two24) zi(j) = zi(j) +two24 +two25
END DO
icrq = kr - ispill
IF (kr >= ispill) GO TO 2310
CALL endget (blk)
CALL getstr (*180,blk)
IF (prec == 2) jstr = 2*jstr - 1
GO TO 100

!     EXTRACT ACTIVE COLUMN VECTOR FROM THE FULL COLUMN VECTOR

180 iac = kr
i   = iac
j   = row
lastpl = -1
190 IF (zi(j)     < 0.0) THEN
  GO TO   260
ELSE IF (zi(j)     == 0.0) THEN
  GO TO  2020
END IF
200 IF (zi(j)-row < 0.0) THEN
  GO TO   210
ELSE IF (zi(j)-row == 0.0) THEN
  GO TO   220
ELSE
  GO TO   250
END IF
210 zi(i) = j
GO TO 230
220 zi(i) = -j
IF (lastpl < 0) lastpl = i - iac
230 i = i + 1
240 j = j + 1
GO TO 270
250 IF (zi(j) < two24) GO TO 240
IF (zi(j) < two25) GO TO 210
zi(j) = zi(j) - two25
GO TO 220
260 j = j - zi(j)
270 IF (j <= nrow) GO TO 190
icrq = i - ispill
IF (i > ispill) GO TO 2310
c = i - iac
cmax  = MAX0(cmax,c)
c5max = maxc(ispill)
nac   = iac + c - 1
IF (lastpl < 0) lastpl = c

!     MAKE SPILL CALCULATIONS

spflg = 0
fc    = c
start = 2
IF (c == 1) start = 0
frstpc = 0
IF (.NOT.spill) GO TO 370
IF (row < lstrow) GO TO 290

!     *3* CURRENT ROW IS LAST ROW OF A SPILL GROUP. DETERMINE IF ANOTHER
!         SPILL GROUP FOLLOWS AND, IF SO, ITS RANGE

280 CONTINUE
start = 0
IF (c > c5max) GO TO 380
spill = .false.
GO TO 420

!     *2* CURRENT ROW IS NEITHER FIRST NOR LAST IN CURRENT SPILL GROUP.
!         TEST FOR PASSIVE COL CONDITION. IF SO, TERMINATE SPILL GROUP.
!         TEST FOR POSSIBLE REDEFINITION OF SPILL GROUP. IF SO, TEST FOR
!         OVERFLOW OF REDEFINITION TABLE,  IF SO, TRY A DIFFERENT
!         STRATEGY FOR DEFINING S AND REDO PREFACE UP TO A LIMIT OF 3
!         TIMES.

290 CONTINUE
IF (IABS(zi(iac+1))-row < clos) GO TO 300
ASSIGN 430 TO iswtch
lstrow = row
spill  = .false.
start  = 0
IF (nspill+2 < buf5) GO TO 350
GO TO 330
300 ASSIGN 460 TO iswtch
IF (c <= zi(sprow)) GO TO 460
jj = nac
310 IF (IABS(zi(jj)) <= lstrow) GO TO 320
jj = jj - 1
GO TO 310
320 sc = jj - iac
m  = sx(fc)
IF (sc <= m) GO TO 460
IF (nspill+2 < buf5) GO TO 340
330 CONTINUE
fcmax = AMAX1(fcmax,FLOAT(cmax))
CALL CLOSE (scra,rew)
CALL CLOSE ( dba,rew)
loop = loop + 1
IF (loop <= 3) GO TO 20
icrq = buf5 - nspill - 2
GO TO 2310
340 s = m
ijkl =  MAX0(iac,jj-(sc-m))
lstrow = IABS(zi(ijkl))
350 IF (zi(nspill) /= 0 .AND. zi(nspill) /= sprow) nspill = nspill + 3
zi(nspill  ) = sprow
zi(nspill+1) = s
zi(nspill+2) = lstrow
IF (row- lstrow < 0.0) THEN
  GO TO   360
ELSE IF (row- lstrow == 0.0) THEN
  GO TO   280
ELSE
  GO TO  2070
END IF
360 CONTINUE
GO TO iswtch, (430,460)

!     *1* CURRENT ROW IS NOT PART OF A SPILL GROUP. TEST FOR CREATION OF
!         A NEW SPILL GROUP

370 CONTINUE
IF (c <= c5max) GO TO 420
380 spill  = .true.
sprow  = row
groups = groups + 1
s  = MIN0(sx(fc),nrow-sprow)
IF (loop == 1) GO TO 410
jj = iac + s - 1
390 IF (IABS(zi(jj)) <= sprow+s) GO TO 400
jj = jj - 1
GO TO 390
400 s  = jj - iac + 1
IF (loop == 3) s = MIN0(s,sx(fcmax))
410 s  = MIN0(s,nrow-sprow)
lstrow = IABS(zi(iac+s-1))
spflg  = s
frstpc = lstrow
savg   = savg + s
GO TO 460

!     TEST FOR CONDITION IN WHICH PASSIVE COLUMNS ARE CREATED

420 col = IABS(zi(iac+1))
IF (row-pcrow < clos .OR. c < clos/2 .OR. col-row < clos) GO TO 460

!     CREATE PASSIVE COLUMNS BY CHANGING THEIR FIRST
!     APPEARANCE IN THE FULL COLUMN VECTOR

430 frstpc = 2
pcrow  = row
pcavg  = pcavg + c - 1
pcsqr  = pcsqr + (c-1)**2
pcmax  = MAX0(pcmax,c-1)
pcgrou = pcgrou + 1
nac    = iac + c - 1
ijkl   = iac + 1
DO  i = ijkl,nac
  jj     = IABS(zi(i))
  IF (zi(jj) <= row) GO TO 440
  zi(jj) = MIN0(andf(zi(jj),two24-1),col)
  CYCLE
  440 zi(jj) = col
END DO

!     WRITE ACTIVE COLUMN VECTOR

460 CONTINUE
CALL WRITE (scra,key,nkey,0)
CALL WRITE (scra,zi(iac),c,1)

!     WRITE ROW OF INPUT MATRIX

ablk( 8) = -1
ablk(12) = row
kr = krow
470 ablk(4) = zi(kr  )
nbrstr  = zi(kr+1)
kr = kr + 2
480 CALL putstr (ablk)
ablk(7) = MIN0(ablk(6),nbrstr)
jstr = ablk(5)
IF (prec == 2) jstr = 2*jstr - 1
nstr = jstr + ablk(7)*nwds - 1
DO  jj = jstr,nstr
  xns(jj) = zr(kr)
  kr = kr + 1
END DO
IF (kr >= iac) GO TO 500
CALL endput (ablk)
IF (ablk(7) == nbrstr) GO TO 470
ablk(4) = ablk(4) + ablk(7)
nbrstr  = nbrstr  - ablk(7)
GO TO 480
500 ablk(8) = 1
CALL endput (ablk)

!     ACCUMULATE TIMING AND STATISTICS INFORMATION

cavg = cavg + c
csqr = csqr + c**2
IF (spill) cspill = cspill + c**2
zi(row) = c
IF (row == nrow) GO TO 600
row = row + 1
GO TO 50

!     HERE WHEN ALL ROWS PROCESSED - CLOSE FILES AND, IF SINGULAR
!     MATRIX, PRINT SINGULAR COLUMNS AND GIVE ALTERNATE RETURN

600 CALL CLOSE (scra,rew)
CALL CLOSE ( dba,rew)
IF (GO) GO TO 620
CALL CLOSE (dbl,rew)
CALL page2 (3)
WRITE  (nout,610) ufm,dbname,(zi(i),i=1,kk)
610 FORMAT (a23,' 3097. SYMMETRIC DECOMPOSITION OF DATA BLOCK ',2A4,  &
    ' ABORTED BECAUSE THE FOLLOWING COLUMNS ARE SINGULAR -', /,(5X,20I6,/))
RETURN 1

!     CALCULATE TIME ESTIMATE, PRINT USER INFORMATION AND
!     CHECK FOR SUFFICIENT TIME TO COMPLETE DECOMPOSITION

620 dens  = FLOAT(dba(7))/10000.
IF (dens <  0.01) dens =  0.01
IF (dens > 99.99) dens = 99.99
IF (groups /=  0) savg = savg/groups
savg  = MAX0(savg,1)
time  = 0.5*tmt(typea)*csqr + 0.5*(tmpstr+tmgstr)*FLOAT(pcsqr) +  &
    tmpstr*FLOAT(cavg)  + tmio*(fnwds+1.0)*cspill/FLOAT(savg)
morcor= nbrwds(cmax) - ispill + 1

cavg  = cavg/nrow
IF (pcgrou /= 0) pcavg = pcavg/pcgrou
CALL tmtogo (ijkl)
jklm  = 1.e-6*time + 1.0
icore = IABS(morcor)
IF (dbc(1) <= 0) GO TO 645
unadd = unuse
IF (morcor > 0) unadd = addi
CALL page2 (4)
WRITE (nout,630,ERR=645) uim, mtype, dbname, nrow,   dens,  &
    jklm, cavg,   pcavg, groups, savg, unadd, icore, cmax,   pcmax, pcgrou, loop
630 FORMAT (a29,' 3023 - PARAMETERS FOR ',2A4,  &
    ' SYMMETRIC DECOMPOSITION OF DATA BLOCK ',2A4,  &
    5H (n =,i6, 5H, d =,f6.2,2H%), /14X,  &
    17H  time estimate = , i7, 17H          c avg = , i6,  &
    17H         pc avg = , i6,18H    spill groups = , i6,  &
    17H          s avg = , i6,      /14X,  &
    a10 ,      7H core = , i9,   15H words  c MAX = , i6,  &
    17H          pcmax = , i6,18H       pc groups = , i6,  &
    17H  preface loops = , i6 )
IF (morcor > 0) WRITE (nout,640)
640 FORMAT (14X,'(FOR OPTIMAL OPERATION)')
645 IF (jklm >= ijkl) GO TO 2320

!     WRITE A END-OF-MATRIX STRING ON THE PASSIVE COLUMN FILE

CALL gopen (scrb,zi(buf2),wrtrew)
bblk(1) = scrb
bblk(2) = typea
bblk(3) = 0
bblk(8) =-1
bblk(12)= 1
CALL putstr(bblk)
bblk(4) = nrow + 1
bblk(7) = 1
bblk(8) = 1
CALL endput (bblk)
CALL CLOSE  (scrb,rew)

!     THE STAGE IS SET AT LAST TO PERFORM THE DECOMPOSITION -
!     SO LETS GET THE SHOW UNDERWAY

CALL gopen (scra,zi(buf1),rdrew )
CALL gopen (scrb,zi(buf2),rdrew )
CALL gopen (dbl ,zi(buf3),wrtrew)
scrc   = scr1
scrd   = scr2
IF (zi(nspill) /= 0) nspill = nspill + 3
zi(nspill) = nrow + 1
splin  = .false.
splout = .false.
spill  = .false.
IF (groups /= 0) spill = .true.
nzzz   = orf(ispill-1,1)
rowone = .false.
dbl(2) = 0
dbl(6) = 0
dbl(7) = lshift(1,nbpw-2 - (nbpw-32))

!     THIS 'NEXT TO SIGN' BIT WILL BE PICKED UP BY WRTTRL. ADD (NBPW-32)
!     SO THAT CRAY, WITH 48-BIT INTEGER, WILL NOT GET INTO TROUBLE

blk(1) = dbl(1)
blk(2) = typea
blk(3) = 1
wa     = nzzz
wb     = wa
prevc  = 0
bblk(8)= -1
CALL getstr (*2080,bblk)
kspill = ispill

!     READ KEY WORDS AND ACTIVE COLUMN VECTOR FOR CURRENT ROW

650 NAME = scra
IF (splin) NAME = scrd
CALL fread (NAME,key,nkey,0)
iac = c*nwds + 1
CALL fread (NAME,zi(iac),c,1)
nac = iac + c - 1
IF (zi(iac) < 0) prevc = 0
IF (splin) GO TO 700

!     READ TERMS FROM THE INPUT MATRIX

ablk(8) = -1
CALL getstr (*2090,ablk)
n = iac - 1
DO  i = 1,n
  zr(i) = 0.
END DO
CALL sdcin (ablk,zi(iac),c,zr,zr)

!     IF DEFINED, MERGE ROW FROM PASSIVE COLUMN FILE

680 IF (row-bblk(4) < 0.0) THEN
  GO TO   710
ELSE IF (row-bblk(4) == 0.0) THEN
  GO TO   690
ELSE
  GO TO  2100
END IF
690 CALL sdcin (bblk,zi(iac),c,zr,zr)
bblk(8) = -1
CALL getstr (*2110,bblk)
GO TO 680

!     READ CURRENT PIVOT ROW FROM SPILL FILE. IF LAST ROW, CLOSE FILE

700 prevc = 0
CALL fread (scrd,zr,c*nwds,1)
IF (row < lstspl) GO TO 710
CALL CLOSE (scrd,rew)

!     IF 1ST ROW OF A NEW SPILL GROUP, OPEN SCRATCH FILE TO WRITE

710 IF (rowone) GO TO 740
IF (splout) GO TO 810
IF (spflg == 0) GO TO 810
splout = .true.
CALL gopen (scrc,zi(buf4),wrtrew)
sprow  = row
s      = spflg
lstrow = frstpc
frstpc = 0

!     IF S WAS REDEFINED, GET NEW DEFINITION

DO  i = kspill,nspill,3
  IF (row-zi(i) < 0.0) THEN
    GO TO   740
  ELSE IF (row-zi(i) == 0.0) THEN
    GO TO   730
  ELSE
    GO TO   720
  END IF
END DO
GO TO 740
730 s = zi(i+1)
lstrow = zi(i+2)
kspill = i + 3

!     WRITE ANY TERMS ALREADY CALCULATED WHICH ARE
!     BEYOND THE RANGE OF THE CURRENT SPILL GROUP

740 IF (.NOT.splout) GO TO 810
n    = 0
ijkl = nac
750 IF (IABS(zi(ijkl)) <= lstrow) GO TO 760
ijkl = ijkl - 1
GO TO 750
760 ijkl = ijkl + 1
IF (ijkl > nac) GO TO 780
DO  i = ijkl,nac
  IF (zi(i) > 0.) n = n + 1
END DO
n = nwds*n*(n+1)/2
780 CALL WRITE (scrc,n,1,0)
CALL WRITE (scrc,zr(nzzz-n),n,1)

!     MOVE WA TO ACCOUNT FOR ANY TERMS JUST WRITTEN

IF (n == 0) GO TO 810
j = nzzz
i = nzzz - n
IF (nzzz-wa == n) GO TO 800
790 j = j - 1
i = i - 1
zr(j) = zr(i)
IF (i > wa) GO TO 790
800 wa = j

!     IF THE PIVOTAL ROW DID NOT COME FROM THE SPILL FILE, IT IS CREATED

810 IF (splin) GO TO 1110
i = iac
l = wa
IF (prec == 2) l = (wa-1)/2 + 1
SELECT CASE ( typea )
  CASE (    1)
    GO TO 820
  CASE (    2)
    GO TO 890
  CASE (    3)
    GO TO 960
  CASE (    4)
    GO TO 1030
END SELECT

!     CREATE PIVOT ROW IN RSP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

820 CONTINUE
IF (zi(iac) < 0) GO TO 850
DO  j = 1,c
  IF (zi(i) < 0) GO TO 830
  zr(j) = zr(j) + zr(l)
  l = l + 1
  830 i = i + 1
END DO
850 CONTINUE
ASSIGN 860 TO khr
IF (zr(1) == 0.0) THEN
  GO TO  1820
END IF
860 IF (ABS(dsr) < 10.) GO TO 870
dsr   = dsr/10.
power = power + 1
GO TO 860
870 IF (ABS(dsr) > 0.1) GO TO 880
dsr   = dsr*10.
power = power - 1
GO TO 870
880 dsr   = dsr*zr(1)
minds = AMIN1(ABS(zr(1)),minds)

!     COUNTING SIGN CHANGES OF THE LEADING PRINCIPLE MINORS IN STURM
!     SEQ.

IF (zr(1) < 0.) sturm = sturm + 1
GO TO 1100

!     CREATE PIVOT ROW IN RDP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

890 CONTINUE
IF (zi(iac) < 0) GO TO 920
DO  j = 1,c
  IF (zi(i) < 0) GO TO 900
  zd(j) = zd(j) + zd(l)
  l = l + 1
  900 i = i + 1
END DO
920 CONTINUE
ASSIGN 930 TO khr
IF (zd(1) == 0.0) THEN
  GO TO  1820
END IF
930 IF (DABS(ddr) < 10.0D0) GO TO 940
ddr   = ddr/10.d0
power = power + 1
GO TO 930
940 IF (DABS(ddr) > 0.1D0) GO TO 950
ddr   = ddr*10.d0
power = power - 1
GO TO 940
950 ddr   = ddr*zd(1)
mindd = DMIN1(DABS(zd(1)),mindd)

!     COUNTING SIGN CHANGES (STURM SEQUENCE PROPERTY)

IF (zd(1) < 0.d0) sturm = sturm + 1
GO TO 1100

!     CREATE PIVOT ROW IN CSP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

960 CONTINUE
IF (zi(iac) < 0) GO TO 990
ci = 2*c - 1
DO  j = 1,ci,2
  IF (zi(i) < 0) GO TO 970
  zr(j  ) = zr(j  ) + zr(l  )
  zr(j+1) = zr(j+1) + zr(l+1)
  l = l + 2
  970 i = i + 1
END DO
990 CONTINUE
SAVE(3) = SQRT(zr(1)**2 + zr(2)**2)
IF (SAVE(3)) 1000,1840,1000
1000 IF (SQRT(dsr**2+dsc**2) < 10.) GO TO 1010
dsr = dsr/10.
dsc = dsc/10.
power = power + 1
GO TO 1000
1010 IF (SQRT(dsr**2+dsc**2) > 0.1) GO TO 1020
dsr = dsr*10.
dsc = dsc*10.
power = power - 1
GO TO 1010
1020 rs  = dsr*zr(1) - dsc*zr(2)
dsc = dsr*zr(2) + dsc*zr(1)
drr = rs
minds = AMIN1(SAVE(3),minds)
GO TO 1100

!     CREATE PIVOT ROW IN CDP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

1030 CONTINUE
IF (zi(iac) < 0) GO TO 1060
ci = 2*c - 1
DO  j = 1,ci,2
  IF (zi(i) < 0) GO TO 1040
  zd(j  ) = zd(j  ) + zd(l  )
  zd(j+1) = zd(j+1) + zd(l+1)
  l = l + 2
  1040 i = i + 1
END DO
1060 CONTINUE

!     IN COMPARING THE SOURCE CODES HERE FOR CSP AND CDP COMPUTATION,
!     IT IS DECIDED TO CHANGE THE ORIGINAL LINES (COMMENTED OUT) TO THE
!     NEW LINES USING DSAVE3 INSTEAD OF RD       BY G.CHAN/UNISYS, 8/84

dsave3 = DSQRT(zd(1)**2 + zd(2)**2)
IF (dsave3 == 0.0) THEN
  GO TO  1840
END IF
1070 IF (DSQRT(ddr**2+ddc**2) < 10.d0) GO TO 1080
ddr   = ddr/10.d0
ddc   = ddc/10.d0
power = power + 1
GO TO 1070
1080 IF (DSQRT(ddr**2+ddc**2) > 0.1D0) GO TO 1090
ddr   = ddr*10.d0
ddc   = ddc*10.d0
power = power - 1
GO TO 1080
1090 rd    = ddr*zd(1) - ddc*zd(2)
ddc   = ddr*zd(2) + ddc*zd(1)
ddr   = rd
mindd = DMIN1(dsave3,mindd)

!     CALCULATE WB

1100 CONTINUE
1110 lasti = 1
IF (start == 0) GO TO 1250
IF (splin ) GO TO 1120
IF (splout) GO TO 1130
ci = c
sc = c
GO TO 1160
1120 ci = c - (start-2)
sc = ci
jj = nac
IF (splout) GO TO 1140
IF ((ci*(ci+1)+2*c)*nwds/2+c > nzzz) GO TO 2120
GO TO 1160
1130 ci = c
sc = lstrow - sprow
jj = MIN0(nac,iac+start+sc-2)
1140 IF (IABS(zi(jj)) <= lstrow) GO TO 1150
jj = jj - 1
GO TO 1140
1150 sc = jj - iac - start + 2
IF (sc > 0) GO TO 1160
sc = 0
wb = wa
GO TO 1180
1160 nterms = sc*(ci-1) - (sc*(sc-1))/2
nwords = nterms*nwds
wb = nzzz - nwords
IF (prec == 2) wb = orf(wb-1,1)
IF (wb < iac+c) GO TO 2060
IF (wb > wa+nwds*prevc) GO TO 2130
1180 CONTINUE
IF (splin .AND. row == lstspl) splin = .false.
lasti = MIN0(start+sc-1,c)
IF (sc == 0) GO TO 1250

!     NOW CALCULATE CONTRIBUTIONS FROM CURRENT PIVOT ROW TO SECOND TERM
!     IN EQUATION (4) IN MEMO CWM-19. NOTE-TERMS ARE CALCULATED ONLY
!     FOR ROW/COL COMBINATIONS IN THE CURRENT SPILL GROUP

SELECT CASE ( typea )
  CASE (    1)
    GO TO 1210
  CASE (    2)
    GO TO 1220
  CASE (    3)
    GO TO 1230
  CASE (    4)
    GO TO 1240
END SELECT
1210 CALL sdcom1 (zi,zi(iac),zr(wa+  prevc),zr(wb))
GO TO 1250
1220 CALL sdcom2 (zi,zi(iac),zr(wa+2*prevc),zr(wb))
GO TO 1250
1230 CALL sdcom3 (zi,zi(iac),zr(wa+2*prevc),zr(wb))
GO TO 1250
1240 CALL sdcom4 (zi,zi(iac),zr(wa+4*prevc),zr(wb))

!     SHIP PIVOT ROW OUT TO EITHER MATRIX OR SPILL FILE

1250 IF (lasti == c) GO TO 1290
IF (.NOT. splout) GO TO 2030

!     PIVOT ROW GOES TO SPILL FILE - SET INDEX WHERE TO BEGIN NEXT AND
!     WRITE ROW AND ACTIVE COLUMNN VECTOR

ijkl   = spflg
ii     = frstpc
spflg  = 0
frstpc = 0
start  = lasti + 1
CALL WRITE (scrc,key,nkey, 0)
CALL WRITE (scrc,zi(iac),c,1)
CALL WRITE (scrc,zr,c*nwds,1)
IF (row < lstrow) GO TO 1440

!     LAST ROW OF CURRENT SPILL GROUP - REWIND FILE AND OPEN IT TO READ.
!     IF ANOTHER SPILL GROUP, SET IT UP

CALL CLOSE (scrc,rew)
jklm   = scrc
scrc   = scrd
scrd   = jklm
CALL gopen (scrd,zi(buf5),rdrew)
lstspl = row
splin  =.true.
splout =.false.
IF (ijkl == 0) GO TO 1280
splout =.true.
sprow  = row
s      = ijkl
lstrow = ii
CALL gopen (scrc,zi(buf4),wrtrew)

!     IF S WAS REDEFINED, GET NEW DEFINITION

DO  i = kspill,nspill,3
  IF (row-zi(i) < 0.0) THEN
    GO TO  1280
  ELSE IF (row-zi(i) == 0.0) THEN
    GO TO  1270
  ELSE
    GO TO  1260
  END IF
END DO
GO TO 1280
1270 s = zi(i+1)
lstrow = zi(i+2)
kspill = i + 3

!     READ ANY TERMS SAVED FROM PREVIOUS SPILL GROUP

1280 IF (row == nrow) GO TO 1500
CALL fread (scrd,n,1,0)
wa = nzzz - n
CALL fread (scrd,zr(wa),n,1)
rowone = .true.
GO TO 650

!     PIVOT ROW GOES TO OUTPUT FILE - IF REQUIRED, CONVERT TO CHOLESKY

1290 IF (row /= dbl(2)+1) GO TO 2040
IF (chlsky == 0) GO TO 1340
IF (rc     == 2) GO TO 2050
IF (prec   == 2) GO TO 1320
IF (zr(1) < 0.) GO TO 1800
zr(1) = SQRT(zr(1))
IF (c == 1) GO TO 1340
DO  i = 2,c
  zr(i) = zr(i)*zr(1)
END DO
GO TO 1340
1320 IF (zd(1) < 0.0D+0) GO TO 1800
zd(1) = DSQRT(zd(1))
IF (c == 1) GO TO 1340
DO  i = 2,c
  zd(i) = zd(i)*zd(1)
END DO

!     WRITE THE ROW WITH PUTSTR/ENDPUT

1340 CALL sdcout (blk,0,zi(iac),c,zr,zr)

!     IF ACTIVE COLUMNS ARE NOW GOING PASSIVE, MERGE ROWS IN CORE
!     WITH THOSE NOW ON THE PC FILE THUS CREATING A NEW PC FILE

IF (frstpc == 0) GO TO 1430
IF (splin .OR. splout) GO TO 2140
CALL gopen (scrc,zi(buf4),wrtrew)
blk(1) = scrc
blk(3) = 0
ijkl = iac + 1
DO  i = ijkl,nac
  1360 IF (IABS(zi(i)) <= bblk(4)) GO TO 1380
  CALL cpystr (bblk,blk,1,0)
  bblk(8) = -1
  CALL getstr (*2150,bblk)
  GO TO 1360
  1380 ci = nac - i + 1
  CALL sdcout (blk,0,zi(i),ci,zr(wb),zr(wb))
  wb = wb + ci*nwds
END DO
icrq = wb - ispill
IF (wb > ispill) GO TO 2310
1400 CALL cpystr (bblk,blk,1,0)
IF (bblk(4) == nrow+1) GO TO 1410
bblk(8) = -1
CALL getstr (*2160,bblk)
GO TO 1400
1410 CALL CLOSE (scrb,rew)
CALL CLOSE (scrc,rew)
i = scrb
scrb = scrc
scrc = i
CALL gopen (scrb,zi(buf2),rdrew)
bblk(1) = scrb
bblk(8) = -1
CALL getstr (*2170,bblk)
blk(1) = dbl(1)
blk(3) = 1

!     ACCUMULATE MCB INFORMATION FOR PIVOT ROW

1430 CONTINUE
nwords = c*nwds
dbl(2) = dbl(2) + 1
dbl(6) = MAX0(dbl(6),nwords)
dbl(7) = dbl(7) + nwords

!     PREPARE TO PROCESS NEXT ROW.

1440 IF (row == nrow) GO TO 1500
prevc  = c - 1
rowone = .false.
wa = wb
GO TO 650

!     CLOSE FILES AND PUT END MESSAGE IN RUN LOG.

1500 subnam(3) = END
CALL conmsg (subnam,5,0)
CALL CLOSE (scra,rew)
CALL CLOSE (scrb,rew)
CALL CLOSE ( dbl,rew)

!     PRINT ROOTS INFORMATION IF THIS IS EIGENVALUE PROBLEM, AND KEEP
!     TWO LARGEST SHIFT POINT DATA IF SEVERAL SHIFT POINT MOVINGS ARE
!     INVOLVED.

IF (shftpt > 0.) WRITE (nout,1510) sturm,shftpt
1510 FORMAT (20X,i5,13H roots below ,1P,e14.6)
IF (sturm /= 0) GO TO 1520
IF (KEEP  <= 0) GO TO 1530
sturm  = KEEP
shftpt = ptshft
GO TO 1530
1520 IF (KEEP > sturm) GO TO 1530
jj     = KEEP
rs     = ptshft
KEEP   = sturm
ptshft = shftpt
sturm  = jj
shftpt = rs
1530 IF (statfl /= 1) RETURN

!     PREPARE AND PRINT STATISTICS REGARDING DECOMPOSITION

IF (2*nrow < buf2) GO TO 1600
CALL page2 (2)
WRITE  (nout,1540) uim
1540 FORMAT (a29,' 2316. INSUFFICIENT CORE TO PREPARE DECOMPOSITION ',  &
    'STATISTICS.')
RETURN

1600 CALL gopen (scra,zi(buf1),rdrew)
CALL gopen ( dbl,zi(buf2),rdrew)
ablk(1) = scra
bblk(1) = dbl(1)
row = 1
DO  i = 1,6
  null(i) = 0
END DO
nn = 2*nrow - 1
epsmax = 0.
n  = 0
DO  j = 1,nn,2
  ablk(8) = -1
  bblk(8) = -1
  CALL fwdrec (*2220,ablk)
  CALL getstr (*2180,ablk)
  CALL getstr (*2190,bblk)
  IF (ablk(4) /= row) GO TO 2200
  IF (bblk(4) /= row) GO TO 2210
  ii = ablk(5)
  jj = bblk(5)
  SELECT CASE ( typea )
    CASE (    1)
      GO TO 1660
    CASE (    2)
      GO TO 1670
    CASE (    3)
      GO TO 1680
    CASE (    4)
      GO TO 1690
  END SELECT
  1660 SAVE(2) = xns(ii)
  SAVE(3) = xns(jj)
  GO TO 1700
  1670 SAVE(2) = xdns(ii)
  SAVE(3) = xdns(jj)
  GO TO 1700
  1680 SAVE(2) = SQRT(xns(ii)**2 + xns(ii+1)**2)
  SAVE(3) = SQRT(xns(jj)**2 + xns(jj+1)**2)
  GO TO 1700
  1690 SAVE(2) = DSQRT(xdns(ii)**2 + xdns(ii+1)**2)
  SAVE(3) = DSQRT(xdns(jj)**2 + xdns(jj+1)**2)
  1700 CALL fwdrec (*2220,ablk)
  CALL fwdrec (*2220,bblk)
  eps = ABS(SAVE(2)/SAVE(3))
  zi(j  ) = row
  zi(j+1) = eps
  IF (SAVE(3) < 0.) n = n + 1
  epsmax = AMAX1(epsmax,eps)
  row = row + 1
END DO
CALL sort (0,0,2,2,zi,2*nrow)
CALL CLOSE (ablk,rew)
CALL CLOSE (bblk,rew)
SAVE(1) = 0.1*epsmax
DO  i = 2,6
  SAVE(i) = 0.1*SAVE(i-1)
END DO
DO  j = 1,nn,2
  IF (zr(j+1) > SAVE(1)) GO TO 1730
  IF (zr(j+1) > SAVE(2)) GO TO 1740
  IF (zr(j+1) > SAVE(3)) GO TO 1750
  IF (zr(j+1) > SAVE(4)) GO TO 1760
  IF (zr(j+1) > SAVE(5)) GO TO 1770
  null(6) = null(6) + 1
  CYCLE
  1730 null(1) = null(1) + 1
  CYCLE
  1740 null(2) = null(2) + 1
  CYCLE
  1750 null(3) = null(3) + 1
  CYCLE
  1760 null(4) = null(4) + 1
  CYCLE
  1770 null(5) = null(5) + 1
END DO
i = MAX0(1,nn-8)
CALL page2 (6)
WRITE  (nout,1790) uim,dbname,n,epsmax,(null(j),j=1,6), (zi(j),j=i,nn,2)
1790 FORMAT (a29,' 2314. STATISTICS FOR SYMMETRIC DECOMPOSITION OF ',  &
    'DATA BLOCK ',2A4,7H follow, /10X,23HNUMBER of uii < 0 = ,i5,  &
    /10X,36HMAXIMUM absolute value of aii/uii = ,1P,e12.5,  &
    /10X,13HN1 thru n6 = ,6I6, /10X,36HROW numbers of 5 largest  aii/uii = ,6I6 )
RETURN

!     DIAGONAL ELEMENT .LT. 0.0 IN CHOLESKY DECOMPOSITION

1800 WRITE  (nout,1810) ufm
1810 FORMAT (a23,' 3181, ATTEMPT TO PERFORM CHOLESKY DECOMPOSITION ON',  &
    ' A NEGATIVE DEFINITE MATRIX IN SUBROUTINE SDCOMP.')
GO TO 2330

!     DIAGONAL ELEMENT .EQ. 0.0

1820 zr(1) = rkhr
IF (typea == 2) zd(1) = rkhr
CALL page2 (3)
WRITE  (nout,1830) uwm,row,rkhr
1830 FORMAT (a25,' 2396, SDCOMP COMPUTED A ZERO ON THE DIAGONAL DURING'  &
    ,      ' DECOMPOSITION AT ROW NUMBER',i6,1H., /5X,  &
    'USE OF DIAG 22 OUTPUT SHOULD PERMIT YOU TO CORRELATE THE',  &
    ' ROW WITH A MODEL D.O.F.', /5X,'A VALUE OF ',e13.6,  &
    ' WILL BE USED IN PLACE OF THE ZERO, HOWEVER', /5X,  &
    ' THE ACCURACY OF THE DECOMPOSITION MAY BE IN DOUBT.')
GO TO khr, (860,930)
1840 CALL CLOSE (scra,rew)
CALL CLOSE (scrb,rew)
CALL CLOSE ( dbl,rew)
CALL CLOSE (scrc,rew)
CALL CLOSE (scrd,rew)
RETURN 1

!     DECOMPOSE A 1X1 MATRIX

1900 itype1 = typea
itype2 = typea
itype3 = typea
power  = 0
i1     = 1
j1     = 1
i2     = 1
j2     = 1
incr1  = 1
incr2  = 1
kk     = 1
null(1)= 1
GO     =.false.
CALL gopen  (dba,zi(buf1),rdrew)
CALL unpack (*600,dba,zr)
CALL CLOSE  (dba,rew)
CALL gopen  (dbl,zi(buf1),wrtrew)
dbl(2) = 0
dbl(6) = 0
SELECT CASE ( typea )
  CASE (    1)
    GO TO 1910
  CASE (    2)
    GO TO 1920
  CASE (    3)
    GO TO 1930
  CASE (    4)
    GO TO 1940
END SELECT
1910 minds = zr(1)
dsr   = zr(1)
IF (zr(1) == 0.0) THEN
  GO TO   600
ELSE
  GO TO  1950
END IF
1920 mindd = zd(1)
ddr   = zd(1)
IF (zd(1) == 0.0) THEN
  GO TO   600
ELSE
  GO TO  1950
END IF
1930 minds = SQRT(zr(1)**2 + zr(2)**2)
dsr   = zr(1)
dsc   = zr(2)
IF (minds == 0) THEN
  GO TO   600
ELSE
  GO TO  1950
END IF
1940 mindd = DSQRT(zd(1)**2 + zd(2)**2)
ddr   = zd(1)
ddc   = zd(2)
IF (mindd == 0) THEN
  GO TO   600
END IF
1950 CALL pack  (zr,dbl,dbl)
CALL CLOSE (dbl,rew)
RETURN

!     VARIOUS ERRORS LAND HERE

2000 kerr = 1045
GO TO  2230
2010 kerr = 1046
GO TO  2230
2020 kerr = 1051
GO TO  2230
2030 kerr = 1310
GO TO  2230
2040 kerr = 1320
GO TO  2230
2050 kerr = 1300
GO TO  2230
2060 kerr = 1288
GO TO  2230
2070 kerr = 1065
GO TO  2230
2080 kerr = 1204
GO TO  2230
2090 kerr = 660
GO TO  2230
2100 kerr = 1215
GO TO  2230
2110 kerr = 1216
GO TO  2230
2120 kerr = 1288
GO TO  2230
2130 kerr = 1170
GO TO  2230
2140 kerr = 1350
GO TO  2230
2150 kerr = 1370
GO TO  2230
2160 kerr = 1340
GO TO  2230
2170 kerr = 1420
GO TO  2230
2180 kerr = 1620
GO TO  2230
2190 kerr = 1630
GO TO  2230
2200 kerr = 1640
GO TO  2230
2210 kerr = 1650
GO TO  2230
2220 kerr = 1407
GO TO  2230
2230 WRITE  (nout,2240) sfm,kerr
2240 FORMAT (a25,' 3130, LOGIC ERROR',i6,' OCCURRED IN SDCOMP.')
j = 66
WRITE  (nout,2250) (key(i),i=1,j)
2250 FORMAT (36H0   contents of / sdcomx / follow -- ,/(1X,10I12))
GO TO 2330

!     ERROR EXITS

2300 ier = -7
ifl = 0
GO TO 2340
2310 ier = -8
ifl = icrq
GO TO 2340
2320 ier = -50
ifl = jklm
GO TO 2340
2330 ier = -37
ifl = 0
2340 CALL mesage (ier,ifl,subnam)
RETURN
END SUBROUTINE sdcompx
