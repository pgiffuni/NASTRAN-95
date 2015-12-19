SUBROUTINE sdcmps (zi,zr,zd)
     
!     SDCMPS PERFORMS THE TRIANGULAR DECOMPOSITION OF A SYMMETRIC
!     MATRIX. THE REAL MATRIX INPUT MAY BE SINGLE OR DOUBLE PRECISION.
!     THE OUTPUT MATRICES HAVE POSITIVE DEFINATE CHECKS AND DIAGONAL
!     SINGULARITY CHECKS
 
!     IF SYSTEM(57) IS .GT.1 - USED FOR -CLOS-,
!                      .LT.0 - STOP AFTER PREPASS
 
 
 INTEGER, INTENT(OUT)                     :: zi(1)
 REAL, INTENT(OUT)                        :: zr(1)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(1)
 EXTERNAL          lshift   ,andf     ,orf
 LOGICAL :: spill    ,splout   ,splin    ,rowone   ,opnscr   ,first
 INTEGER :: ablk     ,andf     ,any      ,begn     ,bblk     ,blk    ,  &
     buf1     ,buf2     ,buf3     ,buf4     ,buf5     ,buf6   ,  &
     c        ,cavg     ,chlsky   ,ci       ,clos     ,cmax   ,  &
col      ,c5max    ,dba      ,dbc      ,dbl      ,END    ,  &
    diagck   ,diaget   ,sys60    ,hicore   ,dbname(2)        ,  &
    eor      ,frstpc   ,groups   ,key(1)   ,orf      ,row    ,  &
    parm     ,pcavg    ,pcgrou   ,pcmax    ,pcrow    ,  &
    pcsqr    ,pdefck   ,power    ,prc      ,prec     ,prevc  , rc       ,rlcmpx
INTEGER :: s        ,savg     ,sc       ,scra     ,scrb     ,scrc   ,  &
    scrd     ,scrdia   ,scrmsg   ,scr1     ,scr2     ,scr3   ,  &
    spflg    ,sprow    ,start    ,statfl   ,stscr    ,  &
    sx       ,sysbuf   ,subnam(3),typea    ,two24    ,two25  ,  &
    wa       ,wb       ,words
REAL :: minds    ,SAVE(4)  , ddrr(2)
DOUBLE PRECISION :: ddia     ,ddc      ,ddr      ,dmant    ,dv     ,  &
    mindd    ,pdefd
CHARACTER (LEN=10) :: unuse    ,addi     ,unadd
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /xmssg /   ufm      ,uwm      ,uim      ,sfm
COMMON /machin/   machx
COMMON /lhpwx /   lhpw(6)  ,mtisa
COMMON /sfact /   dba(7)   ,dbl(7)   ,dbc(7)   ,scr1     ,scr2   ,  &
    lcore    ,ddr      ,ddc      ,power    ,scr3   , mindd    ,chlsky
COMMON /ntime /   nitems   ,tmio     ,tmbpak   ,tmipak   ,tmpak  ,  &
    tmupak   ,tmgstr   ,tmpstr   ,tmt(4)   ,tml(4)
COMMON /system/   ksystm(100)
COMMON /names /   rdnrw    ,rdrew    ,wrt      ,wrtrew   ,rew
COMMON /TYPE  /   prc(2)   ,words(4) ,rlcmpx(4)
COMMON /zzzzzz/   xns(1)
COMMON /sdcomx/   row      ,c        ,spflg    ,start    ,frstpc ,  &
    lastpl   ,lasti    ,sc       ,iac      ,nzzadr ,  &
    wa       ,wb       ,prevc    ,nzzz     ,sprow  ,  &
    s        ,blk(15)  ,ablk(15) ,bblk(20)
COMMON /sdcq  /   nerr(2)  ,noglev   ,buf6     ,scrmsg   ,scrdia ,  &
    stscr    ,pdefck   ,diagck   ,diaget   ,prec   , parm(4)  ,opnscr   ,first
COMMON /packx /   itype1   ,itype2   ,i1       ,j1       ,incr1
COMMON /unpakx/   itype3   ,i2       ,j2       ,incr2
EQUIVALENCE       (nrow,dba(3))  ,(typea,dba(5)) ,  &
    (jstr,blk(5))  ,(col  ,blk(4)) ,(nterms,blk(6)),  &
    (row ,key(1))  ,(dsr   ,ddr  ) ,  &
    (dsc,ddc)      ,(minds,mindd ) ,(ddrr(1),rdia ),  &
    (dv,rv) ,(dmant,rmant),(ddia,rdia),(pdefd,pdefr)
EQUIVALENCE       (ksystm( 1),sysbuf) ,(ksystm( 2),nout) ,  &
    (ksystm(31),hicore) ,(ksystm(40),nbpw) , (ksystm(60),sys60 )
DATA    unuse  ,  addi  /'    UNUSED', 'ADDITIONAL' /
DATA    subnam /  4HSDCM,2HPS, 1H   /,  &
nkey   /  6 / ,begn/ 4HBEGN /,  END/ 4HEND  /,  &
    two24  /  16777216 /, two25 /   33554432    /

!     STATEMENT FUNCTIONS

nbrwds(i) = i + nwds*(i*(i+1))/2
sx(x)     = x - SQRT(AMAX1(x*(x+4.0)-cons,0.0)) - 1.0
maxc(j)   = (SQRT(2.*fnwds*FLOAT(j))-3.0)/fnwds

!     VAX, UNIVAC, AND ALL WORKSTATIONS - OPEN CORE CAN BE INCREASED
!     LOCALLY FOR SDCOMP BY SYS60

x = 1.0
korchg = 0
IF (sys60 == 0 .OR. machx == 2 .OR. nbpw > 36) GO TO 20
korchg = sys60 - hicore
IF (korchg <= 0) GO TO 20
lcore = lcore + korchg
WRITE  (nout,10) uim,sys60
10 FORMAT (a29,' - OPEN CORE FOR SDCOMP IS INCREASED TO',i8,  &
    ' WORDS BY SYSTEM(60)',/)
20 IF (lcore <= 0) CALL mesage (-8,0,subnam)

!     BUFFER ALLOCATION

buf1 = lcore- sysbuf
buf2 = buf1 - sysbuf
buf3 = buf2 - sysbuf
buf4 = buf3 - sysbuf
buf5 = buf4 - sysbuf
buf6 = buf5 - sysbuf

!     INITIALIZATION AS A FUNCTION OF TYPE OF A MATRIX
!     RC   = 1 IF A IS REAL (2 IF A IS COMPLEX - ILLEGAL)
!     PREC = 1 IF A IS SINGLE, 2 IF A IS DOUBLE

rc    = rlcmpx(typea)
IF (rc /= 1) GO TO 1600
statfl= IABS(ksystm(57))
prec  = prc(typea)
nwds  = words(typea)
fnwds = nwds

!     CHECK INPUT PARAMETERS

IF (dba(2) /= dba(3)) GO TO 1600
icrq = nrow + 200 - buf6
IF (icrq > 0) GO TO 1850

!     INITIALIZE POSITIVE DEFINATE CHECKS.  FILES SET IN DRIVER

parm(1) = 0
parm(3) = subnam(1)
parm(4) = subnam(2)
nerr(1) = 0
nerr(2) = 0
IF (pdefck < 0) GO TO 50
i = -diaget
j = 1 - mtisa
IF (prec == 2) GO TO 30
pdefr = 2.0E0**i
rmant = 2.0E0**j
GO TO 50
30 pdefd = 2.0D0**i
dmant = 2.0D0**j
GO TO 50
50 CONTINUE

!     STSCR IS STATUS OF -SCRDIA- FILE AT BUF6
!       0 = NOT OPEN
!       1 = READ
!       2 = WRITE

stscr = 2
CALL gopen (scrdia,zi(buf6),wrtrew)
scra = scr3
scrb = IABS(dbc(1))
noglev = 0
IF (nrow == 1) GO TO 1510

!     GENERAL INITIALIZATION

loop  = 1
ispill= buf6 - MAX0(100,nrow/100)
fcmax = 0.
60 ispill= ispill - (loop-1)*nrow/100
nspill= ispill
krow  = nrow + 1
icrq  = (3-loop)*nrow/100 - ispill
IF (ispill <= 0) GO TO 1850
zi(ispill) = 0
pcgrou= 0
pcavg = 0
pcsqr = 0
pcmax = 0
csqr  = 0.0
savg  = 0
clos  = ALOG(FLOAT(nrow)) + 5.0
IF (statfl > 1) clos = statfl
pcrow = -clos
zi(1) = -nrow
DO  i = 2,nrow
  zi(i) = 0
END DO
CALL fname (dba,dbname)
power = 0
spill = .false.
groups= 0
cons  = 2*ispill/nwds
c5max = maxc(ispill)
dsr   = 1.0
dsc   = 0.
minds = 1.e+25
IF (prec == 1) GO TO 80
ddr   = 1.0D0
ddc   = 0.d0
mindd = 1.d+25
80 CONTINUE
cavg  = 0
cmax  = 0
cspill= 0.0

!     THE FOLLOWING CODE GENERATES THE ACTIVE COLUMN VECTOR FOR EACH
!     ROW, SPILL GROUPS AND TIMING AND USER INFORMATION ABOUT THE
!     DECOMPOSITION

blk(1)  = dba(1)
ablk(1) = scra
ablk(2) = typea
ablk(3) = 0
CALL gopen (dba ,zi(buf1),rdrew )
CALL gopen (scra,zi(buf2),wrtrew)
row = 1
jj  = 0
eor = 1

!     LSTDIA DETERMINES THE LAST DIAGONAL WRITTEN TO SCRATCH FILE

lstdia = 0

!     BEGIN A ROW BY LOCATING THE DIAGONAL ELEMENT

90 blk(8) = -1

!     ANY DETERMINES IF ANY STRINGS SKIPPED PRIOR TO DIAGONAL
!     AND -KK- ALLOWS STRING BEYOND ZERO DIAGONAL TO BE SAVED

any = 0
kr  = krow
100 CALL getstr (*110,blk)
IF (prec == 2) jstr = 2*(jstr-1) + 1
kk  = nterms
any = col
IF (col > row) GO TO 130
kk  = 0
IF (col+nterms-1 >= row) GO TO 140
CALL endget (blk)
GO TO 100

!     NULL COLUMN FOUND.  SAVE COLUMN ID AND SET NOGLEV

110 kk = -1
IF (any /= 0) GO TO 130
IF (lstdia < row) CALL sdcmq (*710,1,0.,0.,0.d0,0.d0,row,zi)
120 IF (blk(8) /= 1) CALL fwdrec (*1680,blk)
row = row + 1
IF (row <= nrow) GO TO 90
GO TO 710

!     ZERO DIAGONAL FOUND.  FILL CORE AND POINTERS

130 col = row
zi(kr  ) = col
zi(kr+1) = 1
zi(kr+2) = 0
IF (nwds == 2) zi(kr+3) = 0
kr = kr + 2 + nwds
nterms = nwds
IF (lstdia >= row) GO TO 140
ddia = 0.0D0
CALL sdcmq (*710,7,0.,0.,0.d0,0.d0,row,zi)
IF (noglev > 1) GO TO 120
CALL WRITE (scrdia,rdia,nwds,eor)
lstdia = row
GO TO 180

!     DIAGONAL TERM IS LOCATED -- COMPLETE ENTRIES IN THE FULL COLUMN
!     VECTOR AND SAVE THE TERMS FROM EACH STRING IN CORE

140 CONTINUE
jstr = jstr + (row-col)*nwds
IF (lstdia >= row) GO TO 150
rdia = xns(jstr)
IF (prec == 2) ddrr(2) = xns(jstr+1)
IF (noglev <= 1) CALL WRITE (scrdia,rdia,nwds,eor)
lstdia = row
150 CONTINUE
nterms = nterms - (row-col)
col = row
160 zi(kr  ) = col
zi(kr+1) = nterms
kr = kr + 2
nstr = jstr + nterms*nwds - 1
DO  jj = jstr,nstr
  zr(kr) = xns(jj)
  kr = kr + 1
END DO
180 CONTINUE
n = col + nterms - 1
DO  j = col,n
  IF (zi(j) < 0.0) THEN
    GO TO   190
  ELSE IF (zi(j) == 0.0) THEN
    GO TO   200
  ELSE
    GO TO   230
  END IF
  190 m = IABS(zi(j))
  zi(j) = row
  IF (m /= 1) zi(j+1) = -(m-1)
  CYCLE
  200 i = j
  210 i = i - 1
  IF (i <= 0) GO TO 1610
  IF (zi(i) < 0.0) THEN
    GO TO   220
  ELSE IF (zi(i) == 0.0) THEN
    GO TO   210
  ELSE
    GO TO  1620
  END IF
  220 m = IABS(zi(i))
  zi(i) = -(j-i)
  zi(j) = row
  left  = m - (j-i+1)
  IF (left > 0) zi(j+1) = -left
  CYCLE
  230 IF (zi(j) > row .AND. zi(j) < two24) zi(j) = zi(j) +two24 +two25
END DO
icrq = kr - ispill
IF (kr >= ispill) GO TO 700

!     CHECK IF ZERO DIAGONAL WAS JUST PROCESSED

IF (kk < 0) THEN
  GO TO   270
ELSE IF (kk == 0) THEN
  GO TO   250
ELSE
  GO TO   260
END IF
250 CALL endget (blk)
CALL getstr (*280,blk)
IF (prec == 2) jstr = 2*jstr - 1
GO TO 160
260 col = any
nterms = kk
kk = 0
GO TO 140

!     EXTRACT ACTIVE COLUMN VECTOR FROM THE FULL COLUMN VECTOR

270 IF (blk(8) /= 1) CALL fwdrec (*1680,blk)
280 iac = kr
i = iac
j = row
lastpl = -1
290 IF (zi(j) < 0.0) THEN
  GO TO   360
ELSE IF (zi(j) == 0.0) THEN
  GO TO  1630
END IF
300 IF (zi(j)-row < 0.0) THEN
  GO TO   310
ELSE IF (zi(j)-row == 0.0) THEN
  GO TO   320
ELSE
  GO TO   350
END IF
310 zi(i) = j
GO TO 330
320 zi(i) = -j
IF (lastpl < 0) lastpl = i - iac
330 i = i + 1
340 j = j + 1
GO TO 370
350 IF (zi(j) < two24) GO TO 340
IF (zi(j) < two25) GO TO 310
zi(j) = zi(j) - two25
GO TO 320
360 j = j - zi(j)
370 IF (j <= nrow) GO TO 290
icrq = i - ispill
IF (i > ispill) GO TO 700
c = i - iac
cmax = MAX0(cmax,c)
nac = iac + c - 1
IF (lastpl < 0) lastpl = c

!     MAKE SPILL CALCULATIONS

spflg = 0
fc    = c
start = 2
IF (c == 1) start = 0
frstpc = 0
IF (.NOT. spill) GO TO 490
IF (row < lstrow) GO TO 410

! *3* CURRENT ROW IS LAST ROW OF A SPILL GROUP. DETERMINE IF ANOTHER
!     SPILL GROUP FOLLOWS AND, IF SO, ITS RANGE

400 CONTINUE
start = 0
IF (c > c5max) GO TO 500
spill = .false.
GO TO 540

! *2* CURRENT ROW IS NEITHER FIRST NOR LAST IN CURRENT SPILL GROUP.
!     TEST FOR PASSIVE COL CONDITION. IF SO, TERMINATE SPILL GROUP.
!     TEST FOR POSSIBLE REDEFINITION OF SPILL GROUP. IF SO, TEST FOR
!     OVERFLOW OF REDEFINITION TABLE,  IF SO, TRY A DIFFERENT STRATEGY
!     FOR DEFINING S AND REDO PREFACE UP TO A LIMIT OF 3 TIMES.

410 CONTINUE
IF (IABS(zi(iac+1))-row < clos) GO TO 420
ASSIGN 550 TO iswtch
lstrow= row
spill = .false.
start = 0
IF (nspill+2 < buf6) GO TO 470
GO TO 450
420 ASSIGN 580 TO iswtch
IF (c <= zi(sprow)) GO TO 580
jj = nac
430 IF (IABS(zi(jj)) <= lstrow) GO TO 440
jj = jj - 1
GO TO 430
440 sc = jj - iac
m  = sx(fc)
IF (sc <= m) GO TO 580
IF (nspill+2 < buf6) GO TO 460
450 CONTINUE
fcmax = AMAX1(fcmax,FLOAT(cmax))
CALL CLOSE (scra,rew)
CALL CLOSE (dba ,rew)
loop = loop + 1
IF (loop <= 3) GO TO 60
icrq = buf6 - nspill - 3
GO TO 1850
460 s = m
ijkl = MAX0(iac,jj - (sc-m))
lstrow = IABS(zi(ijkl))
470 IF (zi(nspill) /= 0 .AND. zi(nspill) /= sprow) nspill = nspill + 3
zi(nspill  ) = sprow
zi(nspill+1) = s
zi(nspill+2) = lstrow
IF (row-lstrow < 0.0) THEN
  GO TO   480
ELSE IF (row-lstrow == 0.0) THEN
  GO TO   400
ELSE
  GO TO  1670
END IF
480 CONTINUE
GO TO iswtch, (550,580)

! *1* CURRENT ROW IS NOT PART OF A SPILL GROUP. TEST FOR
!     CREATION OF A NEW SPILL GROUP

490 CONTINUE
IF (c <= c5max) GO TO 540
500 spill = .true.
sprow = row
groups= groups + 1
s = MIN0(sx(fc),nrow-sprow)
IF (loop == 1) GO TO 530
jj = iac + s - 1
510 IF (IABS(zi(jj)) <= sprow+s) GO TO 520
jj = jj - 1
GO TO 510
520 s = jj - iac + 1
IF (loop == 3) s = MIN0(s,sx(fcmax))
530 s = MIN0(s,nrow-sprow)
lstrow = IABS(zi(iac+s-1))
spflg  = s
frstpc = lstrow
savg   = savg + s
GO TO 580

!     TEST FOR CONDITION IN WHICH PASSIVE COLUMNS ARE CREATED

540 col = IABS(zi(iac+1))
IF (row-pcrow < clos .OR. c < clos/2 .OR. col-row < clos) GO TO 580

!     CREATE PASSIVE COLUMNS BY CHANGING THEIR FIRST
!     APPEARANCE IN THE FULL COLUMN VECTOR

550 frstpc= 2
pcrow = row
pcavg = pcavg + c - 1
pcsqr = pcsqr + (c-1)**2
pcmax = MAX0(pcmax,c-1)
pcgrou= pcgrou + 1
nac   = iac + c - 1
ijkl  = iac + 1
DO  i = ijkl,nac
  jj = IABS(zi(i))
  IF (zi(jj) <= row) GO TO 560
  zi(jj) = MIN0(andf(zi(jj),two24-1),col)
  CYCLE
  560 zi(jj) = col
END DO

!     WRITE ACTIVE COLUMN VECTOR

580 IF (noglev > 1) GO TO 630
CALL WRITE (scra,key,nkey,0)
CALL WRITE (scra,zi(iac),c,1)

!     WRITE ROW OF INPUT MATRIX. -IAC- POINTS TO END OF OUTPUT

ablk(8)  = -1
ablk(12) = row
kr = krow
590 ablk(4)= zi(kr)
nbrstr = zi(kr+1)
kr = kr + 2
600 CALL putstr (ablk)
ablk(7) = MIN0(ablk(6),nbrstr)
jstr = ablk(5)
IF (prec == 2) jstr = 2*jstr - 1
nstr = jstr + ablk(7)*nwds - 1
DO  jj = jstr,nstr
  xns(jj) = zr(kr)
  kr = kr + 1
END DO
IF (kr >= iac) GO TO 620
CALL endput (ablk)
IF (ablk(7) == nbrstr) GO TO 590
ablk(4) = ablk(4) + ablk(7)
nbrstr  = nbrstr  - ablk(7)
GO TO 600
620 ablk(8) = 1
CALL endput (ablk)

!     ACCUMULATE TIMING AND STATISTICS INFORMATION

630 cavg = cavg + c
csqr = csqr + c**2
IF (spill) cspill = cspill + c**2
zi(row) = c
IF (row == nrow) GO TO 710
row = row + 1
GO TO 90

!     HERE WHEN ALL ROWS PROCESSED -- CLOSE FILES AND, IF SINGULAR
!     MATRIX, PRINT SINGULAR COLUMNS AND GIVE ALTERNATE RETURN

700 parm(1) = -8
parm(2) = icrq
noglev  = 2
710 CALL CLOSE (scra,rew)
CALL CLOSE (dba ,rew)
CALL CLOSE (scrdia,rew)

!     CALCULATE TIME ESTIMATE, PRINT USER INFORMATION AND
!     CHECK FOR SUFFICIENT TIME TO COMPLETE DECOMPOSITION

IF (groups /= 0) savg = savg/groups
savg    = MAX0(savg,1)
SAVE(1) = 0.5*tmt(typea)*csqr*1.0E-6
SAVE(2) = 0.5*(tmpstr+tmgstr)*FLOAT(pcsqr)*1.e-6
SAVE(3) = tmpstr*FLOAT(cavg)*1.e-6
SAVE(4) = tmio*(fnwds+1.0)*cspill/FLOAT(savg)*1.0E-6
morcor  = nbrwds(cmax) - ispill + 1

cavg = cavg/nrow
IF (pcgrou /= 0) pcavg = pcavg/pcgrou
CALL tmtogo (ijkl)
jklm = SAVE(1) + SAVE(2) + SAVE(3) + SAVE(4) + 1.0

IF (dbc(1) > 0) CALL page2 (9)
unadd = unuse
IF (morcor > 0) unadd = addi
IF (dbc(1) > 0) WRITE (nout,720)  uim,   dbname,   nrow,  &
    jklm,   cavg,  pcavg, groups,   savg,  &
    unadd,      morcor,   cmax,  pcmax, pcgrou,   loop
720 FORMAT (a29,' 3023 - PARAMETERS FOR SYMMETRIC DECOMPOSITION OF ',  &
    'DATA BLOCK ',2A4,6H ( n = , i5, 2H ) , /  &
    14X, 17H  time estimate = , i7, 17H          c avg = , i6,  &
    17H         pc avg = , i6,18H    spill groups = , i6,  &
    17H          s avg = , i6, /  &
    14X, a10 ,      7H core = , i7, 17H words    c MAX = , i6,  &
    17H          pcmax = , i6,18H       pc groups = , i6,  &
    17H  preface loops = , i6  )
IF (morcor > 0) WRITE (nout,730)
730 FORMAT (15X,'(FOR OPTIMIZED OPERATION)')
IF (dbc(1) > 0) WRITE (nout,740) uim,subnam(1),subnam(2),SAVE
740 FORMAT (a29,' 2378,',a4,a3,' ESTIMATE OF CPU TIME FOR MT =',  &
    1P,e10.3,/18X,'PASSIVE COL. = ',e10.3,14X,'ACTIVE COL. =',  &
    e10.3, /25X,'SPILL = ',e10.3)

!     ESTIMATE FBS TIME AT ONE PASS, 1 LOAD

SAVE(1) = 2.0*FLOAT(nrow)*cavg*(tmt(typea)+tmpstr)*1.e-6
IF (dbc(1) > 0) WRITE (nout,750) SAVE(1)
750 FORMAT (10X,41HESTIMATE for fbs, one pass AND one load =,1P,e10.3)

IF (jklm >= ijkl) GO TO 1840
IF (noglev >  1) GO TO 1880
IF (ksystm(57) < 0) GO TO 1880

!     WRITE A END-OF-MATRIX STRING ON THE PASSIVE COLUMN FILE

CALL gopen (scrb,zi(buf2),wrtrew)
bblk(1) = scrb
bblk(2) = typea
bblk(3) = 0
bblk(8) = -1
CALL putstr (bblk)
bblk(4) = nrow + 1
bblk(7) = 1
bblk(8) = 1
CALL endput (bblk)
CALL CLOSE (scrb,rew)
subnam(3) = begn
CALL conmsg (subnam,3,0)

!     THE STAGE IS SET AT LAST TO PERFORM THE DECOMPOSITION -
!     SO LETS GET THE SHOW UNDERWAY

CALL gopen (scra,zi(buf1),rdrew )
CALL gopen (scrb,zi(buf2),rdrew )
CALL gopen (dbl ,zi(buf3),wrtrew)
CALL gopen (scrdia,zi(buf6),rdrew)
stscr = 1
scrc  = scr1
scrd  = scr2
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
CALL getstr (*1690,bblk)
kspill = ispill

!     READ KEY WORDS AND ACTIVE COLUMN VECTOR FOR CURRENT ROW

800 NAME = scra
IF (splin) NAME = scrd
CALL fread (NAME,key,nkey,0)
iac = c*nwds + 1
CALL fread (NAME,zi(iac),c,1)
nac = iac + c - 1
IF (zi(iac) < 0) prevc = 0
IF (splin) GO TO 840

!     READ TERMS FROM THE INPUT MATRIX

CALL fread (scrdia,rdia,nwds,eor)
ablk(8) = -1
CALL getstr (*1860,ablk)
n = iac - 1
DO  i = 1,n
  zr(i) = 0.
END DO
CALL sdcins (*1830,ablk,zi(iac),c,zr,zd)

!     IF DEFINED, MERGE ROW FROM PASSIVE COLUMN FILE

820 IF (row-bblk(4) < 0.0) THEN
  GO TO   850
ELSE IF (row-bblk(4) == 0.0) THEN
  GO TO   830
ELSE
  GO TO  1700
END IF
830 CALL sdcins (*1830,bblk,zi(iac),c,zr,zd)
bblk(8) = -1
CALL getstr (*1710,bblk)
GO TO 820

!     READ CURRENT PIVOT ROW FROM SPILL FILE. IF LAST ROW, CLOSE FILE

840 prevc = 0
CALL fread (scrd,zr,c*nwds,1)
IF (row < lstspl) GO TO 850
CALL CLOSE (scrd,rew)

!     IF 1ST ROW OF A NEW SPILL GROUP, OPEN SCRATCH FILE TO WRITE

850 IF (rowone) GO TO 880
IF (splout) GO TO 950
IF (spflg == 0) GO TO 950
splout = .true.
CALL gopen (scrc,zi(buf4),wrtrew)
sprow = row
s = spflg
lstrow = frstpc
frstpc = 0

!     IF S WAS REDEFINED, GET NEW DEFINITION

DO  i = kspill,nspill,3
  IF (row-zi(i) < 0.0) THEN
    GO TO   860
  ELSE IF (row-zi(i) == 0.0) THEN
    GO TO   870
  ELSE
    GO TO   880
  END IF
  860 CONTINUE
END DO
GO TO 880
870 s = zi(i+1)
lstrow = zi(i+2)
kspill = i + 3

!     WRITE ANY TERMS ALREADY CALCULATED WHICH ARE
!     BEYOND THE RANGE OF THE CURRENT SPILL GROUP

880 IF (.NOT. splout) GO TO 950
n = 0
ijkl = nac
890 IF (IABS(zi(ijkl)) <= lstrow) GO TO 900
ijkl = ijkl - 1
GO TO 890
900 ijkl = ijkl + 1
IF (ijkl > nac) GO TO 920
DO  i = ijkl,nac
  IF (zi(i) > 0.) n = n + 1
END DO
n = nwds*n*(n+1)/2
920 CALL WRITE (scrc,n,1,0)
CALL WRITE (scrc,zr(nzzz-n),n,1)

!     MOVE WA TO ACCOUNT FOR ANY TERMS JUST WRITTEN

IF (n == 0) GO TO 950
j = nzzz
i = nzzz - n
IF ((nzzz-wa) == n) GO TO 940
930 j = j - 1
i = i - 1
zr(j) = zr(i)
IF (i > wa) GO TO 930
940 wa = j

!     IF THE PIVOTAL ROW DID NOT COME FROM THE SPILL FILE, IT IS CREATED

950 IF (splin) GO TO 1180
i = iac
l = wa
IF (prec  == 2) l = (wa-1)/2 + 1
IF (typea == 2) GO TO 1060

!     CREATE PIVOT ROW IN RSP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

IF (zi(iac) < 0) GO TO 980
DO  j = 1,c
  IF (zi(i) < 0) GO TO 960
  zr(j) = zr(j) + zr(l)
  l = l + 1
  960 i = i + 1
END DO
980 CONTINUE

!     CHECK DIAGONAL AND CORRECT

IF (zr(1) == 0.0) CALL sdcmq (*1870,2,rdia,zr(1),0,0,row,zi)
990 IF (ABS(dsr) < 10.) GO TO 1000
dsr   = dsr/10.
power = power + 1
GO TO 990
1000 IF (ABS(dsr) > 0.1) GO TO 1010
dsr   = dsr*10.
power = power - 1
GO TO 1000
1010 dsr   = dsr*zr(1)
minds = AMIN1(zr(1),minds)

!     PERFORM MATRIX COND. CHECKS - S.P. REAL

IF (zr(1) < 0.0) THEN
  GO TO  1020
ELSE IF (zr(1) == 0.0) THEN
  GO TO  1030
ELSE
  GO TO  1050
END IF
1020 i = 3
GO TO 1040
1030 i = 2
1040 CALL sdcmq (*1870,i,rdia,zr(1),0,0,row,zi)

1050 IF (diagck < 0) GO TO 1170
IF (rdia == 0.0) rdia = zr(1)
IF (rdia == zr(1)) GO TO 1170
rv = ABS(zr(1)/rdia )
IF (rv > 1.001E0) CALL sdcmq (*1870,6,rdia,zr(1),0,0,row,zi)
rv = rmant/rv
IF (rv > pdefr) CALL sdcmq (*1870,4,rdia,zr(1),0,0,row,zi)
GO TO 1170

!     CREATE PIVOT ROW IN RDP, ACCUMULATE DETERMINANT AND MIN DIAGONAL

1060 CONTINUE
IF (zi(iac) < 0) GO TO 1090
DO  j = 1,c
  IF (zi(i) < 0) GO TO 1070
  zd(j) = zd(j) + zd(l)
  l = l + 1
  1070 i = i + 1
END DO
1090 CONTINUE

!     CHECK DIAGONAL AND CORRECT

IF (zd(1) == 0.0D0) CALL sdcmq (*1870,2,0,0,ddia,zd(1),row,zi)
1100 IF (DABS(ddr) < 10.0D0) GO TO 1110
ddr   = ddr/10.d0
power = power + 1
GO TO 1100
1110 IF (DABS(ddr) > 0.1D0) GO TO 1120
ddr   = ddr*10.d0
power = power - 1
GO TO 1110
1120 ddr   = ddr*zd(1)
mindd = DMIN1(zd(1),mindd)

!     PERFORM MATRIX COND. CHECKS - D.P. REAL

IF (zd(1) < 0.0) THEN
  GO TO  1130
ELSE IF (zd(1) == 0.0) THEN
  GO TO  1140
ELSE
  GO TO  1160
END IF
1130 i = 3
GO TO 1150
1140 i = 2
1150 CALL sdcmq (*1870,i,0,0,ddia,zd(1),row,zi)

1160 IF (diagck <   0) GO TO 1170
IF (ddia == 0.0D0) ddia = zd(1)
IF (ddia == 0.0D0) GO TO 1170
dv = DABS(zd(1)/ddia)
IF (dv > 1.001D0) CALL sdcmq (*1870,6,0,0,ddia,zd(1),row,zi)
dv = dmant/dv
IF (dv > pdefd) CALL sdcmq (*1870,4,0,0,ddia,zd(1),row,zi)

!     CALCULATE WB

1170 CONTINUE
1180 lasti = 1
IF (start == 0) GO TO 1260
IF (splin ) GO TO 1190
IF (splout) GO TO 1200
ci = c
sc = c
GO TO 1230
1190 ci = c - (start-2)
sc = ci
jj = nac
IF (splout) GO TO 1210
IF (ci > c5max) GO TO 1720
GO TO 1230
1200 ci = c
sc = lstrow - sprow
jj = MIN0(nac,iac+start+sc-2)
1210 IF (IABS(zi(jj)) <= lstrow) GO TO 1220
jj = jj - 1
GO TO 1210
1220 sc = jj - iac - start + 2
IF (sc > 0) GO TO 1230
sc = 0
wb = wa
GO TO 1240
1230 nterms = sc*(ci-1) - (sc*(sc-1))/2
nwords = nterms*nwds
wb = nzzz - nwords
IF (prec == 2) wb = orf(wb-1,1)
IF (wb < iac+c) GO TO 1660
IF (wb > wa+nwds*prevc) GO TO 1730
1240 CONTINUE
IF (splin .AND. row == lstspl) splin = .false.
lasti = MIN0(start+sc-1,c)
IF (sc == 0) GO TO 1260

!     NOW CALCULATE CONTIBUTIONS FROM CURRENT PIVOT ROW TO
!     SECOND TERM IN EQUATION (4) IN MEMO CWM-19. NOTE-TERMS ARE
!     CALCULATED ONLY FOR ROW/COL COMBINATIONS IN THE CURRENT SPILL
!     GROUP

IF (typea == 2) GO TO 1250
CALL sdcom1 (zi,zi(iac),zr(wa+prevc),zr(wb))
GO TO 1260
1250 CALL sdcom2 (zi,zi(iac),zr(wa+2*prevc),zr(wb))

!     SHIP PIVOT ROW OUT TO EITHER MATRIX OR SPILL FILE

1260 IF (lasti == c) GO TO 1300
IF (.NOT. splout) GO TO 1640

!     PIVOT ROW GOES TO SPILL FILE - SET INDEX WHERE TO BEGIN NEXT AND
!                                    WRITE ROW AND ACTIVE COLUMN VECTOR

ijkl  = spflg
ii    = frstpc
spflg = 0
frstpc= 0
start = lasti + 1
CALL WRITE (scrc,key,nkey, 0)
CALL WRITE (scrc,zi(iac),c,1)
CALL WRITE (scrc,zr,c*nwds,1)
IF (row < lstrow) GO TO 1410

!     LAST ROW OF CURRENT SPILL GROUP - REWIND FILE AND OPEN IT TO READ.
!                                      IF ANOTHER SPILL GROUP, SET IT UP

CALL CLOSE (scrc,rew)
jklm = scrc
scrc = scrd
scrd = jklm
CALL gopen (scrd,zi(buf5),rdrew)
lstspl = row
splin  = .true.
splout = .false.
IF (ijkl == 0) GO TO 1290
splout = .true.
sprow  = row
s      = ijkl
lstrow = ii
CALL gopen (scrc,zi(buf4),wrtrew)

!     IF S WAS REDEFINED, GET NEW DEFINITION

DO  i = kspill,nspill,3
  IF (row-zi(i) < 0.0) THEN
    GO TO  1270
  ELSE IF (row-zi(i) == 0.0) THEN
    GO TO  1280
  ELSE
    GO TO  1290
  END IF
  1270 CONTINUE
END DO
GO TO 1290
1280 s = zi(i+1)
lstrow = zi(i+2)
kspill = i + 3

!     READ ANY TERMS SAVED FROM PREVIOUS SPILL GROUP

1290 IF (row == nrow) GO TO 1500
CALL fread (scrd,n,1,0)
wa = nzzz - n
CALL fread (scrd,zr(wa),n,1)
rowone = .true.
GO TO 800

!     PIVOT ROW GOES TO OUTPUT FILE - IF REQUIRED, CONVERT TO CHOLESKY

1300 IF (row /= dbl(2)+1) GO TO 1650
IF (chlsky == 0) GO TO 1340
IF (prec   == 2) GO TO 1320
IF (zr(1) < 0.) CALL sdcmq (*1870,3,rdia,zr(1),0,0,row,zi)
zr(1) = SQRT(zr(1))
IF (c == 1) GO TO 1340
DO  i = 2,c
  zr(i) = zr(i)*zr(1)
END DO
GO TO 1340
1320 IF (zd(1) < 0.0D0) CALL sdcmq (*1870,3,0,0,ddia,zd(1),row,zi)
zd(1) = DSQRT(zd(1))
IF (c == 1) GO TO 1340
DO  i = 2,c
  zd(i) = zd(i)*zd(1)
END DO

!     WRITE THE ROW WITH PUTSTR/ENDPUT

1340 CALL sdcout (blk,0,zi(iac),c,zr,zr)

!     IF ACTIVE COLUMNS ARE NOW GOING PASSIVE, MERGE ROWS IN CORE
!     WITH THOSE NOW ON THE PC FILE THUS CREATING A NEW PC FILE

IF (frstpc == 0) GO TO 1400
IF (splin .OR. splout) GO TO 1740
CALL gopen (scrc,zi(buf4),wrtrew)
blk(1) = scrc
blk(3) = 0
ijkl   = iac + 1
DO  i = ijkl,nac
  1350 IF (IABS(zi(i)) <= bblk(4)) GO TO 1360
  CALL cpystr (bblk,blk,1,0)
  bblk(8) = -1
  CALL getstr (*1750,bblk)
  GO TO 1350
  1360 ci = nac - i + 1
  CALL sdcout (blk,0,zi(i),ci,zr(wb),zr(wb))
  wb = wb + ci*nwds
END DO
icrq = wb - ispill
IF (wb > ispill) GO TO 1850
1380 CALL cpystr (bblk,blk,1,0)
IF (bblk(4) == nrow+1) GO TO 1390
bblk(8) = -1
CALL getstr (*1760,bblk)
GO TO 1380
1390 CALL CLOSE (scrb,rew)
CALL CLOSE (scrc,rew)
i    = scrb
scrb = scrc
scrc = i
CALL gopen (scrb,zi(buf2),rdrew)
bblk(1) = scrb
bblk(8) = -1
CALL getstr (*1770,bblk)
blk(1) = dbl(1)
blk(3) = 1

!     ACCUMULATE MCB INFORMATION FOR PIVOT ROW

1400 CONTINUE
nwords = c*nwds
dbl(2) = dbl(2) + 1
dbl(6) = MAX0(dbl(6),nwords)
dbl(7) = dbl(7) + nwords

!     PREPARE TO PROCESS NEXT ROW.

1410 IF (row == nrow) GO TO 1500
prevc = c - 1
rowone= .false.
wa    = wb
GO TO 800

!     CLOSE FILES AND PUT END MESSAGE IN RUN LOG.

1500 subnam(3) = END
CALL conmsg (subnam,3,0)
GO TO 1870

!     DECOMPOSE A 1X1 MATRIX

1510 itype1= typea
itype2= typea
itype3= typea
power = 0
i1    = 1
j1    = 1
i2    = 1
j2    = 1
incr1 = 1
incr2 = 1
CALL gopen (dba,zi(buf1),rdrew)
parm(2) = dba(1)
CALL unpack (*1570,dba,zr)
CALL CLOSE (dba,rew)
CALL gopen (dbl,zi(buf1),wrtrew)
dbl(2) = 0
dbl(6) = 0
IF (typea == 2) GO TO 1520
minds = zr(1)
dsr   = zr(1)
IF (zr(1) < 0.0) THEN
  GO TO  1530
ELSE IF (zr(1) == 0.0) THEN
  GO TO  1540
ELSE
  GO TO  1560
END IF
1520 mindd = zd(1)
ddr   = zd(1)
IF (zd(1) < 0.0) THEN
  GO TO  1530
ELSE IF (zd(1) == 0.0) THEN
  GO TO  1540
ELSE
  GO TO  1560
END IF

1530 i = 3
GO TO 1550
1540 i = 2
1550 CALL sdcmq (*1870,i,zr,zr,zd,zd,1,zi)
1560 CALL pack (zr,dbl,dbl)
CALL CLOSE (dbl,rew)
GO TO 1880

!     1X1 NULL COLUMN

1570 CALL sdcmq (*1870,1,0.,0.,0.d0,0.d0,1,zi)
GO TO 1870

!     VARIOUS ERRORS LAND HERE

1600 CALL mesage (-7,dba(2),subnam)
1610 kerr = 1045
GO TO  1800
1620 kerr = 1046
GO TO  1800
1630 kerr = 1051
GO TO  1800
1640 kerr = 1311
GO TO  1800
1650 kerr = 1320
GO TO  1800
1660 kerr = 1288
GO TO  1800
1670 kerr = 1065
GO TO  1800
1680 kerr = 1034
GO TO  1800
1690 kerr = 1204
GO TO  1800
1700 kerr = 1215
GO TO  1800
1710 kerr = 1216
GO TO  1800
1720 kerr = 1288
GO TO  1800
1730 kerr = 1289
GO TO  1800
1740 kerr = 1330
GO TO  1800
1750 kerr = 1333
GO TO  1800
1760 kerr = 1340
GO TO  1800
1770 kerr = 1344
GO TO  1800
1800 WRITE  (nout,1810) sfm,kerr
1810 FORMAT (a25,' 2379, LOGIC ERROR',i6,' IN SDCMPS.')
j = 66
WRITE  (nout,1820) (key(i),i=1,j)
1820 FORMAT (36H0   contents of / sdcomx / follow --, /,(1X,10I12))
1830 parm(1) = -37
parm(2) = 0
parm(3) = subnam(1)
parm(4) = subnam(2)
GO TO 1870

!     INSUFFICIENT TIME

1840 parm(1) = -50
parm(2) = ijkl
GO TO 1870

!     INSUFFICIENT CORE

1850 parm(1) = -8
parm(2) = icrq
GO TO 1870

!     UNEXPECTED NULL COLUMN

1860 dv = 0.0
CALL sdcmq (*1870,5,rv,rv,dv,dv,row,zi)

1870 CALL CLOSE (dba, rew)
CALL CLOSE (scra,rew)
CALL CLOSE (scrb,rew)
CALL CLOSE (dbl ,rew)
1880 CALL CLOSE (scrdia,rew)
IF (nerr(1)+nerr(2) <= 0) GO TO 1890
CALL gopen (scrmsg,zi(buf6),wrt)
bblk(2) = 0
bblk(3) = 0
bblk(4) = 0
CALL WRITE (scrmsg,bblk(2),3,1)
CALL CLOSE (scrmsg,rew)
1890 CONTINUE
IF (korchg > 0) lcore = lcore - korchg
RETURN
END SUBROUTINE sdcmps
