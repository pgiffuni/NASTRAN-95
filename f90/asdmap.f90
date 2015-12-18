SUBROUTINE asdmap
     
!     THIS ROUTINE PROCESSES THE SUBSTRUCTURE COMMAND DATA DECK
 
!     IT CREATES A SET OF SUBSTRUCTURE DATA ON THE FRONT OF THE CASE
!     FILE AND GENERATES DMAP ALTERS FOR THE XALTER FILE. THE ALTERS ARE
!     PLACED FIRST ON THE SCRATCH FILE AND THEN COPIED TO THE PROBLEM
!     TAPE
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        andf        ,orf        ,rshift      ,lshift    , complf
 LOGICAL :: alter       ,altfl      ,first       ,ifin      ,  &
     solve       ,opsof      ,pass2       ,recov     , reject      ,skip
 REAL :: fact        ,xx
 DIMENSION       alts(18)    ,card(20)   ,cdata(30)   ,fname(7)  ,  &
     comnd(2,25) ,dmap(18,60),extra(3,200),ii(9)     ,  &
     nasub(2,100),ocard(200) ,itemp(200)  ,asd1(2)   ,  &
     asd2(2)     ,phs(3)     ,exdef(2,14) ,itmn(5)   ,  &
     ncasec(2)   ,ncases(2)  ,nhead(16)   ,subnam(2) ,  &
     dvec(3)     ,dbvar(6)   ,dbval(2,6,5),r3var(5)  ,  &
     r3val(5,5)  ,var(3,200) ,ivar(3,200) ,z(1)      , corey(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm         ,uwm        ,uim         ,sfm
 COMMON /machin/ mchn
 COMMON /ginox / idum(161)   ,iginob
 COMMON /BLANK / xx
 COMMON /asdbd / irdm        ,nrdm       ,ixtra       ,nxtra     ,  &
     ioct        ,noct       ,iptbs       ,nptbs     ,  &
     iph         ,nph        ,idat(1248)
 COMMON /zzzzzz/ corex(1)
 COMMON /output/ iotit(68)   ,ihead(20)
 COMMON /sofcom/ nsof        ,nname(10)  ,length(10)  ,stat      ,  &
     paswd(2)    ,first      ,opsof
 COMMON /system/ sys(90)     ,lpch
 EQUIVALENCE     (ibuf,sys(1))           ,(outt,sys(2))          ,  &
     (nogo,sys(3))           ,(intp,sys(4))          ,  &
     (nlpp,sys(9))           ,(nlines,sys(12))       ,  &
     (iprec,sys(55))         ,(bandit,sys(77))
 EQUIVALENCE     (var(1,1),ivar(1,1))    ,(itemp(1),ocard(1))    ,  &
     (corex(1),corey(1),ndbs),(corey(2),z(1))
 DATA   alt1   / 4HALTE /,            alt2   / 4HR    /,  &
asd1   / 4HASDM,4HBEGN /,     asd2   / 4HASDM,4HEND /,  &
    BLANK  / 4H     /,            case   / 4HCASE /,  &
    disk   / 4HDISK /,            dolsn  / 4H$    /,  &
    dry    / 4HDRY  /,            drygo  / 4HDRYG /,  &
    ends   / 4HENDS /,            eqsn   / 4H=    /
DATA   exdef  / -1,0,  4HTAPE,4H    ,-1,0   , 4HINTE,4H    , 0,0,  &
    4HNORE,4H    ,4HALL ,4H     , 4HWHOL,4HESOF,  &
    8*4HXXXX,-1,0,   -1 ,0               /
DATA   idg    / 2H0    /,            idry   / 2H-1   /,  &
    inpt   / 4HINPT /,            iopen  / 1      /,  &
    istp   / 2H1    /,            item   / 4HITEM /
DATA   itmn   / 4HITM1,4HITM2,4HITM3,4HITM4,  4HITM5 /
DATA   jnew   / 4HNEW  /,            kbeg   / 4HBEGI /,  &
    kwd    / 4HK    /,            gorun  / 4HGO   /,  &
    lpar   / 4H(    /,            mach   / 4HMACH /,  &
    mwd    / 4HM    /,            NAME   / 4HNAME /,  &
    nano   / 4HNANO /,            nbrec  / 4HBREC /,  &
    ncasec / 4HCASE,4HCC   /,     ncases / 4HCASE,4HSS   /, NEW    / 4HNEW  /
DATA   nhead  / 2*4H    ,4HN a ,4HS t ,4HR a ,4HN  s,4H u b,4H s t  &
    , 4H r u,4H c t,4H u r,4H e  ,4HD e ,4HC k ,4H e c, 4H h o /
DATA   nh1    / 4H1    /,            nh1a   / 4H1A   /,  &
    nphase / 3      /,            nprec  / 4HPREC /,  &
    nrec   / 4HRECO /,            nsave  / 4HSAVE /,  &
    nsol   / 4HSOL  /,            nstp   / 4HNSTP /,  &
    poit   / 4HPOIT /,            nxalt  / 4HXALT /,  &
    nxcsa  / 4HXCSA /,            nxl2   / 4HER   /,  &
    oper   / 4HOPER /,            opti   / 4HOPTI /,  &
    pass   / 4HPASS /,            pass2  / .false./
DATA   phs    / 4HE1  ,4HE2  , 4HE3         /
DATA   posi   / 4HPOSI /,            ptape  / 4HNPTP /,  &
    pwd    / 4HP    /,            run    / 4HRUN  /,  &
    scrt   / 301    /,            sof    / 4HSOF  /,  &
    step   / 4HSTEP /,            titl   / 4HTITL /
DATA   mskp   / 4HMSKP /
DATA   papp   / 4HPAPP /,            pawd   / 4HPA   /,  &
    pitm   / 4HPITM /,            poap   / 4HPOAP /,  &
    pove   / 4HPOVE /,            pvec   / 4HPVEC /
DATA   bwd    / 4HB    /,            k4wd   / 4HK4   /
DATA   outp   / 4HOUTP /,            rang   / 4HRANG /
DATA   dvec   / 4HDVEC,4HUDVF, 4HUDVT       /
DATA   ndbvar / 6      /
DATA   dbvar  / 4HGORL,4HPVEC, 4HUVEC,4HPFTL, 4HOVEC,4HOVC2 /
DATA   dbval  / 4HGEOM,4H4   , 4HPGG ,4H    , 4HUGV ,4H    &
    , 4H    ,4H    , 4HOUGV,4H1   , 4HOUGV,4H    &
    , 4HGEOM,4H4   , 4HPGG ,4H    , 4HUGV ,4H    &
    , 4H    ,4H    , 4HOUGV,4H1   , 4HOUGV,4H    &
    , 4HLAMA,4H    , 4H    ,4H    , 4HPHIG,4H    &
    , 4H    ,4H    , 4HOPHI,4HG1  , 4HOPHI,4HG   &
    , 4HGEOM,4H4   , 4HPPF ,4H    , 4HUGV ,4H    &
    , 4HPPF ,4H    , 4HOUGV,4H1   , 4HOUGV,4H    &
    , 4HGEOM,4H4   , 4HPPT ,4H    , 4HUGV ,4H    &
    , 4HTOL ,4H    , 4HOUGV,4H1   , 4HOUGV,4H     /
DATA   nr3var / 5      /
DATA   r3var  / 4HUAPH,4HPGVC, 4HPSVC,4HDYNT, 4HQVEC /
DATA   r3val  / 4HULV ,4HPGS , 4HPSS ,4H    , 4HQG  ,&
                4HULV ,4HPGS , 4HPSS ,4H    , 4HQG  ,&
                4HPHIA,4H    , 4H    ,4HLAMA, 4HQG  ,&
                4HUDVF,4H    , 4H    ,4HPPF , 4HQPC ,&
                4HUDVT,4HPPT , 4HPST ,4HTOL , 4HQP   /
DATA   ncom   / 25 /
DATA   comnd  / 4HSUBS ,1    , 4HRUN  ,2  &
    , 4HENDD ,2    , 4HCOMB ,3 , 4HREDU ,4    , 4HSOLV ,5  &
    , 4HRECO ,6    , 4HMREC ,6 , 4HBREC ,7    , 4HMRED ,9  &
    , 4HCRED ,13   , 4HDEST ,10 , 4HEDIT ,10   , 4HEQUI ,10  &
    , 4HSOFP ,10   , 4HDELE ,10 , 4HRENA ,10   , 4HSOFI ,11  &
    , 4HSOFO ,11   , 4HREST ,11 , 4HDUMP ,11   , 4HCHEC ,11  &
    , 4HCOMP ,11   , 4HAPPE ,11 , 4HPLOT ,12   /
DATA   subnam / 4HASDM,4HAP  /


CALL conmsg (asd1,2,0)
DO  i  = 64,68
  iotit(i) = BLANK
END DO
DO  i  = 1,16
  ihead(i) = nhead(i)
END DO
CALL page
nz   = korsz(z(1))
buf1 = nz - ibuf + 1
buf2 = buf1 - ibuf
buf3 = buf2 - ibuf

!     INITIALIZE THE CASE CONTROL FILE

CALL OPEN (*2620,case,z(buf2),1)
CALL CLOSE (case,1)
iopen = 1
nopen = buf3 - 1
IF (nopen <= 100) CALL mesage (-8,100-nopen,subnam)
first = .true.
skip  = .false.
isopt = 0

!     SET NUMBER OF POSSIBLE COMMANDS HERE

!     SET LAST WORD INDICATER

i6777 = rshift(complf(0),1)

!     READ FIRST CARD AFTER CEND

ASSIGN 70 TO iread
GO TO 50
30 IF (skip) GO TO 60
IF (nlines >= nlpp) CALL page
nlines = nlines + 1
WRITE (outt,40) card
40 FORMAT (1H ,4X,20A4)
50 CALL xread (*2600,card)
CALL xrcard (ocard,200,card)
IF (ocard(1) > 0 .AND. ocard(2) == BLANK) GO TO 30
IF (ocard(1) == 0) GO TO 30
IF (ocard(2) == titl .OR. ocard(2) == kbeg) GO TO 2600
60 skip = .false.
GO TO iread, (70,90,330,630)
!                      90?  NOT ASSIGNED BY ANYBODY    G.CHAN  4/93

70 IF (ocard(1) > 0 .AND. ocard(2) == comnd(1,1)) GO TO 100

!     NO SUBSTRUCTURE CARD

WRITE (outt,80) ufm
nlines = nlines + 2

80 FORMAT (a23,' 6001. SUBSTRUCTURE DATA IS REQUIRED WITH THIS ', 'APPROACH')
nogo  = 1
90 phase = 2
alter = .false.
skip  = .true.
icom  = 1
IF (ocard(2) == ends) GO TO 2200
GO TO 130

!     PROCESS SUBSTRUCTURE CARD

100 cname = comnd(1,1)
j = ocard(1)*2
DO  i = 1,nphase
  IF (ocard(j+1) /= phs(i)) CYCLE
  phase = i
  alter = .true.
  icom  = 1
  GO TO 130
END DO

!     NO PHASE IS DEFINED

WRITE (outt,120) uwm
nlines = nlines +2
120 FORMAT (a25,' 6002, INCORRECT PHASE DATA')
alter = .false.
nogo  = 1
icom  = 1
phase = 2

!     FOUND PHASE. TURN BANDIT OFF IF PHASE IS 2

130 IF (phase == 2) bandit = -1
j    = 2
iapp = IABS(sys(21))
IF (iapp /= 2) alter = .false.
iap2 = sys(69)/10
IF (iap2 == 1) alter = .false.
sol  = 1
IF (.NOT. alter) GO TO 200
FILE = ptape
kalt = 0
kfile= 0
CALL OPEN (*2620,ptape,z(buf1),0)
140 CALL skpfil (ptape,1)
kfile = kfile + 1
CALL READ (*2620,*150,ptape,fname,7,1,nwords)
150 CONTINUE
IF (fname(1) /= nxalt) GO TO 160
kalt = kfile

GO TO 140
160 IF (fname(1) /= nxcsa) GO TO 140
altfl = .false.
sol = 1
IF (iapp == 3) GO TO 180
CALL READ (*2620,*180,ptape,ii,6,0,nwds)
sol = ii(5)
180 CALL REWIND (ptape)
IF (kalt /= 0) altfl = .true.
IF (altfl) GO TO 190
CALL skpfil (ptape,kfile)
GO TO 200
190 CALL skpfil (ptape,kalt)
CALL fwdrec (*2620,ptape)
CALL READ (*2620,*200,ptape,alts,2,1,nwds)

!     NO XALTER FILE

!     OPEN CASE FILE FOR SUBSTRUCTURE DATA OR TITLE

200 IF (phase == 3) GO TO 300
FILE = case
CALL OPEN (*2620,case,z(buf2),1)
CALL WRITE (case,ncases,2,1)
FILE = scrt

!     SET UP INITAL VALUES

300 iac   = 0
istep = 0
ndbs  = 0
dryflg= 1
obits = 55
IF (sol == 1) obits = 5
IF (sol == 2) obits = 7
IF (sol == 3) obits = 3
IF (sol == 8) obits = 55
IF (sol == 9) obits = 55
newbt = obits
recov = .false.
solve = .false.
iapp  =  sys(21)
IF (iapp == 3) alter = .false.
IF (.NOT.alter ) GO TO 310
CALL OPEN (*2620,scrt,z(buf3),1)
ii(1) = nxalt
ii(2) = nxl2
CALL WRITE (scrt,ii,2,1)
310 CONTINUE
nsof = 0
isof = 1
nname(1) = inpt
stat = 1
length(1)= 100
paswd(1) = BLANK
paswd(2) = BLANK

!     READ PASSWORD AND SOF DECLARATIONS

inex = 0
320 ASSIGN 330 TO iread
GO TO 30
330 IF (ocard(2) /= pass) GO TO 340
k = 4
IF (ocard(5) == eqsn) k = 6
paswd(1) = ocard(k)
paswd(2) = ocard(k+1)
GO TO 320
340 IF (ocard(2) /= sof) GO TO 380
k = 4
IF (ocard(5) /= lpar) GO TO 350
k = 9
isof = ocard(7)
350 IF (isof < 0 .OR. isof > 10) GO TO 370
nsof = nsof + 1
IF (ocard(k+1) == eqsn) k = k + 2
IF (ocard(k+4) == jnew .OR. ocard(k+5) == jnew) stat = 0
nname (isof) = ocard(k  )
length(isof) = ocard(k+3)
IF (ocard(k+2) == -1) GO TO 320
length(isof) = 100
IF (nlines+3 > nlpp) CALL page
nlines = nlines + 3
IF (.NOT.skip) WRITE (outt,40) card
WRITE  (outt,360) uwm,isof
360 FORMAT (a25,', SOF(',i2,') FILESIZE NOT SPECIFIED. DEFAULT OF ',  &
    '100K WORDS WILL BE ALLOCATED',/)
ASSIGN 330 TO iread
IF (skip) GO TO 60
GO TO 50
370 WRITE (outt,790) ufm
nlines = nlines + 1
nogo = 1
GO TO 320
380 IF (inex == 1) GO TO 640
inex = 1
skip = .true.
icnext = 1

!     START PROCESSING SUBSTRUCTURE COMMAND CARDS HERE
!     TOP OF COMMAND LOOP

400 icom = icnext
IF ( ocard(2) == ends) GO TO 2100
DO  l = 1,30
  cdata(l) = ocard(l)
END DO
420 cname = comnd(1,icom)
jcom  = comnd(2,icom)
IF (icom == 6 .AND. sol > 3) jcom = 8
reject = .false.
SELECT CASE ( jcom )
  CASE (    1)
    GO TO 430
  CASE (    2)
    GO TO 440
  CASE (    3)
    GO TO 450
  CASE (    4)
    GO TO 460
  CASE (    5)
    GO TO 470
  CASE (    6)
    GO TO 480
  CASE (    7)
    GO TO 490
  CASE (    8)
    GO TO 500
  CASE (    9)
    GO TO 510
  CASE (   10)
    GO TO 520
  CASE (   11)
    GO TO 530
  CASE (   12)
    GO TO 540
  CASE (   13)
    GO TO 550
END SELECT
430 CALL ascm01 (cname,phase,sol,nogo)
GO TO 600
440 CALL ascm02 (cname,phase,sol,nogo)
GO TO 600
450 CALL ascm03 (cname,phase,sol,nogo)
GO TO 600
460 CALL ascm04 (cname,phase,sol,nogo)
GO TO 600
470 CALL ascm05 (cname,phase,sol,nogo)
GO TO 600
480 CALL ascm06 (cname,phase,sol,nogo)
GO TO 600
490 CALL ascm07 (cname,phase,sol,nogo)
GO TO 600
500 CALL ascm08 (cname,phase,sol,nogo)
GO TO 600
510 CALL ascm09 (cname,phase,sol,nogo)
GO TO 600
520 CALL ascm10 (cname,phase,sol,nogo)
GO TO 600
530 CALL ascm11 (cname,phase,sol,nogo)
GO TO 600
540 CALL ascm12 (cname,phase,sol,nogo)
GO TO 600
550 CALL ascm13 (cname,phase,sol,nogo)
600 jx = 0
istep = istep + 1

!     TRANSFER RAW DMAP TO WORKING AREA

m = irdm - 1
DO  j = 1,nrdm
  DO  i = 1,18
    m = m + 1
    dmap(i,j) = idat(m)
  END DO
END DO

!     READ IN EXTRAS, FIND IN OPTION LIST, STOP AT NEXT COMMAND

620 ASSIGN 630 TO iread
GO TO 30
630 IF (ocard(2) == pass .OR. ocard(2) == sof) GO TO 330
640 IF (reject) GO TO 2090
IF (itemp(2) == dolsn) GO TO 620
IF (itemp(2) ==  ends) GO TO 810
IF (itemp(2) /=  opti) GO TO 660
newbt = 0
i2 = 4
IF (itemp(5) == eqsn) i2 = 6
DO  i = 1,6
  j = 2*i + i2 - 2
  IF (itemp(j) == kwd) newbt = orf(newbt,1)
  IF (itemp(j) == mwd) newbt = orf(newbt,2)
  IF (itemp(j) == pwd) newbt = orf(newbt,4)
  IF (itemp(j) == pawd) newbt = orf(newbt,8)
  IF (itemp(j) == bwd) newbt = orf(newbt,16)
  IF (itemp(j) == k4wd) newbt = orf(newbt,32)
END DO
IF (andf(newbt,12) == 12) GO TO 780
IF (istep <= 1) obits = newbt
GO TO 620
660 CONTINUE
IF (nxtra == 0) GO TO 760
m = ixtra - 1
DO  i = 1,nxtra
  m = m + 1
  IF (itemp(2) == idat(m)) GO TO 680
END DO

!     CARD IS NOT AN EXTRA

GO TO 760

!     FOUND AN EXTRA, STORE SEQUENTIALLY AS PAIRS OF TWO WORD ITEMS

680 jx = jx + 1
extra(1,jx) = itemp(2)
i2 = 4
IF (itemp(5) == eqsn) i2 = 6
extra(2,jx) = itemp(i2  )
extra(3,jx) = itemp(i2+1)

!     SPECIAL OUTPUT EXTRA

IF (itemp(2) /= outp) GO TO 700
extra(2,jx) = -1
extra(3,jx) = 0
690 IF (itemp(i2) /= -1 .OR. itemp(i2+1) <= 0 .OR. itemp(i2+1) > 31)  &
    GO TO 620
j = lshift(1,itemp(i2+1)-1)
extra(3,jx) = orf(extra(3,jx),j)
i2 = i2 + 2
GO TO 690
700 CONTINUE
IF (itemp(2) /= rang) GO TO 720
jx = jx + 1
extra(1,jx) = rang
i2 = i2 + 2
IF (itemp(i2) == -1 .OR. itemp(i2) == -2) GO TO 710
extra(2,jx) = extra(2,jx-1)
extra(3,jx) = extra(3,jx-1)
extra(3,jx-1) = 0
GO TO 620
710 extra(2,jx) = itemp(i2)
extra(3,jx) = itemp(i2+1)
GO TO 620
720 CONTINUE
IF (extra(1,jx) /= run) GO TO 620
extra(3,jx) = BLANK
extra(2,jx) = idry
IF (itemp(i2) ==   dry) GO TO 620
IF (itemp(i2) /= gorun) GO TO 730
extra(2,jx) = idg
GO TO 620
730 IF (itemp(i2) /= step) GO TO 740
extra(2,jx) = istp
GO TO 620
740 IF (itemp(i2) /= drygo) GO TO 750
dryflg = 0
GO TO 620
750 jx = jx - 1
GO TO 620

!      CHECK AND SET IF COMMAND CARD

760 DO  i = 2,ncom
  icnext = i
  IF (ocard(2) == comnd(1,i)) GO TO 800
END DO
780 WRITE (outt,790) ufm
nlines = nlines +2
790 FORMAT (a23,' 6003. ILLEGAL COMMANDS OR OPTIONS DEFINED ON NEXT ', 'CARD')
nogo = 1
GO TO 620

!     FOR PHASE 3 RECOVERY, CHANGE RECO TO BREC

800 IF (phase /= 3 .OR. comnd(1,icnext) /= nrec) GO TO 810
ocard(2) = nbrec
icnext = 9
810 CONTINUE

GO TO ( 820,1030,2200,1100,1200,1300,1400,1400,1700,1200,  &
    1200,1500,1500,1500,1500,1500,1500,1730,1730,1730,  &
    1730,1730,1730,1730,1900), icom

!     SUBSTRUCTURE   PHASES
!         PHASE 1
!      VARIABLES,     NO.     TYPE      POSITION      DEFINITION
!                     1,2,3    I           1    ALTE,R.F. REMOVE NUMBERS
!                     4,5,6    I           4    ALTE,R.F. REMOVE NUMBERS
!                     7,8,9                     ALTE,R.F. REMOVE NUMBERS
!                     10,11,12                  SAVE,-1, PLOT SET ID
!                 1   13,14,15                  RUN ,-1, RUN  FLAG
!                 1   16,17,18                  NAME, SUBS NAME

820 SELECT CASE ( phase )
  CASE (    1)
    GO TO 830
  CASE (    2)
    GO TO 900
  CASE (    3)
    GO TO 1000
END SELECT
830 nvar = 24
nout = 0
DO  i = 1,nvar
  var(i,1) = 0
END DO
var(1,8) = pitm
var(2,8) = pvec
var(3,8) = BLANK
IF (andf(obits,8) /= 0) var(2,8) = papp
DO  i = 1,jx
  DO  j = 1,3
    var(j,i+4) = extra(j,i)
  END DO
END DO
nx   = jx + 4
inam = 0
irun = 0
isav = 0

!     CHECK FOR REQUIRED NAME

DO  i = 5,nx
  IF (var(1,i) ==  NAME) inam = i
  IF (var(1,i) ==   run) irun = i
  IF (var(1,i) == nsave) isav = i
END DO

!     NO NAME DEFINED IS A LEVEL 3 ERROR

IF (inam <= 0) GO TO 2640
IF (irun /= 0) GO TO 870
irun = nx + 1
var(1,irun) =  run
var(2,irun) = istp
var(3,irun) = BLANK
nx = nx + 1
870 CONTINUE
IF (isav /= 0) GO TO 880
var(1,nx+1) = nsave
var(2,nx+1) = -1
var(3,nx+1) = 0
880 CONTINUE
m = iph - 1
DO  i = 1,4
  m = m + 2
  var(1,i) = alt1
  var(2,i) = idat(m-1)
  var(3,i) = idat(m  )
END DO
GO TO 2000

!     PHASE 2 PROCESS

900 IF (jx > 0) GO TO 910
jx = 1
extra(1,1) = run
extra(2,1) = istp
extra(3,1) = BLANK
910 var(1,1) = alt1
var(2,1) = 4
IF (sol == 1) var(2,1) = 5
var(3,1) = 0

DO  j = 1,jx
  DO  i = 1,3
    var(i,j+1) = extra(i,j)
  END DO
END DO
nvar = 3*(1+jx)
nout = 0
DO  i = 1,5
  dmap(1,i) = -1
END DO
GO TO 2000

!     PHASE 3 PROCESSING
!     NORMALLY THIS IS A RESTART, IF NOT THE DATA WILL BE REGENERATED

1000 nvar = 6
var(1,1) = alt1
var(2,1) = idat(iph)
var(3,1) = idat(iph+1)
var(1,2) = run
var(2,2) = istp
var(3,2) = BLANK
IF (jx < 1) GO TO 1010
IF (extra(1,1) == run) var(2,2) = extra(2,1)
1010 nout = 0
DO  i = 1,5
  dmap(1,i) = -1
END DO
GO TO 2000

!     RUN COMMAND (SOMETIMES AN EXTRA)

1030 i2 = 4
IF (cdata(5) == eqsn) i2 = 6
var(1,1) = cdata(2)
var(2,1) = istp
var(3,1) = BLANK
IF (cdata(i2) == step) GO TO 1040
var(2,1) = idry
IF (cdata(i2) == drygo) dryflg = 0
1040 IF (dryflg == 0) GO TO 2080
nvar = 3
nout = 0
GO TO 2000

!     COMBINE OPERATION, USES SUBROUTINE COMBO

1100 CALL combo (cdata,jx,extra,iac,nasub,ns,var(1,3),ier)
nvar = 3*(5+jx+3*ns)
var(1,1) = ns
var(2,1) = 0
var(3,1) = 0
var(1,2) = nstp
var(2,2) =-1
var(3,2) = istep
nvar = nvar + 3
var(nvar+1,1) = pitm
var(nvar+2,1) = pvec
var(nvar+3,1) = BLANK
IF (andf(obits,8) /= 0) var(nvar+2,1) = papp
nvar = nvar + 3
nout = nvar
IF (ier == 0) THEN
  GO TO  2000
ELSE
  GO TO  2640
END IF

!     REDUCE, MREDUCE, CREDUCE OPERATIONS - VARIABLES TO BE SET ARE

!                STEP - STEP NO.
!                NONA - NO. OF SUBSTRUCTURE A
!                NONB - NO. OF SUBSTRUCTURE B
!                NAMA - NAME OF SUBSTRUCTURE A
!                NAMB - NAME OF SUBSTRUCTURE B
!                PREC - PRECISION FLAG
!                PITM - LOAD ITEM
!                POIT - LOAD TRANSFORMATION ITEM

1200 CALL redu (cdata,jx,extra,iac,nasub,nvar,var(1,2),iprec,ier)
var(1,1) = step
var(2,1) = -1
var(3,1) = istep
nvar     = nvar+3
var(nvar+1,1) = pitm
var(nvar+2,1) = pvec
var(nvar+3,1) = BLANK
var(nvar+4,1) = poit
var(nvar+5,1) = pove
var(nvar+6,1) = BLANK
IF (andf(obits,8) == 0) GO TO 1210
var(nvar+2,1) = papp
var(nvar+5,1) = poap
1210 nvar = nvar + 6
nout = nvar
IF (ier == 0) THEN
  GO TO  2000
ELSE
  GO TO  2090
END IF

!     SOLVE OPERATION - VARIABLES ARE SUBSTRUCTURE NAME AND ALTER NO S

1300 nvar = 33
nout = nvar
i2   = 4
IF (cdata(5) == eqsn) i2 = 6
IF (cdata(1)*2 < i2) GO TO 2660
var(1,8) = NAME
var(2,8) = cdata(i2  )
var(3,8) = cdata(i2+1)
nsolv1   = cdata(i2  )
nsolv2   = cdata(i2+1)

!     FIND STRUCTURE NUMBER

ns = iac
IF (ns == 0) GO TO 1320
DO  i = 1,ns
  IF (cdata(i2) == nasub(1,i) .AND. cdata(i2+1) == nasub(2,i)) GO TO 1330
END DO
1320 CONTINUE
ns = ns+1
nasub(1,ns) = cdata(i2  )
nasub(2,ns) = cdata(i2+1)
i = ns
1330 var(1, 9) = nano
var(2, 9) = -1
var(3, 9) = i
var(1,10) = step
var(2,10) = -1
var(3,10) = istep
IF (jcom == 8) GO TO 1340
var(1,11) = nsol
var(2,11) = BLANK
var(3,11) = BLANK
GO TO 1350
1340 var(1,11) = dvec(1)
var(2,11) = dvec(2)
var(3,11) = BLANK
IF (sol == 9) var(2,11) = dvec(3)
1350 CONTINUE
IF (sol == 1) GO TO 1360
IF (sol == 2) var(2,11) = nh1a
IF (sol == 3) var(2,11) = nh1
nvar = 36
nout = 36
var(1,12) = mskp
var(2,12) = BLANK
var(3,12) = BLANK
1360 CONTINUE
iac = ns
m   = iph - 1
DO  i = 1,7
  m = m + 2
  var(1,i) = alt1
  var(2,i) = idat(m-1)
  var(3,i) = idat(m  )
END DO
solve = .true.
GO TO 2000

!     RECOVERY PHASE2 - VARIABLES ARE SOLUTION STRUCTURE NAME,
!                       PRINT, NAME AND/OR SAVE, NAME+ALTER

1400 i2 = 4
IF (cdata(5) == eqsn) i2 = 6
IF (cdata(1)*2 < i2) GO TO 2660
var(1,1) = ncases(1)
var(2,1) = ncases(2)
var(3,1) = BLANK
var(1,2) = NAME
var(2,2) = cdata(i2  )
var(3,2) = cdata(i2+1)
isol = sol
IF (sol  > 3) isol = isol - 4
IF (icom == 8) isol = 3
DO  i = 1,ndbvar
  var(1,i+2) = dbvar(i)
  var(2,i+2) = dbval(1,i,isol)
  IF (icom == 8 .AND. i < 4) var(2,i+2) = BLANK
  var(3,i+2) = dbval(2,i,isol)
  IF (icom == 8 .AND. i < 4) var(3,i+2) = BLANK
END DO
var(1,9) = nsol
var(2,9) = -1
var(3,9) = sol
IF (icom == 8) var(3,9) = 3
var(1,10) = step
var(2,10) = -1
var(3,10) = istep
IF (jx <= 0) GO TO 1430
DO  i = 1,jx
  DO  k = 1,3
    var(k,i+10) = extra(k,i)
  END DO
END DO
1430 IF (solve) GO TO 1440

!     SAVE OPTION BITS AND SET TO ZERO

obits   = 0
var(1,4)= 0
var(2,1)= ncasec(2)
1440 recov   = .true.
nvar    = 3*jx + 30
nout    = nvar
GO TO 2000

!     UTILITY COMMANDS - USE SOFOUT MODULE TO MANIPULATE SOF FILE(S).
!     DESTROY, EDITOUT, EQUIV, PRINT, DELETE, AND RENAME

1500 nvar = 0
i2   = 4
kwds = 1

!     DECODE AND STORE COMMAND DATA FROM HEADER CARD

1510 kwds = kwds + 1
IF (cdata(i2+1) == lpar .OR. cdata(i2+1) == eqsn) i2 = i2 + 2
var(2,kwds) = cdata(i2  )
var(3,kwds) = cdata(i2+1)
i2 = i2 + 2
IF (cdata(i2) == i6777 .OR. cdata(i2+1) == i6777) GO TO 1520
IF (var(2,kwds) == -1) i2 = i2 + 1
IF (kwds < 8) GO TO 1510

!     INSERT VARIABLE NAMES

1520 j = icom - 11
var(1,1) = oper
var(2,1) = cname
var(3,1) = BLANK
jopt = 0
SELECT CASE ( j )
  CASE (    1)
    GO TO 1530
  CASE (    2)
    GO TO 1540
  CASE (    3)
    GO TO 1550
  CASE (    4)
    GO TO 1580
  CASE (    5)
    GO TO 1590
  CASE (    6)
    GO TO 1550
END SELECT

!      DESTROY NAME

1530 var(1,2) = NAME
nvar = 2
GO TO 1620

!     EDITOUT(CODE) = NAME

1540 GO TO 1580

!     EQUIV A,B   +PREFIX = B CARD

1550 var(1,2) = NAME
var(1,3) = NEW
nvar = 3
IF (j    == 6) GO TO 1620
IF (kwds < 2) GO TO 2640
IF (jx   < 1) GO TO 1560
var(1,4) = extra(1,1)
var(2,4) = extra(2,1)
var(3,4) = extra(3,1)
nvar = 4
GO TO 1620
1560 WRITE  (outt,1570) uwm
1570 FORMAT (a25,' 6004, NO PREFIX DEFINED AFTER EQUIVALENCE.')
nlines = nlines + 2
GO TO 1620

!     PRINT(CODE) = NAME,ITM1,ITM2,ITM3,ITM4,ITM5

1580 IF (var(2,2) /= -1) GO TO 1590
var(1,2) = opti
jopt = 1
i2   = 3
GO TO 1600
1590 i2   = 2
1600 var(1,i2) = NAME
ns = kwds - i2
DO  i = 1,ns
  j  = i2 + i
  var(1,j) = itmn(i)
END DO
nvar = kwds
1620 nout = 0
IF (jopt == 1) GO TO 1630
nvar = nvar + 1
var(1,nvar) = opti
var(2,nvar) = -1
var(3,nvar) = 32
IF (icom == 15) var(3,nvar) = 0
1630 nvar = 3*nvar
GO TO 2000

!     RECOVERY, PHASE 3.  VARIABLES ARE NAME, SOL, STEP, PREC, UAPH,
!     PGVC, PSVC, DYNT, QVEC

1700 i2 = 4
IF (cdata(i2) == eqsn) i2 = i2 + 2
m = iph - 1
DO  i = 1,3
  m = m + 2
  var(1,i) = alt1
  var(2,i) = idat(m-1)
  var(3,i) = idat(m  )
END DO
isol = sol
IF (sol > 3) isol = isol - 4
DO  i = 1,nr3var
  var(1,i+3) = r3var(i)
  var(2,i+3) = r3val(i,isol)
  var(3,i+3) = BLANK
END DO
var(1, 9) = NAME
var(2, 9) = cdata(i2  )
var(3, 9) = cdata(i2+1)
var(1,10) = nsol
var(2,10) = -1
var(3,10) = sol
var(1,11) = step
var(2,11) = -1
var(3,11) = istep
var(1,12) = nprec
var(2,12) = -1
var(3,12) = iprec
nvar      = 36
nout      = 0
GO TO 2000

!     EXIO OPERATIONS -
!     SOFIN, SOFOUT, RESTORE, DUMP, CHECK, COMPRESS AND APPEND

1730 nvar = 42
nout = 0

DO  i = 1,14
  var(1,i) = 100 + i
  var(2,i) = exdef(1,i)
  var(3,i) = exdef(2,i)
END DO

!     DECODE COMMAND CARD

var(2,5) = cdata(2)
var(3,5) = cdata(3)
i2 = 4
IF (cdata(5) /= lpar) GO TO 1750
var(2,4) = cdata(6)
var(3,4) = cdata(7)
i2 = 8
1750 IF (cdata(i2+1) == eqsn) i2 = i2 + 2
IF (cdata(i2) == i6777 .OR. cdata(i2+1) == lpar) GO TO 1830
var(2,3) = cdata(i2  )
var(3,3) = cdata(i2+1)
IF (cdata(i2+2) == i6777) GO TO 1760
var(2,2) = cdata(i2+2)
var(3,2) = cdata(i2+3)

!     SET EXTRAS

1760 nn = 0
IF (jx == 0) GO TO 1820
DO  i = 1,jx
  IF (extra(1,i) /= mach) GO TO 1770
  k = 1
  GO TO 1800
  1770 IF (extra(1,i) /= posi) GO TO 1780
  k = 6
  GO TO 1800
  1780 IF (extra(1,i) /= item) GO TO 1790
  k = 7
  GO TO 1800
  1790 IF (extra(1,i) /= NAME) CYCLE
  nn = nn + 1
  k  = nn + 7
  1800 var(2,k) = extra(2,i)
  var(3,k) = extra(3,i)
END DO

!     SET DISK FIELD FOR COMPRESS ETC

1820 IF (icom >= 24) var(2,2) = disk
GO TO 2000
1830 WRITE  (outt,1840) ufm
1840 FORMAT (a23,' 6008, ILLEGAL INPUT ON THE PREVIOUS COMMAND.', /5X,  &
    'MISSING FILE NAME FOR IO OPERATION')
nlines = nlines + 3
GO TO 2090

!     PLOT COMMAND
!          FORMAT
!     PLOT NAME

1900 nvar = 6
nout = 0
i2   = 4
IF (cdata(i2) == eqsn) i2 = 6
IF (cdata(1)*2 <  i2) GO TO 2640
var(1,1) = NAME
var(2,1) = cdata(i2  )
var(3,1) = cdata(i2+1)
var(1,2) = step
var(2,2) = -1
var(3,2) = istep

!     PROCESS VARIABLE CHARACTERS IF DMAP IS TO BE GENERATED

2000 IF (.NOT.alter) GO TO 2080
CALL aspro (dmap,ivar,nvar,obits,sol)

!     RESET OPTION BITS IF DUMMY VALUE WAS USED

obits = newbt
IF (nogo >= 1) GO TO 2080

!     WRITE  DMAP ON  SCRATCH FILE

DO  i = 1,nrdm
  
!     GO TO SPECIAL CODE IF AN ALTER CARD
  
  IF (dmap(1,i) /= alt1) GO TO 2060
  
  ii(1) = dmap(2,i)
  ii(2) = dmap(3,i)
  IF (.NOT.altfl) GO TO 2050
  IF (ii(2)  == 0) ii(2)  = -ii(1)
  2010 IF (alts(2) == 0) alts(2)= -alts(1)
  
  IF (alts(1) > IABS(ii(2))) GO TO 2050
  
!     OVERLAPPING DMAP
  
  IF (IABS(alts(2)) >= ii(1)) GO TO 2660
  
!     ALTERS ENCOUNTERED BEFORE NEW ALTERS
  
  IF (alts(2) < 0) alts(2) = 0
  CALL WRITE (scrt,alts,2,1)
  FILE = ptape
  2020 CALL READ (*2040,*2030,ptape,alts,18,1,nwds)
  
!     DMAP DATA ENCOUNTERED
  
  CALL WRITE (scrt,alts,18,1)
  GO TO 2020
  
!     MORE ALTERS ENCOUNTERED
  
  2030 IF (nwds == 2) GO TO 2010
  
!     END OF USER ALTERS
  
  2040 altfl = .false.
  
!     INSERT NEW DMAP ALTERS
  
  2050 IF (ii(2) < 0) ii(2) = 0
  
  CALL WRITE (scrt,ii,2,1)
  CYCLE
  
!     WRITE ORDINARY DMAP DATA HERE
  
  2060 CALL WRITE (scrt,dmap(1,i),18,1)
END DO

!     WRITE COMMAND AND VARIABLE DATA ON CASE CONTROL FILE

2080 IF (phase == 3) GO TO 2090
ii(1) = cname
ii(2) = nout
CALL WRITE (case,ii,2,0)
CALL WRITE (case,ivar,nout,1)


2090 IF (itemp(2) == ends) GO TO 2100

reject = .false.

GO TO 400

!     ENDSUBS ENCOUNTERED,  STOP PROCESS
!     ENSURE THAT A RECOVER ALWAYS EXISTS FOLLOWING A SOLVE

2100 IF (phase /= 2 .OR. .NOT.solve .OR. recov) GO TO 2110

!     CONSTRUCT A DUMMY INPUT CARD

cdata(1) = 4
cdata(2) = nrec
cdata(3) = BLANK
cdata(4) = nsolv1
cdata(5) = nsolv2
skip = .true.
icom = 7
GO TO 420

!     CHECK SOF AND PASSWORD DECLARATIONS

2110 IF (paswd(1) /= BLANK) IF (nsof) 2120,2120,2140
2120 CALL page2 (2)
WRITE  (outt,2130) ufm
2130 FORMAT (a23,' 6011, SOF DATA PASSWORD MISSING')
nogo = 1
2140 CONTINUE
iblksz = ibuf - 4
IF (mchn == 3 .OR. mchn == 4) iblksz = iginob
fact = 1000.0/iblksz
DO  i = 1,10
  length(i) = length(i)*fact
  jx = length(i)/2
  length(i) = 2*jx
END DO

!     INITIALIZE DIRECT ACCESS FILES FOR IBM 360/370 MACHINES

IF (mchn == 2) CALL sofioi
2200 CALL page2 (1)
WRITE  (outt,2210)
2210 FORMAT (7X,7HENDSUBS)
IF (.NOT. alter) GO TO 2540

!     WRAP UP DMAP
!     PUT LABEL ON END OF ALTER DECK

cname = comnd(1,3)
CALL ascm02 (cname,phase,sol,nogo)
m = irdm + 18
CALL WRITE (scrt,idat(m),18,1)

!     REPEAT ALTER IF DRYGO IS ON

IF (dryflg /= 0) GO TO 2230
DO  i = 1,3
  m = m + 18
  CALL WRITE (scrt,idat(m),18,1)
END DO

!     JUMP TO FINISH OF RIGID FORMAT

2230 IF (phase /= 3) CALL WRITE (scrt,idat(91),18,1)
ifile = ptape
ofile = scrt
ifin  = .false.
pass2 = .false.
IF (.NOT.altfl) GO TO 2240
IF (alts(2) < 0) alts(2)=0
CALL WRITE (scrt,alts,2,1)
GO TO 2270
2240 CALL eof (scrt)
CALL READ (*2620,*2250,ptape,ii,9,1,nw)

!     COPY REMAINDER OF PROBLEM TAPE TO SCRATCH FILE

2250 ifin = .false.
2260 IF (ii(1) == nxcsa) ifin = .true.
CALL WRITE (ofile,ii,nw,1)


irec = 1
2270 CALL READ (*2290,*2280,ifile,z(iopen),nopen,0,nwds)

CALL WRITE (ofile,z(iopen),nopen,0)
GO TO 2270

!     SET ALTER FLAG ON SOL RECORD OF XCSA FILE

2280 IF (ifin .AND. pass2 .AND. irec == 1) z(iopen+2) = 1
CALL WRITE (ofile,z(iopen),nwds,1)
irec = irec + 1
GO TO 2270
2290 CONTINUE
CALL eof (ofile)
IF (ifin) GO TO 2300
CALL READ (*2620,*2260,ifile,ii,9,1,nw)
GO TO 2260
2300 CALL CLOSE (ifile,1)
CALL CLOSE (ofile,3)
IF (pass2) GO TO 2530

!     PRINT OR PUNCH ALTER DECK HERE

!     DIAG 23 REQUESTS PRINT
!     DIAG 24 REQUESTS PUNCH

CALL sswtch (23,kprt)
CALL sswtch (24,kpch)
IF (kprt == 0 .AND. kpch == 0) GO TO 2510
icard = 0
CALL OPEN (*2620,scrt,z(buf3),0)

2310 CONTINUE
CALL page
WRITE  (outt,2320)
2320 FORMAT (5X,'ALTER DECK ECHO')
nlines = nlines + 1
2330 IF (nlines >= nlpp .AND. kprt /= 0) GO TO 2310
CALL READ (*2500,*2360,scrt,card,18,1,nw)

!     DMAP CARD

nc = 18
IF (kprt /= 0) WRITE (outt,2340) icard,(card(i),i=1,nc)
IF (kprt /= 0) nlines = nlines + 1
IF (kpch /= 0) WRITE (lpch,2350) (card(i),i=1,nc)
icard = icard + 1
2340 FORMAT (4X,i5,4X,18A4)
2350 FORMAT (18A4)
GO TO 2330
2360 IF (icard > 0) GO TO 2370
icard = 1
GO TO 2330

!      ALTER CARD

2370 IF (card(2) <= 0) GO TO 2400
IF (kprt /= 0) WRITE (outt,2380) icard,alt1,alt2,(card(i),i=1,2)
IF (kprt /= 0) nlines = nlines + 1
IF (kpch /= 0) WRITE (lpch,2390) alt1,alt2,card(1),card(2)
icard = icard + 1
2380 FORMAT (5X,i4,4X,2A4,i8,1H, ,i3)
2390 FORMAT (2A4,i8,1H,,i3)
GO TO 2330
2400 IF (kprt /= 0) WRITE (outt,2410)icard,alt1,alt2,card(1)
IF (kprt /= 0) nlines = nlines + 1
IF (kpch /= 0) WRITE (lpch,2420) alt1,alt2,card(1)
icard = icard + 1
2410 FORMAT (5X,i4,4X,2A4,i8)
2420 FORMAT (2A4,i8)
GO TO 2330

!     END OF FILE

2500 CALL CLOSE (scrt,0)
2510 CONTINUE
CALL OPEN (*2620,scrt,z(buf3),0)
CALL OPEN (*2620,ptape,z(buf1),0)

!     COPY SCRATCH TO PROB.TAPE, FIRST POSITION PTAPE TO XALTER OR
!     XCSA FILE

iskp = kfile
IF (kalt /= 0) iskp = kalt
CALL skpfil (ptape,iskp)
CALL CLOSE (ptape,2)
CALL OPEN (*2620,ptape,z(buf1),3)
CALL READ (*2620,*2520,scrt,ii,9,1,nw)
2520 pass2 = .true.
ifin  = .false.
ifile = scrt
ofile = ptape
GO TO 2260
2530 CONTINUE

CALL CLOSE (scrt,1)

!     CLOSE CASE CONTROL

2540 IF (phase /= 3) CALL CLOSE (case,2)
CALL conmsg (asd2,2,0)
RETURN

!     USER FATAL MESSAGES

2600 WRITE  (outt,2610) ufm
2610 FORMAT (a23,' 6017, MISSING ENDSUBS CARD.')
CALL mesage (-37,0,subnam)

!     SYSTEM ERROR MESSAGES

2620 WRITE  (outt,2630) sfm,FILE
2630 FORMAT (a25,' 6007, IMPROPER FILE SETUP FOR ',a4)
nlines = nlines + 2
GO TO 2680
2640 WRITE  (outt,2650) ufm,cname
2650 FORMAT (a23,' 6005, ILLEGAL OR MISSING DATA FOR THE PREVIOUS ',  &
    'COMMAND - ',a4)
nlines = nlines + 2
nogo   = 1
GO TO 2090
2660 WRITE  (outt,2670) ufm,alts(1),alts(2),ii
2670 FORMAT (a23,' 6006, DMAP ALTERS  ',2I8, /5X,  &
    'INTERFERE WITH SUBSTRUCTURE ALTERS  ',2I4)
nlines = nlines + 3
nogo   = 1
GO TO 2090
2680 WRITE  (outt,2690) ufm
2690 FORMAT (a23,' 6009, UNRECOVERABLE ERROR CONDITIONS IN SUBROUTINE',  &
    ' ASDMAP')
nlines = nlines + 2
nogo   = 3
CALL CLOSE (scrt ,1)
CALL CLOSE (case ,1)
CALL CLOSE (ptape,1)
RETURN
END SUBROUTINE asdmap
