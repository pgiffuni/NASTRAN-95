SUBROUTINE sdr2a
     
!     SDR2A PROCESSES THE CASE CONTROL DATA BLOCK. DEPENDING ON THE
!     RIGID FORMAT AND THE VARIOUS OUTPUT REQUESTS, SDR2A SETS FLAGS
!     AND PARAMETERS TO CONTROL OPERATION OF THE REMAINDER OF THE PHASES
!     OF SDR2
 
 EXTERNAL        lshift,rshift
 LOGICAL :: axic  ,ddrmm ,strain
 INTEGER :: z     ,casecc,cstm  ,FILE  ,buf1  ,buf2  ,buf3  ,  &
     buf4  ,buf5  ,all   ,any   ,displ ,spcf  ,stress,  &
     eldef ,any1  ,setno ,zi    ,rshift,strnfl,force ,  &
     prevs ,prevf ,two   ,plots ,ret   ,pass  ,sysbuf,  &
     app   ,sta   ,rei   ,ds0   ,ds1   ,frq   ,trn   ,  &
     bk0   ,bk1   ,cei   ,pla   ,branch,sort2 ,vel   , acc   ,tloads
 DIMENSION       isystm(175)
 COMMON /machin/ mach  ,ihalf ,jhalf
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / app(2),sort2,istrn  ,strnfl,idummy(5)    ,strain
 COMMON /sdr2x1/ ieigen,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym
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
COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
COMMON /system/ sysbuf,opte  ,nogo  ,intap ,mpcn  ,spcn  ,method,  &
    loadnn,symm  ,stftmp,page  ,line  ,tline ,maxlin,  &
    date(3)      ,time  ,echo  ,plots ,ddd(6),mn
COMMON /two   / two(32)
EQUIVALENCE     (sysbuf, isystm)
DATA    mmreig/ 4HMMRE /


!     CHECK FOR STRAIN OPTION

strain = .false.
IF (istrn >= 0) strain = .true.

!     PERFORM BUFFER ALLOCATION.

buf1 = korsz(z) - sysbuf - 2
buf2 = buf1 - sysbuf
buf3 = buf2 - sysbuf
buf4 = buf3 - sysbuf
buf5 = buf4 - sysbuf

!     SET PARAMETER FOR APPROACH.

n = 2*nrigds - 1

!     FIRST CHECK FOR SPECIAL APPROACH FOR DYNAMIC-DATA-RECOVERY-MATRIX-
!     METHOD.  IF APPROACH IS -MMREIG- THEN DDRMM FLAG IS SET TO INSURE
!     ENOUGH OUTPUTS UNDER CERTAIN CONDITIONS.

ddrmm = .false.
IF (app(1) /= mmreig) GO TO 7
ddrmm = .true.
i = 3
GO TO 20

7 DO  i = 1,n,2
  IF (sta(i) == app(1)) GO TO 20
END DO
CALL mesage (-30,75,app)
20 branch = (i+1)/2

!    OPEN CASE CONTROL. SKIP HEADER RECORD.
!    IF DIFF. STIFF. PHASE 1 OR BUCKLING PHASE 1, SKIP 1ST CASECC RECORD

CALL gopen (casecc,z(buf1),rdrew)
IF (app(1) == ds1(1) .OR. app(1) == bk1(1)) CALL skprec (casecc,1)
kwdcc = 0

!     INITIALIZE VARIOUS OUTPUT REQUEST FLAGS.

all   = 0
any   = 0
displ = 0
vel   = 0
acc   = 0
spcf  = 0
loads = 0
stress= 0
force = 0
tloads= 0
eldef = 0
ii    = 0
prevs = 0
prevf = 0

!     READ A RECORD IN CASE CONTROL.
!     IF REQUEST FOR STRESSES IS PRESENT, TURN ON STRESS FLAG.
!     IF REQUEST FOR FORCES   IS PRESENT, TURN ON FORCE  FLAG.
!     -ANY- FLAG = STRESS .OR. FORCE.
!     -ALL- FLAG = ANY REQUEST FOR ALL STRESSES OR FORCES.
!     IF ANY.NE.0 .AND ALL.EQ.0, BUILD LIST OF UNIQUE ELEMENT IDS.

40 CALL READ (*220,*50,casecc,z,buf5-1,1,ncc)
CALL mesage (+8,0,nam)
all  = 1
50 any1 = 0
kwdcc= MAX0(kwdcc,ncc)
mset = MAX0(mset,kwdcc+1)

!     SET DMAP FLAG FOR USE IN DISP R.F. 1

IF (istrn >= 0 .OR. strnfl >= 0) GO TO 55
j = 180
IF (z(j) /= 0) strnfl = 1
55 istr = 23
IF (strain) istr = 180
IF (z(istr) < 0.0) THEN
  GO TO    60
ELSE IF (z(istr) == 0.0) THEN
  GO TO    80
ELSE
  GO TO    70
END IF
60 all = 1
70 stress = 1
any1 = 1
80 IF (z(ielf) < 0.0) THEN
  GO TO    90
ELSE IF (z(ielf) == 0.0) THEN
  GO TO   110
ELSE
  GO TO   100
END IF
90 all = 1
100 force = 1
any1  = 1
110 IF (all /= 0 .OR. any1 == 0) GO TO 200

!     INITIALIZE TO PROCESS STRESS OUTPUT REQUEST.
!     BUILD MASTER SET LIST ONLY IF CURRENT SET ID IS NEW

ASSIGN 190 TO pass
setno = z(istr)
IF (setno == prevs) GO TO 190
prevs = setno

!     IF REQUEST PRESENT, LOCATE SET DEFINITION IN CASE CONTROL DATA.

120 IF (setno == 0) GO TO pass, (190,200)
isetno = ilsym + z(ilsym) + 1
130 iset = isetno + 2
nset = z(isetno+1) + iset - 1
IF (z(isetno) == setno) GO TO 140
isetno = nset + 1
IF (isetno < ncc) GO TO 130
all = 1
GO TO 200

!     PICK UP ELEMENT IDS IN SET. SAVE IN UNIQUE LIST.

140 i = iset
150 IF (i  ==  nset) GO TO 170
IF (z(i+1) > 0) GO TO 170
zi= z(i  )
n =-z(i+1)
i = i + 1
ASSIGN 160 TO ret
GO TO 260
160 zi = zi + 1
IF (zi > n) GO TO 180
ii =ii + 1
IF (ii > buf2) GO TO 280
z(ii) = zi
GO TO 160
170 zi = z(i)
ASSIGN 180 TO ret
GO TO 260
180 i = i + 1
IF (i <= nset) GO TO 150
GO TO pass, (190,200)

!     INITIALIZE TO PROCESS FORCE OUTPUT REQUEST.
!     BUILD MASTER SET LIST ONLY IF CURRENT SET ID IS NEW

190 setno = z(ielf)
IF (setno == prevf) GO TO 200
prevf = setno
ASSIGN 200 TO pass
GO TO 120

!     TURN ON FLAGS FOR OTHER OUTPUT REQUESTS.

200 IF (z(iloads) /= 0) loads = 1
IF (z(ispcf ) /= 0) spcf  = 1
IF (z(idispl) /= 0) displ = 1
IF (z(ivel  ) /= 0) vel   = 1
IF (z(iacc  ) /= 0) acc   = 1
IF (z(ieldef) /= 0) eldef = 1
IF (z(itload) /= 0) tloads= 1
IF (z(iloads+2) < 0 .OR. z(ispcf +2) < 0 .OR.  &
    z(idispl+2) < 0 .OR. z(ivel  +2) < 0 .OR.  &
    z(iacc  +2) < 0 .OR. z(istr  +2) < 0 .OR.  &
    z(ielf  +2) < 0 .OR. app(     1) == trn(1)) sort2 = 1
any = stress + force

!     CONICAL SHELL PROBLEM

axic = .false.
IF (mn == 0) GO TO 210
nrings = isystm(161)
nharms = mn
axic = .true.
210 CONTINUE

!     RETURN TO READ ANOTHER RECORD IN CASE CONTROL (UNLESS DIFF STIFF
!     PHASE 0 OR BUCKLING PHASE 0)

IF (app(1) /= ds0(1) .AND. app(1) /= bk0(1)) GO TO 40

!     IF ALL .EQ. 0, SORT LIST OF ELEMENT IDS AND MOVE LIST TO END OF
!     CORE. AND THROW AWAY ANY DUPLICATE.

220 IF (all /= 0 .OR. any == 0) GO TO 240
kn = ii - mset + 1
CALL sort (0,0,1,-1,z(mset),kn)
jj = buf2 - 1
230 z(jj) = z(ii)
235 ii = ii - 1
IF (z(ii) == z(jj)) GO TO 235
jj = jj - 1
IF (ii >= mset) GO TO 230
mset  = jj + 1
knset = buf2 - mset
GO TO 250
240 mset = buf2 - 1

!     CLOSE CASE CONTROL AND RETURN

250 CALL CLOSE (casecc,clsrew)
IF (app(1) /= bk1(1)) RETURN
eldef = 0
tloads= 0
RETURN


!     SEARCH LIST OF ELEM ID. IF CURRENT ID IS IN LIST RETURN
!     OTHERWISE ADD ID TO LIST


!     ADD ELEM ID TO LIST. NO NEED TO CHECK DUPLICATE ID HERE

260 IF (ii == 0) ii = mset - 1
ii = ii + 1
IF (ii < buf2) GO TO 290
280 all = 1
GO TO 200
290 z(ii) = zi
GO TO ret, (160,180)

END SUBROUTINE sdr2a
