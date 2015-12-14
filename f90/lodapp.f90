SUBROUTINE lodapp
     
!     THIS MODULE APPENDS NEW LOAD VECTORS (PAPP AND POAP) TO THE
!     SUBSTRUCTURE  -NAME-.  THE NEW VECTORS ARE MERGED WITH ALREADY
!     EXISTING (PVEC AND POVE) MATRICES OR ARE SIMPLY RECOPIED AS
!     THE NEW PVEC AND POVE ITEMS. LOAP DATA IS ALSO MERGED WITH THE
!     LODS DATA OR IS SIMPLY COPIED AS THE NEW LODS ITEM.
 
!     SOF ITEMS -
 
!     LOAP - APPENDED LOAD SET IDENTIFICATION TABLE
!     PAPP - APPENDED LOAD MATRICES (G-SET)
!     POAP - APPENDED LOAD MATRICES (O-SET)
!     LODS - LOAD SET IDENTIFICATION TABLE **BECOMES THE NEW LODS**
!     PVEC - LOAD MATRICES (G-SET)         **BECOMES THE NEW PVEC**
!     POVE - LOAD MATRICES (O-SET)         **BECOMES THE NEW POVE**
 
 EXTERNAL        rshift     ,andf
 LOGICAL :: lpapp      ,lpvec      ,lpoap      ,lpove      ,  &
     lmerg      ,llsub      ,mdiup      ,ditup
 INTEGER :: rshift     ,andf       ,buf
 DIMENSION       iz(1)      ,nn(2)      ,mcbloc(7)  ,adump(4000),  &
     nprog(2)   ,NAME(2)    ,namell(2)  ,icorx(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm        ,uwm        ,uim        ,sfm
 COMMON /BLANK / buf(3)
 COMMON /system/ isbuff     ,lp
 COMMON /packx / itypin     ,itypot     ,ifirst     ,ilast      , incr
 COMMON /parmeg/ mcbk(7)    ,mcbk11(7)  ,mcbk21(7)  ,mcbk12(7)  ,  &
     mcbk22(7)  ,lcore      ,irule
 COMMON /zzzzzz/ rz(1)
 COMMON /names / ird        ,irdrew     ,iwrt       ,iwrtrw     ,  &
     irew       ,inorew     ,iefnrw     ,irsp       ,  &
     irdp       ,icsp       ,icdp       ,isqure     ,  &
     irect      ,idiag      ,iupper     ,ilower     , isym
 COMMON /sof   / ditdum(6)  ,iodum(8)   ,mdidum(4)  ,nxtdum(15) ,  &
     ditup      ,mdiup
 EQUIVALENCE     (icorx(1)  ,rz(1))
 EQUIVALENCE     (rz(1)     ,iz(1))
 DATA    ipapp , ipoap      /101        ,102        /
 DATA    iscr1 , iscr2      ,iscr3      ,iscr4      ,iscr5      ,  &
     iscr6 , iscr7      ,iscr8      /  &
     301   , 302        ,303        ,304        ,305        ,  &
     306   , 307        ,308        /
 DATA    npapp , npoap      ,npvec      ,npove      , nloap , nlods      /  &
     4HPAPP, 4HPOAP     ,4HPVEC     ,4HPOVE     , 4HLOAP, 4HLODS     /
 DATA    nprog              /4HLODA     ,4HPP       /
 DATA    blnk  / 4H         /
 
!     INITIALIZE PARAMETERS
 
 CALL tmtogo (itime1)
 NAME(1) = buf(1)
 NAME(2) = buf(2)
 idry    = buf(3)
 ncore = korsz(iz(1))
 
!     INITIALIZE OPEN CORE - THERE ARE NIZ WORDS AVAILABLE
 
 ib1  = ncore - (isbuff+1)
 ib2  = ib1 - isbuff - 1
 ib3  = ib2 - isbuff
 ibuf1= ib3 - isbuff
 niz  = ibuf1 - 1
 nstart = 1
 
!     TEST CORE
 
 nchave = niz
 IF (nchave <= 0) GO TO 7001
 CALL sofopn (iz(ib1),iz(ib2),iz(ib3))
 
!     CHECK STATUS OF SUBSTRUCTURE BEING REFERENCED - NAME
 
 CALL sfetch (NAME,nlods,3,igo)
 IF (igo  == 4) GO TO 7002
 IF (idry < 0) GO TO 1001
 
!     CHECK LOCATION OF THE PAPP VECTOR - EITHER ON FILE IPAPP OR SOF
 
 iz(1) = ipapp
 CALL rdtrl (iz(1))
 IF (iz(1) > 0) GO TO 10
 lpapp = .false.
 iuapp = iscr1
 nitem = npapp
 CALL mtrxi (iuapp,NAME,nitem,0,itest)
 IF (itest /= 1) GO TO 7003
 lpapp = .true.
 GO TO 20
 10 lpapp = .true.
 iuapp = ipapp
 20 CONTINUE
 
!     CHECK STATUS OF THE POAP VECTOR
!     FIRST GET THE NAME OF THE LOWER LEVEL SUBSTRUCTURE WHERE
!     THE POAP ITEM TO BE USED IS LOCATED
 
 lpoap = .false.
 llsub = .false.
 CALL fndlvl (NAME,namell)
 IF (namell(1) == blnk) GO TO 7002
 IF (NAME(1) /= namell(1) .OR. NAME(2) /= namell(2)) llsub = .true.
 IF (.NOT.llsub) GO TO 40
 iz(1) = ipoap
 CALL rdtrl (iz(1))
 IF (iz(1) > 0) GO TO 30
 iuoap = iscr2
 nitem = npoap
 CALL mtrxi (iuoap,namell,nitem,0,itest)
 IF (itest /= 1) GO TO 40
 lpoap = .true.
 GO TO 40
 30 lpoap = .true.
 iuoap = ipoap
 40 CONTINUE
 
!     ESTABLISH TYPE OF CASE BEING RUN, I.E. THE CASE NO. BEING DEFINED
!     BY  NN(1) AND NN(2).
 
 nn(1) = 0
 nn(2) = 0
 
!     CHECK STATUS OF PVEC AND POVE VECTORS
 
!     1) PVEC
 
 lpvec = .true.
 iz(1) = 0
 CALL softrl (NAME,npvec,iz(1))
 IF (iz(1) /= 1) lpvec = .false.
 
!     2) POVE
 
 lpove = .true.
 iz(1) = 0
 IF (llsub) CALL softrl (namell,npove,iz(1))
 IF (iz(1) /= 1) lpove = .false.
 
!     KNOWING THE STATUS OF PAPP, PVEC, POAP, POVE DEFINE CASE NO.
 
 IF (lpapp) GO TO 50
 nn(1) = 4
 IF (lpvec) nn(1) = 3
 GO TO 60
 50 nn(1) = 2
 IF (lpvec) nn(1) = 1
 60 CONTINUE
 IF (lpoap) GO TO 70
 nn(2) = 4
 IF (lpove) nn(2) = 3
 GO TO 80
 70 nn(2) = 2
 IF (lpove) nn(2) = 1
 80 igo = nn(2)
 
!     KNOWING NN(1) AND NN(2) THE CASE IS DEFINED
 
 IF (nn(1) == 1) THEN
    SELECT CASE ( igo )
     CASE (    1)
       GO TO 1001
     CASE (    2)
       GO TO 7004
     CASE (    3)
       GO TO 7004
     CASE (    4)
       GO TO 1001
   END SELECT
 END IF
 IF (nn(1) == 2) THEN
    SELECT CASE ( igo )
     CASE (    1)
       GO TO 7004
     CASE (    2)
       GO TO 1001
     CASE (    3)
       GO TO 7004
     CASE (    4)
       GO TO 1001
   END SELECT
 END IF
 IF (nn(1) == 3) GO TO 7004
 IF (nn(1) == 4) GO TO 7004
 
!     READ IN LOAP DATA
 
 1001 irw = 1
 CALL sfetch (NAME,nloap,irw,itloap)
 IF (itloap > 1) GO TO 7003
 CALL suread (iz(nstart),-1,nwds,ichk)
 nl = iz(nstart+2)
 ns = iz(nstart+3)
 nfini  = 4 + ns*3 + nl
 nstart = 5 + ns*2
 nas = 1
 naf = nfini
 nchave = nfini
 nsubs  = ns
 IF (nchave > niz) GO TO 7001
 nbasn  = 4 + nsubs*2 + 1
 DO  iloop = 1,nsubs
   CALL suread (iz(nbasn),-1,nwds,ichk)
   nbasn  = iz(nbasn) + 1 + nbasn
 END DO
 nstart = naf + 1
 nps    = nstart
 lmerg  = .true.
 
!     IF DRY RUN (IDRY .LT. 0) CHECK FOR LODS ITEM
 
 IF (idry < 0) GO TO 1002
 IF (lpvec) GO TO 1002
 
!     SIMPLE COPY OF NEW APPENDED LOADS TO SOF
 
!     NEW  LODS  ITEM
 
 nitem = nlods
 lmerg = .false.
 itest = 3
 irw   = 2
 CALL sfetch (NAME,nitem,irw,itest)
 IF (itest /= 3) GO TO 7005
 nwds = 4 + 2*iz(nas+3)
 CALL suwrt (iz(nas),nwds,2)
 nbasn = nas + nwds
 DO  n = 1,nsubs
   nwds = iz(nbasn) + 1
   CALL suwrt (iz(nbasn),nwds,2)
   nbasn = nbasn + nwds
 END DO
 
!     END OF ITEM CALL TO SUWRT
 
 CALL suwrt (0,0,3)
 
!     NEW  PVEC  ITEM
 
 CALL mtrxo (iuapp,NAME,npvec,0,itest)
 
!     NEW  POVE  ITEM  IF ANY
 
 CALL bug (nprog(1),101,lpoap,1)
 IF (lpoap) CALL mtrxo (iuoap,namell,npove,0,itest)
 
!     MODULE IS FINISHED WITH THE DIRECT COPY CASE
 
 GO TO 7000
 1002 CONTINUE
 
!     ITS BEEN DETERMINED THAT A MERGE OPERATION WILL TAKE PLACE.  THE
!     ONLY CHECK NOW IS TO SEE IF A LODS ITEM EXISTS.
 
 irw  = 1
 nitem = nlods
 CALL sfetch (NAME,nitem,irw,itlods)
 IF (itlods /= 1 .AND. idry > 0) GO TO 7003
 IF (itlods /= 1 .AND. idry < 0) GO TO 9999
 CALL suread (iz(nstart),-1,nwds,ichk)
 ncnt = nwds + naf
 nl   = iz(nps+2)
 ns   = iz(nps+3)
 nfini= nps + 3 + 3*ns + nl
 npf  = nfini
 nstart = npf + 1
 nchave = nfini
 IF (nchave > niz) GO TO 7001
 nbasn = nps + 3 + 2*nsubs + 1
 DO  iloop = 1,nsubs
   CALL suread (iz(nbasn),-1,nwds,ichk)
   nbasn = nbasn + iz(nbasn) + 1
 END DO
 nldsa = iz(nas+2)
 nldsp = iz(nps+2)
 nloads= nldsa + nldsp
 nbasa = nas + 3 + 2*nsubs + 1
 nbasp = nps + 3 + 2*nsubs + 1
 nlbasa= iz(nbasa)
 nlbasp= iz(nbasp)
 
!     CHECK FOR DUPLICATE LOAD IDS IN THE  LOAP  AND  LODS  ITEMS.
 
 DO  l = 1,nsubs
   IF (nlbasp == 0 .OR. nlbasa == 0) GO TO 104
   DO  m = 1,nlbasa
     DO  n = 1,nlbasp
       IF( iz(nbasa+m) /= iz(nbasp+n).OR.iz(nbasa+m) == 0 ) CYCLE
       WRITE  (lp,6955) ufm,iz(nbasa+m),NAME
     6955 FORMAT (a23,' 6955, DUPLICATE LOAD IDS DURING APPEND OPERATION.',  &
         '  LOAD ID NO.',i9,' SUBSTRUCTURE ',2A4)
     idry = -2
   END DO
 END DO
 104 nbasa = nbasa + nlbasa + 1
 nbasp = nbasp + nlbasp + 1
 nlbasa = iz(nbasa)
 nlbasp = iz(nbasp)
END DO

!     END OF RUN IF A DRY RUN(IDRY .LT. 0)

IF (idry < 0) GO TO 9999

!     CALCULATE LENGTH OF THE MERG AND  N E W  LODS TABLE AND THEIR
!     LOCATIONS IN OPEN CORE

lmergt = 2*nsubs
nmergs = npf + 1
nmergf = npf + lmergt
lnewlt = 4 + 3*nsubs + nloads
nnews  = nmergf + 1
nnewf  = nmergf + lnewlt

!     CREATE THE NEW LODS TABLE IN OPEN CORE - GROUP  0  FIRST

iz(nnews  ) = iz(nas  )
iz(nnews+1) = iz(nas+1)
iz(nnews+2) = nloads
iz(nnews+3) = nsubs
nloop = 2*nsubs
nnew1 = nnews + 3
ndel1 = nas + 3
DO  ns1 = 1,nloop
  iz(nnew1+ns1) = iz(ndel1+ns1)
END DO

!     COMPLETION OF THE NEW LODS TABLE - GROUPS  1  THRU  NSUBS  --  AND
!     CREATION OF THE MERGE TABLE

nbasn  = nnew1 + nloop + 1
nbasa  = nas + 3 + 2*nsubs + 1
nbasp  = nps + 3 + 2*nsubs + 1
nloada = iz(nbasa)
nloadp = iz(nbasp)
nloadn = nloada + nloadp
nmergn = nmergs
imergn = 1

!     ZERO THE MERG TABLE LOCATION

DO  i = 1,lmergt
  iz(npf+i) = 0
END DO
DO  iloop = 1,nsubs
  iz(nbasn) = nloadn
  IF (nloadp == 0) GO TO 2108
  DO  n = 1,nloadp
    iz(nbasn+n) = iz(nbasp+n)
  END DO
  2108 CONTINUE
  nbasn = nbasn + nloadp
  IF (nloada == 0) GO TO 2107
  DO  n = 1,nloada
    iz(nbasn+n) = iz(nbasa+n)
  END DO
  2107 CONTINUE
  
!     LOCATION IN THE MERGE TABLE OF THE  1(S)
  
  imergn = imergn + nloadp
  iz(nmergn) = imergn
  imergn = imergn + nloada
  iz(nmergn+1) = nloada
  nmergn = nmergn + 2
  nbasn  = nbasn + nloada + 1
  IF (iloop == nsubs) CYCLE
  nbasa  = nbasa + nloada + 1
  nbasp  = nbasp + nloadp + 1
  nloada = iz(nbasa)
  nloadp = iz(nbasp)
  nloadn = nloada + nloadp
END DO

!     END OF GENERATION OF NEW LODS ITEM AND CREATION OF MERGE TABLE

!     CALCULATE BEGINNING LOCATION OF MERGE VECTOR IN OPEN CORE

nmrvcs = nnewf + 1
nmergn = nmergs
nmrvcn = nmrvcs - 2
lvect  = iz(nnews+2)
nmrvcf = nnewf + lvect
nchave = nmrvcf
IF (nchave > niz) GO TO 7001

!     FILL THE MERGE VECTOR WITH  1(S)  ACCORDING TO THE MERGE TABLE

!     1) ZERO FIRST

DO  i = 1,lvect
  rz(nmrvcs-1+i) = 0.
END DO

!     2) NOW FILL

DO  iloop = 1,nsubs
  idrc1 = iz(nmergn)
  idrc2 = iz(nmergn+1)
  IF (idrc2 == 0) GO TO 2111
  DO  n = 1,idrc2
    rz(nmrvcn+idrc1+n) = 1.0
  END DO
  2111 CONTINUE
  nmergn = nmergn + 2
END DO

!     WRITE THE MERGE VECTOR ON SCRATCH  5  USING  PACK-SEE COMMON PACKX
!     THIS IS A COLUMN PARTITIONING VECTOR (REFERRED TO AS A ROW VECTOR
!     BY MERGE)

itypin = 1
itypot = 1
ifirst = 1
ilast  = lvect
incr   = 1
CALL gopen (iscr5,iz(ibuf1),iwrtrw)

!     ZERO THE TRAILER INFO. LOCATIONS

DO  i = 1,7
  mcbloc(i) = 0
END DO
mcbloc(1) = iscr5
mcbloc(3) = lvect
mcbloc(4) = irect
mcbloc(5) = irsp
CALL pack (rz(nmrvcs),iscr5,mcbloc(1))
CALL CLOSE (iscr5,irew)
CALL wrttrl (mcbloc(1))
idump = -iscr5
CALL dmpfil (idump,adump,4000)

!     READ IN THE  PVEC  AND  POVE(IF EXISTS)  USING MTRXI

!     1) PVEC

iuvec = iscr3
CALL mtrxi (iuvec,NAME,npvec,0,ichk)

!     2) POVE

iuove = iscr4
IF (lpove) CALL mtrxi (iuove,namell,npove,0,ichk)

!     SET UP TO CALL MERGE FOR  PAPP  AND  PVEC

idump = -iuvec
CALL dmpfil (idump,adump,4000)
idump = -iuapp
CALL dmpfil (idump,adump,4000)
icore = nmrvcf + 1
iz(icore) = iscr5
CALL rdtrl (iz(icore))

!     SETUP NULL ROW PARTITIONING VECTOR USING  ISCR8
!     THIS IS A ROW PARTITIONING VECTOR REFERRED TO AS A COLUMN VECTOR
!     BY MERGE)

mcbk11(1) = iuvec
CALL rdtrl (mcbk11(1))
mcbk12(1) = iuapp
CALL rdtrl (mcbk12(1))
DO  k = 1,7
  mcbk21(k) = 0
  mcbk22(k) = 0
END DO
iz(icore+ 7) = iscr8
iz(icore+ 8) = 0
iz(icore+ 9) = mcbk11(3)
iz(icore+10) = irect
iz(icore+11) = irsp
iz(icore+12) = 0
iz(icore+13) = 0
ncnt = icore +13

CALL gopen (iscr8,iz(ibuf1),iwrtrw)
itypin = 1
itypot = 1
ifirst = 1
ilast  = 1
incr   = 1
CALL pack (0,iscr8,iz(icore+7))
CALL CLOSE (iscr8,irew)
CALL wrttrl (iz(icore+7))
lcore = ib3 - icore - 15
i = icore + 15
irule   = 0
mcbk(1) = iscr6
mcbk(2) = mcbk11(2)+mcbk12(2)
mcbk(3) = mcbk11(3)
mcbk(4) = irect
mcbk(5) = mcbk11(5)
mcbk(6) = 0
mcbk(7) = 0
CALL merge (iz(icore),iz(icore+7),iz(i))
CALL wrttrl (mcbk(1))

!     SETUP TO MERGE  POVE  AND  POAP(IF THEY EXIST)

idump = -iuove
CALL dmpfil (idump,adump,4000)
idump = -iuoap
CALL dmpfil (idump,adump,4000)
IF (.NOT. lpove) GO TO 1005
mcbk11(1) = iuove
CALL rdtrl (mcbk11(1))
mcbk12(1) = iuoap
CALL rdtrl (mcbk12(1))
DO  k = 1,7
  mcbk21(k) = 0
  mcbk22(k) = 0
END DO
irule   = 0
mcbk(1) = iscr7
mcbk(2) = mcbk11(2)+mcbk12(2)
mcbk(3) = mcbk11(3)
mcbk(4) = irect
mcbk(5) = mcbk11(5)
mcbk(6) = 0
mcbk(7) = 0
CALL merge (iz(icore),iz(icore+7),iz(i))
CALL wrttrl (mcbk(1))

!     CHECK TIME REMAINING AND RETURN WITH USER FATAL MESSAGE IF NOT
!     ENOUGH REMAINING

1005 CALL tmtogo (itime2)
IF (itime2 <= 0) GO TO 7007

!     CALCULATE TIME USED IN REACHING THIS LOCATION

itused = itime1 - itime2

!     CONTINUE IF ITIME2 IS GREATER THAN ITUSED

IF (itime2 < itused) GO TO 7007

!     WRITE NEW LODS ITEM TO SOF

!     1) DELETE OLD LODS ITEM

CALL DELETE (NAME,nlods,ichk)

!     DELETE LODS ITEMS ON ANY SUBSTRUCTURE SECONDARY TO NAME - THIS
!     WILL ALLOW THE NEW LODS ITEM TO BE COPIED DURING FUTURE EQUIV
!     OPERATIONS

ii = itcode (nlods)
CALL fdsub (NAME,ind)
CALL fmdi (ind,imdi)
ips = andf(icorx(imdi+1),1023)
IF (ips /= 0) GO TO 118
117 iss = andf(rshift(icorx(imdi+1),10),1023)
IF (iss == 0) GO TO 118
CALL fmdi (iss,imdi)
iblk = andf(icorx(imdi+ii),65535)
IF (iblk /= 0 .AND. iblk /= 65535) CALL retblk (iblk)
icorx(imdi+ii) = 0
mdiup = .true.
GO TO 117

!     2) BEGIN WRITING

118 ichk = 3
irw  = 2
CALL sfetch (NAME,nlods,irw,ichk)
nwords = 4 + 2*iz(nnews+3)
CALL suwrt (iz(nnews),nwords,2)
nbasn = nnews + nwords
DO  n = 1,nsubs
  nwords = iz(nbasn) + 1
  CALL suwrt (iz(nbasn),nwords,2)
  nbasn = nbasn + nwords
END DO
CALL suwrt (0,0,3)

!     WRITE NEW  PVEC  AND  POVE(IF IT EXISTS)  TO SOF

CALL DELETE (NAME,npvec,ichk)
CALL mtrxo (iscr6,NAME,npvec,0,ichk)
IF (lpove) CALL DELETE (namell,npove,ichk)
IF (lpove) CALL mtrxo (iscr7,namell,npove,0,ichk)

WRITE  (lp,6900) uim,NAME
6900 FORMAT (a29,' 6900, LOADS HAVE BEEN SUCCESSFULLY APPENDED FOR ',  &
    'SUBSTRUCTURE ',2A4)
GO TO 9999
7000 WRITE  (lp,6901) uim,NAME
6901 FORMAT (a29,' 6901, ADDITIONAL LOADS HAVE BEEN SUCCESSFULLY ',  &
    'MERGED FOR SUBSTRUCTURE ',2A4)
GO TO 9999
7001 WRITE  (lp,6951) ufm,nchave
6951 FORMAT (a23,' 6951, INSUFFICIENT CORE TO LOAD TABLES', /5X,  &
    'IN MODULE LODAPP, CORE =',i8)
CALL mesage (-8,nprog,0)

7002 WRITE  (lp,6952) sfm,NAME
6952 FORMAT (a25,' 6952, REQUESTED SUBSTRUCTURE ',2A4, ' DOES NOT EXIST')
idry = -2
GO TO 9999
7003 WRITE  (lp,6101) sfm,nitem,NAME
6101 FORMAT (a25,' 6101, REQUESTED SOF ITEM DOES NOT EXIST.  ITEM ',a4,  &
    ' SUBSTRUCTURE ',2A4)
idry = -2
GO TO 9999
7004 WRITE  (lp,6953) sfm,NAME
6953 FORMAT (a25,' 6953, A WRONG COMBINATION OF LOAD VECTORS EXISTS ',  &
    'FOR SUBSTRUCTURE ',2A4)
idry = -2
GO TO 9999
7005 WRITE  (lp,6954) sfm,nitem,NAME
6954 FORMAT (a25,' 6954, THE ,A4,62H ITEM EXISTS BUT HAS NO ',  &
    'ASSOCIATED PVEC ITEM FOR SUBSTRUCTURE ',2A4)
idry = -2
GO TO 9999
7007 WRITE  (lp,6956) ufm,itime2
6956 FORMAT (a23,' 6956, INSUFFICIENT TIME REMAINING FOR MODULE ',  &
    'LODAPP, TIME LEFT =',i8)
idry = -2
9999 CONTINUE
CALL sofcls

!     RETURN VALUE OF DRY PARAMETER

buf(3) = idry
RETURN
END SUBROUTINE lodapp
