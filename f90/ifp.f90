SUBROUTINE ifp
     
 IMPLICIT INTEGER (a-z)
 EXTERNAL        orf,andf
 LOGICAL :: badfor,baddat,abort,eofflg,cf,cl,iax,lharm
 INTEGER :: fnm(2,16),ii(16),ifle(16),STATUS(16),iblkda(2),  &
     jr(20),iend(3),itrl(7),ifpna1(2),nm(2),kap(4),  &
     ap(12),kkl(40),ifpna2(2),ooo(5),mentry(40),  &
     NAME(2),inam(2),itype(2),nentry(80),itype1(2), itype2(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / endara(40)
 COMMON /system/ n1,nout,abort,n2(17),iapp,n3(5),axiccc,junk(10),  &
     axifcc,dum(30),isubs
 COMMON /two   / two(32)
 COMMON /ifpdta/ id,n,k,kx,ky,i(100),m(100),mf(100),m1(100),  &
     m1f(100),kn,baddat,badfor,nopen,nparam,iax,nax,  &
     iaxf,naxf,lharm,knt,kslot1,kslot2,kslot3,kslot4, kslot5,gc(7),ll(6)
 COMMON /zzzzzz/ ibuff(1)
 
!     NCDS    = LENGTH OF T1
!     NCDSMX  = NO. OF CARD NAMES IN NASTRAN
!     T3(1,K) = THE GINO OUTPUT FILE NUMBER
!     T3(2,K) = THE APPROACH ACCEPTANCE FLAG
!     T4(1,K) = THE CONICAL SHELL PROBLEM FLAG
!     T4(2,K) = USED AS INTERNAL STORAGE WITHIN IFP
!     T5(1,K) = THE MIN NO. OF WORDS ALLOWED PER CARD
!               (MINUS MEANS OPEN-ENDED CARD)
!     T5(2,K) = THE MAX NO. OF WORDS ALLOWED PER CARD
!     T6(1,K) = THE FORMAT CHECK POINTER INTO F( )
!     T6(2,K) = FIELD 2 UNIQUENESS CHECK FLAG
!     T7(1,K) = LOCATE CODE
!     T7(2,K) = TRAILER BIT POSITION
!     F(T6(1,K)) = THE START OF THE FORMAT ACCEPTANCE STRING
 
!     T1(1,K),T1(2,K) = THE BCD CARD NAMES
 
 COMMON /ifpx0 / lbd,lcc,ib(18)
 COMMON /ifpx1 / ncds,t1(2,370)
 COMMON /ifpx2 / t3(2,370)
 COMMON /ifpx3 / t4(2,370)
 COMMON /ifpx4 / t5(2,370)
 COMMON /ifpx5 / t6(2,370)
 COMMON /ifpx6 / t7(2,370)
 COMMON /ifpx7 / f(1469)
 EQUIVALENCE     (n3(3),iumfed), (n2(9),line)
 DATA    ncdsmx/ 359 /
 DATA    nfls  / 16  /
 DATA    fnm(1, 1),fnm(2, 1) / 4HGEOM,4H1   /
 DATA    fnm(1, 2),fnm(2, 2) / 4HEPT ,4H    /
 DATA    fnm(1, 3),fnm(2, 3) / 4HMPT ,4H    /
 DATA    fnm(1, 4),fnm(2, 4) / 4HEDT ,4H    /
 DATA    fnm(1, 5),fnm(2, 5) / 4HDIT ,4H    /
 DATA    fnm(1, 6),fnm(2, 6) / 4HPVT ,4H    /
 DATA    fnm(1, 7),fnm(2, 7) / 4HDYNA,4HMICS/
 DATA    fnm(1, 8),fnm(2, 8) / 4HGEOM,4H2   /
 DATA    fnm(1, 9),fnm(2, 9) / 4HGEOM,4H3   /
 DATA    fnm(1,10),fnm(2,10) / 4HGEOM,4H4   /
 DATA    fnm(1,11),fnm(2,11) / 4HGEOM,4H5   /
 DATA    fnm(1,12),fnm(2,12) / 4HPOOL,4H    /
 DATA    fnm(1,13),fnm(2,13) / 4HFORC,4HE   /
 DATA    fnm(1,14),fnm(2,14) / 4HMATP,4HOOL /
 DATA    fnm(1,15),fnm(2,15) / 4HAXIC,4H    /
 DATA    fnm(1,16),fnm(2,16) / 4HIFPF,4HILE /
 DATA    ifle   / 201,202,203,204,205,4HNPTP,207,208,209,210,211,  &
     4HPOOL ,213,214,215,216   /
 DATA    iend   , eofz   /3*2147483647,4HZZZZ    /
 DATA    kkl    / 48, 49, 50, 67, 71, 75, 68, 72, 76, 11, 10*0  ,  &
     45, 46, 44, 41,250,260, 39, 42,121, 34, 37, 43, 31, 7*0/
 DATA    iblkda / 4HBULK, 4HDATA /,  ooo    / 1HA,1HB,1HC,1HD,1HE/
 DATA    BLANK  / 1H   /, kap    / 0,-1,1,-1/
DATA    ifpna1 / 4HIFP ,4HBEGN/, ifpna2/ 4HIFP ,4HEND /
DATA    iparm  , ivary /4H1PAR , 4H1VAR/
DATA    icount , jcount, kcount/ 3*0   /
DATA    it1k   , it2k  , jt1k  ,jt2k   , kt1k  ,kt2k /  &
    1H , 1H , 1H , 1H , 1H , 1H   /
DATA    ap     / 4HDMAP,4H     , 4H    , 4HDISP,4HLACE , 4HMENT,  &
    4HHEAT,4H     , 4H    , 4HAERO,4H     , 4H    /
DATA    mentry / 3001  , 3701  ,  3901 ,  1201 ,   401,  &
    801  , 1301  ,   501 ,   901 ,  5201, 10*0,  &
    202  ,  302  ,   402 ,   502 ,  2202,  &
    5302  ,  802  ,  1002 ,  2102 ,  1302, 1402  , 1702  ,  1802 ,   7*0 /
DATA    NAME   / 4HIFP , 4H    /
DATA    nentry / 4HCROD, 4H    , 4HCTUB, 4HE   , 4HCVIS, 4HC   ,  &
    4HCMAS, 4HS3  , 4HCDAM, 4HP3  , 4HCELA, 4HS3  ,  &
    4HCMAS, 4HS4  , 4HCDAM, 4HP4  , 4HCELA, 4HS4  , 4HPLOT, 4HEL  , 20*0  ,  &
    4HPDAM, 4HP   , 4HPELA, 4HS   , 4HPMAS, 4HS   ,  &
    4HPQDM, 4HEM  , 4HPQDM, 4HEM1 , 4HPQDM, 4HEM2 ,  &
    4HPQUA, 4HD2  , 4HPSHE, 4HAR  , 4HPTOR, 4HDRG ,  &
    4HPTRI, 4HA2  , 4HPTRM, 4HEM  , 4HPTWI, 4HST  , 4HPVIS, 4HC   , 14*0  /
DATA    itype1 / 4HELEM, 4HENT /
DATA    itype2 / 4HPROP, 4HERTY/

!     ============================================================
!     REMEMBER TO CHECK FOR THE LONGEST LINK IN OVERLAY STRUCTURE.
!     ============================================================

!     INITIALIZE COMMON BLOCKS CIFS1P, 2P, 3P, 4P, AND CIFS5P

CALL cifsdd

DO  j = 1,16
  STATUS(j) = 1
END DO
STATUS( 6) = 3
STATUS(12) = 3
lm     = 100
curfil = 0
kick   = 0
ipvs   = 0
eofflg = .false.
baddat = .false.
badfor = .false.
nparam = 0
kn     = 0
iax    = .false.
nax    =-1
iaxf   = 0
naxf   =-1
lharm  = .true.
kslot1 = 0
kslot2 = 0
kslot3 = 0
kslot4 = 0
kslot5 = 0
CALL conmsg (ifpna1,2,0)
iap    = IABS(iapp)
jap    = kap(iap)
knt    =-1
iaxic  = axiccc
iaxif  = axifcc
axiccc = 0
axifcc = 0
DO  j = 1,nfls
  ii(j)  = 0
END DO
DO  j = 1,40
  endara(j) = 0
END DO
nopen = korsz(ibuff) - 3*n1
CALL sswtch (42,l42)
IF (nopen >= 0) GO TO 100
CALL page2 (2)
WRITE  (nout,40) sfm
40 FORMAT (a25,' 303, NO OPEN CORE FOR IFP.')
abort =.true.
RETURN

!     OPEN NPTP AND LOCATE BULK DATA

100 kfil = ifle(6)
CALL OPEN (*130,kfil,ibuff(n1+1),0)
110 CALL skpfil (kfil,1)
CALL READ (*1390,*160,kfil,jr,2,1,kdum)
IF (jr(1) == iblkda(1) .AND. jr(2) == iblkda(2)) GO TO 180
kick = kick + 1
IF (kick < 5) GO TO 110
CALL page2 (2)
WRITE  (nout,120) sfm,jr(1),jr(2)
120 FORMAT (a25,' 304, IFP NOT READING NPTP. FILE BEING READ = ',2A4)
GO TO 150
130 CALL page2 (2)
WRITE  (nout,140) sfm,kfil
140 FORMAT (a25,' 305, IFP CANNOT OPEN GINO FILE',i10)
150 abort =.true.
GO TO 1850
160 CALL page2 (2)
WRITE  (nout,170) sfm
170 FORMAT (a25,' 306, READ LOGICAL RECORD ERROR')
GO TO 150
180 CALL READ (*1380,*160,ifle(6),jr,20,1,kdum)
knt = knt + 1

!     CHECK FOR 1PARM OR 1VARY CARDS

IF (jr(1) == iparm .OR. jr(1) == ivary) CALL ifppvc (*190,ipvs,jr)
IF (l42 == 0) CALL rcard2 (m1,m1f,nw,jr)
IF (l42 /= 0) CALL rcard  (m1,m1f,nw,jr)
IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1430
GO TO 220
190 CALL CLOSE (ifle(6),1)
GO TO 1900

!     READ AND DECODE ONE PHYSICAL CARD

200 IF (eofflg) GO TO 1410
CALL READ (*1460,*160,ifle(6),jr,20,1,kdum)
knt = knt + 1
IF (l42 == 0) CALL rcard2 (m1,m1f,nw,jr)
IF (l42 /= 0) CALL rcard  (m1,m1f,nw,jr)
IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 200
210 IF (eofflg) GO TO 1460

!     IDENTIFY CARD NAME

220 DO  j = 1,ncdsmx
  k = j
  IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k)) GO TO 280
END DO
IF (kt1k /= t1(1,k) .OR. kt2k /= t1(2,k)) GO TO 240
kcount = kcount + 1
IF (kcount-7 < 0) THEN
  GO TO   250
ELSE IF (kcount-7 == 0) THEN
  GO TO   270
ELSE
  GO TO   200
END IF
240 kt1k = t1(1,k)
kt2k = t1(2,k)
250 CALL page2 (2)
WRITE  (nout,260) ufm,m1(1),m1(2)
260 FORMAT (a23,' 307, ILLEGAL NAME FOR BULK DATA CARD ',2A4 )
abort =.true.
GO TO 200
270 CALL page2 (3)
WRITE  (nout,1150)
GO TO 200
280 kcount = 0
cl =.false.
cf =.true.
kx = k  - 100
ky = kx - 100

!     CHECK APPROACH ACCEPTABILITY

IF (t3(2,k)*jap+1 < 0.0) THEN
  GO TO   300
ELSE IF (t3(2,k)*jap+1 == 0.0) THEN
  GO TO   320
ELSE
  GO TO   340
END IF
300 WRITE  (nout,310) ufm,t1(1,k),t1(2,k),ap(3*iap-2),ap(3*iap-1), ap(3*iap)
310 FORMAT (a23,' 308, CARD ',2A4,' NOT ALLOWED IN ',3A4,' APPROACH.')
CALL page2 (2)
abort =.true.
GO TO 340
320 WRITE  (nout,330) uwm,t1(1,k),t1(2,k),ap(3*iap-2),ap(3*iap-1), ap(3*iap)
330 FORMAT (a25,' 309, CARD ',2A4,' IMPROPER IN ',3A4,' APPROACH.')
CALL page2(2)
340 IF (.NOT.iax .OR. t4(1,k) >= 0) GO TO 400
CALL page2 (2)
WRITE  (nout,350) ufm,t1(1,k),t1(2,k)
350 FORMAT (a23,' 310, CARD ',2A4,' NOT ALLOWED IN SAME DECK WITH ',  &
    'AXIC CARD.')
abort =.true.

!     ESTABLISH PROPER OUTPUT FILES FOR THIS CARD

400 indx = t3(1,k)
IF (indx == curfil .OR. indx == 6) GO TO 420
IF (curfil == 0 .OR. STATUS(curfil) == 1) GO TO 410
CALL CLOSE (ifle(curfil),2)
STATUS(curfil) = 3
410 kfil = ifle(indx)
CALL OPEN (*130,kfil,ibuff,STATUS(indx))
curfil = indx
STATUS(curfil) = -STATUS(curfil)
IF (STATUS(curfil) /= -1) GO TO 420
CALL WRITE (ifle(curfil),fnm(1,curfil),2,1)
ii(curfil) = 1
STATUS(curfil) = -3
420 id = m1(3)
430 jf = nw - 2
DO  l = jf,lm
  m(l) = 0
END DO
DO  l = 1,jf
  m(l) = m1(l+2)
END DO

!     TEST UNIQUENESS OF FIELD 2 IF APPLICABLE

IF (m1(1) == 0 .AND. m1(2) == 0 .OR.cf .OR. t6(2,k) /= 1) GO TO 480
IF (id == m(1)) GO TO 460
id = m(1)
GO TO 480
460 knt1 = knt + 1
CALL page2 (2)
WRITE  (nout,470) ufm,t1(1,k),t1(2,k),m(1),knt1
470 FORMAT (a23,' 311, NON-UNIQUE FIELD 2 ON BULK DATA CARD ',2A4,i8,  &
    10X,'H SORTED CARD COUNT =',i7)
abort =.true.
480 DO  l = 1,lm
  mf(l) = 0
END DO
lf = 0
DO  l = 1,jf
  
!     =========================================
!     THIS SHOULD BE CHANGED WHEN RCARD CHANGES
  
  IF (m1f(l+1) < 0) GO TO 540
!     ========================================
  lf = lf + 1
  mf(l) = m1f(l+1)
END DO
GO TO 540
510 IF (eofflg) GO TO 1420

!     READ ANOTHER CARD (TO BE PROCESSED NEXT)

knt = knt + 1
CALL READ (*550,*160,ifle(6),jr,20,1,kdum)
IF (l42 == 0) CALL rcard2 (m1,m1f,nw,jr)
IF (l42 /= 0) CALL rcard  (m1,m1f,nw,jr)
IF (m1(1) /= 0  .OR. m1(2) /= 0) GO TO 580

!     CHECK FOR TOO MANY CONTINUATIONS

IF (t6(1,k) < 0 .AND. lf > 4) GO TO 600
IF (jf+nw-2-lm > 0) GO TO 560
k1 = nw - 2
DO  l = 1,k1
  k2 = jf + l
  m(k2) = m1(l+2)
END DO
jf = jf + nw - 2
DO  l = 1,k1
  
!     =========================================
!     THIS SHOULD BE CHANGED WHEN RCARD CHANGES
  
  IF (m1f(l+1) < 0) GO TO 540
!     =========================================
  lf = lf + 1
  mf(lf) = m1f(l+1)
END DO
540 mf(lf+1) = -32767
GO TO 510
550 eofflg =.true.
m1(1)  = eofz
m1(2)  = eofz
GO TO 590
560 WRITE  (nout,570) ufm,t1(1,k),t1(2,k),m(1),knt
570 FORMAT (a23,' 312, TOO MANY CONTINUATIONS FOR BULK DATA CARD ',  &
    2A4,i8,6X,'SORTED CARD COUNT =',i7)
CALL page2 (2)
abort =.true.
GO TO 510
580 IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k)) GO TO 600
590 cl =.true.
600 IF (.NOT.cf .OR. t6(2,k) == 2) GO TO 610
kkk = t3(1,k)
ii(kkk) = ii(kkk) + 1
cf =.false.
IF (kkk == 6 .OR. kkk == 12) GO TO 640
itrl(1) = t7(1,k)
itrl(2) = t7(2,k)
itrl(3) = k
CALL WRITE (ifle(curfil),itrl,3,0)

!     CHECK FOR MIN-MAX NO. OF WORDS

610 IF (t5(1,k) < 0) GO TO 640
l = jf
IF (t5(1,k)-l < 0.0) THEN
  GO TO   630
ELSE IF (t5(1,k)-l == 0.0) THEN
  GO TO   690
ELSE
  GO TO   650
END IF
620 l = l + 4
630 IF (t5(2,k)-l < 0.0) THEN
  GO TO   650
ELSE IF (t5(2,k)-l == 0.0) THEN
  GO TO   690
ELSE
  GO TO   620
END IF
640 l =-t5(1,k)
IF (jf >= l .AND. jf <= t5(2,k)) GO TO 690
650 WRITE  (nout,660) ufm,t1(1,k),t1(2,k),m(1),knt
660 FORMAT (a23,' 313, ILLEGAL NUMBER OF WORDS ON BULK DATA CARD ',  &
    2A4,i8,6X,'SORTED CARD COUNT =',i7)
WRITE  (nout,670) t5(1,k),t5(2,k),k,l,jf
670 FORMAT ('   T5(1&2,K),K,L,JF =',5I4)
CALL page2 (2)
abort =.true.
IF (t6(1,k) < 0.0) THEN
  GO TO   710
END IF
680 IF (.NOT.cl) GO TO 430
IF (t6(2,k) == 2) GO TO 210
CALL WRITE (ifle(curfil),m,0,1)
IF (t4(2,k) > 0) GO TO 210
ii(kkk) = ii(kkk) - 1
CALL bckrec (ifle(curfil))
GO TO 210

!     CHECK FOR PROPER FORMAT

690 IF (t6(1,k) < 0) GO TO 710
l  = t6(1,k)
l1 = 0
DO  k1 = 1,lf
  l1 = l1 + 1
  IF (mf(k1) == 3) l1 = l1 + 1
  k2 = l + k1 - 1
  IF (f(k2) == mf(k1) .OR. f(k2) == 5) CYCLE
  IF (mf(k1) == 1    .AND. m(l1) == 0) CYCLE
  IF (mf(k1) /= 0 .OR. f(k2) /= 1 .AND. f(k2) /= 2) GO TO 1350
END DO
710 n = 0
baddat =.false.
badfor =.false.
IF (ipvs /= 0) CALL ifpmdc

!     CALL SECONDARY ROUTINE TO EXAMINE EACH TYPE OF CARD

kb = (k-1)/20 + 1
IF (kb > 18) GO TO 1060
SELECT CASE ( kb )
  CASE (    1)
    GO TO  810
  CASE (    2)
    GO TO  820
  CASE (    3)
    GO TO  830
  CASE (    4)
    GO TO  840
  CASE (    5)
    GO TO  850
  CASE (    6)
    GO TO  860
  CASE (    7)
    GO TO  870
  CASE (    8)
    GO TO  880
  CASE (    9)
    GO TO  890
  CASE (   10)
    GO TO  900
  CASE (   11)
    GO TO 910
  CASE (   12)
    GO TO  920
  CASE (   13)
    GO TO  930
  CASE (   14)
    GO TO  940
  CASE (   15)
    GO TO  950
  CASE (   16)
    GO TO  960
  CASE (   17)
    GO TO  970
  CASE (   18)
    GO TO  980
END SELECT
810 kb = k
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1030
  CASE (    2)
    GO TO 1030
  CASE (    3)
    GO TO 1050
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1030
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1030
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
820 kb = k - 20
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1010
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1030
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1010
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
830 kb = k - 40
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1010
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1010
  CASE (   13)
    GO TO 1010
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
840 kb = k - 60
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1010
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1010
  CASE (   13)
    GO TO 1010
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1040
  CASE (   20)
    GO TO 1040
END SELECT
850 kb = k - 80
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1030
  CASE (    3)
    GO TO 1030
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1020
  CASE (    6)
    GO TO 1020
  CASE (    7)
    GO TO 1020
  CASE (    8)
    GO TO 1050
  CASE (    9)
    GO TO 1020
  CASE (   10)
    GO TO 1040
  CASE (   11)
    GO TO 1040
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1020
  CASE (   14)
    GO TO 1020
  CASE (   15)
    GO TO 1020
  CASE (   16)
    GO TO 1020
  CASE (   17)
    GO TO 1020
  CASE (   18)
    GO TO 1040
  CASE (   19)
    GO TO 1050
  CASE (   20)
    GO TO 1050
END SELECT
860 kb = k - 100
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1050
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1050
  CASE (    4)
    GO TO 1040
  CASE (    5)
    GO TO 1040
  CASE (    6)
    GO TO 1050
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1050
  CASE (    9)
    GO TO 1050
  CASE (   10)
    GO TO 1050
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1050
  CASE (   14)
    GO TO 1050
  CASE (   15)
    GO TO 1050
  CASE (   16)
    GO TO 1050
  CASE (   17)
    GO TO 1050
  CASE (   18)
    GO TO 1050
  CASE (   19)
    GO TO 1020
  CASE (   20)
    GO TO 1020
END SELECT
870 kb = k - 120
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1030
  CASE (    4)
    GO TO 1040
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1030
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1030
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1020
  CASE (   14)
    GO TO 1020
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1030
  CASE (   19)
    GO TO 1030
  CASE (   20)
    GO TO 1020
END SELECT
880 kb = k - 140
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1020
  CASE (    2)
    GO TO 1010
  CASE (    3)
    GO TO 1030
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1030
  CASE (    6)
    GO TO 1030
  CASE (    7)
    GO TO 1030
  CASE (    8)
    GO TO 1030
  CASE (    9)
    GO TO 1030
  CASE (   10)
    GO TO 1030
  CASE (   11)
    GO TO 1030
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1030
  CASE (   14)
    GO TO 1030
  CASE (   15)
    GO TO 1030
  CASE (   16)
    GO TO 1030
  CASE (   17)
    GO TO 1030
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1050
  CASE (   20)
    GO TO 1050
END SELECT
890 kb = k - 160
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1050
  CASE (    2)
    GO TO 1020
  CASE (    3)
    GO TO 1050
  CASE (    4)
    GO TO 1050
  CASE (    5)
    GO TO 1050
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1050
  CASE (    9)
    GO TO 1050
  CASE (   10)
    GO TO 1050
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1050
  CASE (   14)
    GO TO 1050
  CASE (   15)
    GO TO 1050
  CASE (   16)
    GO TO 1050
  CASE (   17)
    GO TO 1050
  CASE (   18)
    GO TO 1050
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
900 kb = k - 180
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1030
  CASE (    3)
    GO TO 1030
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1030
  CASE (    6)
    GO TO 1050
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1020
  CASE (    9)
    GO TO 1040
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1020
  CASE (   12)
    GO TO 1020
  CASE (   13)
    GO TO 1050
  CASE (   14)
    GO TO 1050
  CASE (   15)
    GO TO 1040
  CASE (   16)
    GO TO 1040
  CASE (   17)
    GO TO 1060
  CASE (   18)
    GO TO 1050
  CASE (   19)
    GO TO 1040
  CASE (   20)
    GO TO 1020
END SELECT
910 kb = k - 200
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1040
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1040
  CASE (    4)
    GO TO 1040
  CASE (    5)
    GO TO 1040
  CASE (    6)
    GO TO 1040
  CASE (    7)
    GO TO 1040
  CASE (    8)
    GO TO 1040
  CASE (    9)
    GO TO 1040
  CASE (   10)
    GO TO 1040
  CASE (   11)
    GO TO 1040
  CASE (   12)
    GO TO 1040
  CASE (   13)
    GO TO 1040
  CASE (   14)
    GO TO 1040
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1030
  CASE (   17)
    GO TO 1040
  CASE (   18)
    GO TO 1040
  CASE (   19)
    GO TO 1040
  CASE (   20)
    GO TO 1040
END SELECT
920 kb = k - 220
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1040
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1010
  CASE (    6)
    GO TO 1010
  CASE (    7)
    GO TO 1010
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1010
  CASE (   13)
    GO TO 1010
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1040
  CASE (   20)
    GO TO 1010
END SELECT
930 kb = k - 240
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1050
  CASE (    6)
    GO TO 1050
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1050
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1010
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1050
  CASE (   14)
    GO TO 1050
  CASE (   15)
    GO TO 1050
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
940 kb = k - 260
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1020
  CASE (    2)
    GO TO 1020
  CASE (    3)
    GO TO 1050
  CASE (    4)
    GO TO 1050
  CASE (    5)
    GO TO 1050
  CASE (    6)
    GO TO 1050
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1010
  CASE (    9)
    GO TO 1050
  CASE (   10)
    GO TO 1050
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1030
  CASE (   14)
    GO TO 1030
  CASE (   15)
    GO TO 1050
  CASE (   16)
    GO TO 1050
  CASE (   17)
    GO TO 1050
  CASE (   18)
    GO TO 1050
  CASE (   19)
    GO TO 1030
  CASE (   20)
    GO TO 1020
END SELECT
950 kb = k - 280
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1020
  CASE (    2)
    GO TO 1020
  CASE (    3)
    GO TO 1020
  CASE (    4)
    GO TO 1030
  CASE (    5)
    GO TO 1030
  CASE (    6)
    GO TO 1030
  CASE (    7)
    GO TO 1030
  CASE (    8)
    GO TO 1030
  CASE (    9)
    GO TO 1010
  CASE (   10)
    GO TO 1030
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1010
  CASE (   13)
    GO TO 1010
  CASE (   14)
    GO TO 1010
  CASE (   15)
    GO TO 1040
  CASE (   16)
    GO TO 1040
  CASE (   17)
    GO TO 1030
  CASE (   18)
    GO TO 1030
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
960 kb = k - 300
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1050
  CASE (    2)
    GO TO 1050
  CASE (    3)
    GO TO 1050
  CASE (    4)
    GO TO 1050
  CASE (    5)
    GO TO 1050
  CASE (    6)
    GO TO 1050
  CASE (    7)
    GO TO 1050
  CASE (    8)
    GO TO 1050
  CASE (    9)
    GO TO 1050
  CASE (   10)
    GO TO 1050
  CASE (   11)
    GO TO 1050
  CASE (   12)
    GO TO 1050
  CASE (   13)
    GO TO 1050
  CASE (   14)
    GO TO 1050
  CASE (   15)
    GO TO 1010
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1010
END SELECT
970 kb = k - 320
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1040
  CASE (    2)
    GO TO 1040
  CASE (    3)
    GO TO 1040
  CASE (    4)
    GO TO 1040
  CASE (    5)
    GO TO 1040
  CASE (    6)
    GO TO 1040
  CASE (    7)
    GO TO 1040
  CASE (    8)
    GO TO 1040
  CASE (    9)
    GO TO 1030
  CASE (   10)
    GO TO 1030
  CASE (   11)
    GO TO 1010
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1040
  CASE (   14)
    GO TO 1040
  CASE (   15)
    GO TO 1040
  CASE (   16)
    GO TO 1040
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1050
  CASE (   19)
    GO TO 1050
  CASE (   20)
    GO TO 1010
END SELECT
980 kb = k - 340
SELECT CASE ( kb )
  CASE (    1)
    GO TO 1010
  CASE (    2)
    GO TO 1010
  CASE (    3)
    GO TO 1010
  CASE (    4)
    GO TO 1010
  CASE (    5)
    GO TO 1030
  CASE (    6)
    GO TO 1030
  CASE (    7)
    GO TO 1030
  CASE (    8)
    GO TO 1030
  CASE (    9)
    GO TO 1030
  CASE (   10)
    GO TO 1030
  CASE (   11)
    GO TO 1030
  CASE (   12)
    GO TO 1030
  CASE (   13)
    GO TO 1030
  CASE (   14)
    GO TO 1020
  CASE (   15)
    GO TO 1060
  CASE (   16)
    GO TO 1010
  CASE (   17)
    GO TO 1010
  CASE (   18)
    GO TO 1010
  CASE (   19)
    GO TO 1010
  CASE (   20)
    GO TO 1060
END SELECT
1010 CALL ifs1p (*1360,*680,*1100)
GO TO 1230
1020 CALL ifs2p (*1360,*680,*1100)
GO TO 1230
1030 CALL ifs3p (*1360,*680,*1100)
GO TO 1230
1040 CALL ifs4p (*1360,*680,*1100)
GO TO 1230
1050 CALL ifs5p (*1360,*680,*1100)
GO TO 1230
1060 CALL page2 (2)
WRITE  (nout,1070) sfm,k
1070 FORMAT (a25,' 314, INVALID CALL FROM IFP.  K =',i10)
abort =.true.
GO TO 1850

1100 IF (.NOT.badfor) GO TO 1160
IF (it1k /= t1(1,k) .OR. it2k /= t1(2,k)) GO TO 1110
icount = icount + 1
IF (icount-7 < 0) THEN
  GO TO  1120
ELSE IF (icount-7 == 0) THEN
  GO TO  1140
ELSE
  GO TO  1170
END IF
1110 it1k = t1(1,k)
it2k = t1(2,k)
1120 CALL page2 (2)
IF (id == 0) id = m(1)
WRITE  (nout,1130) ufm,t1(1,k),t1(2,k),id,knt
1130 FORMAT (a23,' 315, FORMAT ERROR ON BULK DATA CARD ',2A4,i8,17X,  &
    'SORTED CARD COUNT =',i7)
GO TO 1170
1140 CALL page2 (3)
WRITE  (nout,1150)
1150 FORMAT (31X,'.', /29X,'MORE', /31X,'.')
GO TO 1170
1160 IF (.NOT.baddat) icount = 0
1170 IF (.NOT.baddat) GO TO 1220
IF (jt1k /= t1(1,k) .OR. jt2k /= t1(2,k)) GO TO 1180
jcount = jcount + 1
IF (jcount-7 < 0) THEN
  GO TO  1190
ELSE IF (jcount-7 == 0) THEN
  GO TO  1210
ELSE
  GO TO  1230
END IF
1180 jt1k = t1(1,k)
jt2k = t1(2,k)
1190 CALL page2 (2)
IF (id == 0) id = m(1)
WRITE  (nout,1200) ufm,t1(1,k),t1(2,k),id,knt
1200 FORMAT (a23,' 316, ILLEGAL DATA ON BULK DATA CARD ',2A4,i8,17X,  &
    'SORTED CARD COUNT =',i7)
GO TO 1230
1210 CALL page2 (3)
WRITE (nout,1150)
GO TO 1230
1220 IF (.NOT.badfor) jcount = 0
1230 IF (.NOT.badfor .AND. .NOT.baddat) GO TO 1300
n = 0
abort =.true.
GO TO 1340

!     WRITE OUT CARD DATA ON APPROPRIATE IFP OUTPUT FILE

1300 IF (n == 0) GO TO 1340
t4(2,k) = t4(2,k) + n
DO  l = 1,40
  IF (k == kkl(l)) GO TO 1320
END DO
GO TO 1330
1320 CALL WRITE (ifle(curfil),i,n,0)
GO TO 1340
1330 CONTINUE
IF (indx /= 6 .AND. .NOT.abort .OR. indx == 15)  &
    CALL WRITE (ifle(curfil),i,n,0)
1340 IF (kn == 0) GO TO 680
kn = 0
GO TO 430
1350 badfor =.true.
1360 IF (.NOT.badfor) GO TO 1370
CALL page2 (2)
WRITE (nout,1130) ufm,t1(1,k),t1(2,k),m(1),knt
abort = .true.
1370 IF (.NOT.baddat) GO TO 1380
CALL page2 (2)
WRITE (nout,1200) ufm,t1(1,k),t1(2,k),m(1),knt
abort =.true.
GO TO 680
1380 IF (iapp  == 1) GO TO 1850
IF (isubs /= 0) GO TO 1850
1390 WRITE  (nout,1400) sfm
1400 FORMAT (a25,' 319, IFP READING EOF ON NPTP.')
CALL page2 (2)
abort =.true.
GO TO 1850
1410 kerror = 1410
GO TO 1440
1420 kerror = 1420
GO TO 1440
1430 kerror = 1430
1440 CALL page2 (6)
WRITE  (nout,1450) sfm,kerror,(jr(l),l=1,20),knt
1450 FORMAT (a25,' 320, IFP ERROR',i5, /5X,'LAST CARD PROCESSED IS -',  &
    20A4,1H-, /5X,'SORTED CARD COUNT =',i7)
abort =.true.
GO TO 1850
1460 IF (curfil /= 0) CALL CLOSE (ifle(curfil),2)
DO  l = 1,nfls
  IF (l == 6 .OR. l == 12 .OR. STATUS(l) == 1) CYCLE
  kfil = ifle(l)
  CALL OPEN  (*130,kfil,ibuff,3)
  CALL WRITE (ifle(l),iend,3,1)
  ii(l) = ii(l) + 1
  CALL CLOSE (ifle(l),1)
END DO

!     CHECK TO SEE IF ALL MULTI-ENTRY CARD DATA (CROD, CTUBE, ETC.)
!     ARE SORTED ON THEIR ELEMENT/PROPERTY IDS

DO  l = 1,40
  IF (endara(l) < 0) GO TO 1490
END DO

!     EITHER NO MULTI-ENTRY CARD DATA EXIST OR, IF THEY DO,
!     THEY ARE ALL SORTED ON THEIR ELEMENT/PROPERTY IDS

GO TO 1700

!     NOT ALL MULTI-ENTRY CARD DATA ARE SORTED ON THEIR
!     ELEMENT/PROPERTY IDS.

!     CLOSE SCRATCH FILE (FILE 6) AT CURRENT POSITION WITHOUT REWIND
!     AND WITHOUT END-OF-FILE.

1490 CALL CLOSE (ifle(6),2)

!     READ DATA FROM GEOM2/EPT FILE, SORT ALL MULTI-ENTRY CARD DATA ON
!     THEIR ELEMENT/PROPERTY IDS AND WRITE THE RESULTING DATA ON
!     SCRATCH FILE (FILE 16)

!     NOTE.  GEOM2 IS IFLE(8) AND EPT IS IFLE(2)

DO  nnn = 1,2
  IF (nnn == 2) GO TO 1500
  ifile   = ifle(8)
  inam(1) = fnm(1,8)
  inam(2) = fnm(2,8)
  itype(1)= itype1(1)
  itype(2)= itype1(2)
  jmin    = 1
  jmax    = 20
  GO TO 1510
  1500 ifile   = ifle(2)
  inam(1) = fnm(1,2)
  inam(2) = fnm(2,2)
  itype(1)= itype2(1)
  itype(2)= itype2(2)
  jmin    = 21
  jmax    = 40
  1510 DO  l = jmin,jmax
    IF (endara(l) < 0) GO TO 1530
  END DO
  CYCLE
  1530 ileft = nopen - nparam - 2
  istrt = 2*n1  + nparam + 2
  CALL gopen (ifile,ibuff,0)
  kfil = ifle(16)
  CALL OPEN  (*130,ifle(16),ibuff(n1+1),1)
  CALL WRITE (ifle(16),inam,2,1)
  INDEX = jmin
  1540 CALL READ (*1680,*1670,ifile,ibuff(istrt),3,0,iflag)
  CALL WRITE (ifle(16),ibuff(istrt),3,0)
  IF (INDEX > jmax) GO TO 1560
  DO  l = jmin,jmax
    IF (ibuff(istrt) == mentry(l) .AND. endara(l) < 0) GO TO 1580
  END DO
  1560 CALL READ  (*1660,*1570,ifile,ibuff(istrt),ileft,0,iflag)
  CALL WRITE (ifle(16),ibuff(istrt),ileft,0)
  GO TO 1560
  1570 CALL WRITE (ifle(16),ibuff(istrt),iflag,1)
  GO TO 1540
  1580 INDEX = INDEX + 1
  CALL page2 (3)
  WRITE  (nout,1590) uim,nentry(2*l-1),nentry(2*l),itype
  1590 FORMAT (a29,' 334, ',2A4,' MULTI-ENTRY CARD DATA ARE NOT SORTED ',  &
      'ON THEIR ',2A4,' IDS.', /5X, 'SUBROUTINE IFP WILL SORT THE DATA.')
  ifail = 0
  1600 CALL READ (*1660,*1610,ifile,ibuff(istrt),ileft,0,iflag)
  ifail = ifail + 1
  GO TO 1600
  1610 IF (ifail == 0) GO TO 1630
  nwds = (ifail-1)*ileft + iflag
  CALL page2 (4)
  WRITE  (nout,1620) ufm,nentry(2*l-1),nentry(2*l),nwds
  1620 FORMAT (a23,' 333, UNABLE TO SORT ',2A4,' MULTI-ENTRY CARD DATA ',  &
      'IN SUBROUTINE IFP DUE TO INSUFFICIENT CORE.', /5X,  &
      'ADDITIONAL CORE REQUIRED =',i10,7H  words)
  CALL mesage (-61,0,0)
  1630 nwds = 4
  IF (l == 10 .OR. l == 33) nwds = 3
  IF (l == 21 .OR. l == 23) nwds = 2
  CALL sort  (0,0,nwds,1,ibuff(istrt),iflag)
  CALL WRITE (ifle(16),ibuff(istrt),iflag,1)
  
!     CHECK SORTED MULTI-ENTRY CARD DATA FOR NON-UNIQUE
!     ELEMENT/PROPERTY IDS
  
  irept = -10000000
  nidsm1= iflag/nwds - 1
  DO  kk = 1,nidsm1
    eid   = ibuff(istrt+kk*nwds)
    eidm1 = ibuff(istrt+kk*nwds-nwds)
    IF (eid /= eidm1) CYCLE
    IF (eid == irept) CYCLE
    irept = eid
    abort = .true.
    CALL page2 (2)
    WRITE  (nout,1640) ufm,itype,eid,nentry(2*l-1),nentry(2*l)
    1640 FORMAT (a23,' 335, NON-UNIQUE ',2A4,' ID',i9,' ENCOUNTERED IN ',  &
        2A4,' MULTI-ENTRY CARD DATA.')
  END DO
  GO TO 1540
  1660 CALL mesage (-2,ifile,NAME)
  1670 CALL mesage (-3,ifile,NAME)
  1680 CALL CLOSE  (ifile,1)
  CALL CLOSE  (ifle(16),1)
  
!     COPY DATA BACK FROM SCRATCH FILE (FILE 16) TO GEOM2/EPT FILE
  
  kfil = ifle(16)
  CALL OPEN (*130,ifle(16),ibuff,0)
  kfil = ifile
  CALL OPEN (*130,ifile,ibuff(n1+1),1)
  CALL cpyfil (ifle(16),ifile,ibuff(istrt),ileft,iflag)
  CALL CLOSE (ifle(16),1)
  CALL CLOSE (ifile,1)
END DO

!     RE-OPEN SCRATCH FILE (FILE 6) TO WRITE WITHOUT REWIND

kfil = ifle(6)
CALL OPEN (*130,ifle(6),ibuff(n1+1),3)

!     WRITE TRAILERS

1700 DO  j = 1,nfls
  IF (j == 6 .OR. j == 12) CYCLE
  DO  l = 2,7
    itrl(l) = 0
  END DO
  itrl(1) = ifle(j)
  IF (ii(j) <= 2 .OR. abort) GO TO 1730
  DO  l = 1,ncdsmx
    IF (t3(1,l) /= j .OR. t4(2,l) <= 0) CYCLE
    kt721 = andf(t7(2,l),511)
    k1 = (kt721-1)/16 + 2
    k2 = kt721 - (k1-2)*16 + 16
    itrl(k1) = orf(itrl(k1),two(k2))
  END DO
  1730 CALL wrttrl (itrl)
END DO

!     WRITE PARAM CARDS ON NPTP

kfil = ifle(16)
CALL ifppar
IF (nparam <= 0 .OR. abort) GO TO 1850
CALL OPEN (*130,kfil,ibuff,1)
itrl(1) = kfil
itrl(2) = nparam
CALL wrttrl (itrl(1))
CALL WRITE (kfil,fnm(1,6),2,1)
CALL WRITE (kfil,ibuff(2*n1+1),nparam,1)
ipm   = 1
ipn   = 2*n1 + ipm
GO TO 1840
1800 ipn   = 2*n1 + ipm
nm(1) = ibuff(ipn  )
nm(2) = ibuff(ipn+1)
jpm   = 1
1810 jpn   = 2*n1 + jpm
IF (nm(1) /= ibuff(jpn) .OR. nm(2) /= ibuff(jpn+1)) GO TO 1830
CALL page2 (2)
WRITE  (nout,1820) ufm,nm(1),nm(2)
1820 FORMAT (a23,' 321, NON-UNIQUE PARAM NAME - ',2A4,1H-)
abort =.true.
1830 jpm   = jpm + 4
IF (ibuff(jpn+2) > 2) jpm = jpm + 1
IF (ibuff(jpn+2) > 5) jpm = jpm + 2
IF (jpm < ipm) GO TO 1810
1840 ipm = ipm + 4
IF (ibuff(ipn+2) > 2) ipm = ipm + 1
IF (ibuff(ipn+2) > 5) ipm = ipm + 2
IF (ipm < nparam) GO TO 1800
CALL eof (kfil)
CALL CLOSE (kfil,1)
1850 CALL CLOSE (ifle(6),1)

!     CHECK FOR PROPERTY ID UNIQUENESS IN EPT FILE AND PROPERTY ID
!     SPECIFIED IN GEOM2 ELEMENTS

CALL sswtch (34,jj1)
IF (jj1 == 1) GO TO 1900
kfil = ifle(2)
itrl(1) = kfil
CALL rdtrl (itrl)
j = itrl(2) + itrl(3) + itrl(4) + itrl(5) + itrl(6) + itrl(7)
jj1 = 1
IF (itrl(1) < 0 .OR. j == 0) GO TO 1860
jj1 = 0
CALL OPEN (*130,kfil,ibuff,0)
1860 kfil = ifle(8)
itrl(1) = kfil
CALL rdtrl (itrl)
j = itrl(2) + itrl(3) + itrl(4) + itrl(5) + itrl(6) + itrl(7)
IF (itrl(1) < 0 .OR. j == 0) GO TO 1880
CALL OPEN (*130,kfil,ibuff(n1+1),0)
jj = n1*2 + 1
CALL pidck (ifle(2),kfil,jj1,ibuff(jj))
CALL CLOSE (kfil,1)
IF (jj1 < 0) THEN
  GO TO  1880
ELSE IF (jj1 == 0) THEN
  GO TO  1870
ELSE
  GO TO  1890
END IF
1870 IF (ibuff(jj) == 0) GO TO 1880
jj1 = jj + ibuff(jj) + 1

!     CHECK FOR MATERIAL ID UNIQUENESS IN MPT FILE
!     AND MATERIAL ID SPECIFIED IN PROPERTY CARDS

kfil = ifle(3)
itrl(1) = kfil
CALL rdtrl (itrl)
j = itrl(2) + itrl(3) + itrl(4) + itrl(5) + itrl(6) + itrl(7)
ibuff(jj1) = 1
IF (itrl(1) < 0 .OR. j == 0) ibuff(jj1) = 0
IF (ibuff(jj1) == 1) CALL OPEN (*130,kfil,ibuff(n1+1),0)
CALL matck (kfil,ifle(2),ibuff(jj),ibuff(jj1))
IF (ibuff(jj1) /= 0) CALL CLOSE (kfil,1)
1880 CALL CLOSE (ifle(2),1)

!     CHECK COORDINATE ID'S AND THEIR REFERENCES FROM
!     OTHER BULK DATA CARDS

1890 jj = nopen + n1 - 2
!                + N1 - 2 = 2*N1 - (N1+2)
CALL cidck (ibuff(n1+2),ibuff,jj)
1900 CONTINUE

!     CHECK FOR ERRORS IN AXISYMMETRIC DATA

IF (iax) axiccc = 1
axifcc = iaxf
IF (axiccc <= 0 .OR. axifcc <= 0) GO TO 1920
axiccc = 0
axifcc = 0
abort  = .true.
CALL page2 (2)
WRITE (nout,1910) ufm
1910 FORMAT (a23,' 337, BOTH AXIC AND AXIF CARDS USED IN BULK DATA.')
GO TO 1980
1920 IF (axiccc <= 0) GO TO 1950
IF (iaxic  > 0) GO TO 1980
axiccc = 0

!     SUPPRESS ABORT IF IT IS A UMFEDIT RUN

1930 IF (iumfed /= 0) GO TO 1980
abort = .true.
CALL page2 (2)
WRITE (nout,1940) ufm
1940 FORMAT (a23,' 338, AXISYMMETRIC CARD REQUIRED IN CASE CONTROL')
GO TO 1980
1950 IF (axifcc <= 0) GO TO 1960
IF (iaxif > 0 .OR. axifcc == 2) GO TO 1980
axifcc = 0
GO TO 1930
1960 IF (iaxic <= 0 .AND. iaxif <= 0) GO TO 1980
axiccc = 0
axifcc = 0

!     SUPPRESS ABORT IF IT IS A UMFEDIT RUN

IF (iumfed /= 0) GO TO 1980
abort = .true.
CALL page2 (2)
WRITE  (nout,1970) ufm
1970 FORMAT (a23,' 339, ILLEGAL USE OF AXISYMMETRIC CARD IN CASE ',  &
    'CONTROL DECK.')

1980 IF (iapp >= 0) GO TO 1990

!     CHECK CERTAIN RESTART FLAGS BASED ON BULK DATA

mn = lbd + 1

!     TURN ON TEMPMX$ IF MATERIALS USE TEMPS

IF (t4(2,91)+t4(2,102)+t4(2,189) == 0) GO TO 1990
IF (andf(ib(1),two(28)) == 0 .AND. andf(ib(5),two(32)) == 0 .AND.  &
    andf(ib(4),two( 6)) == 0 .AND. andf(ib(3),two(32)) == 0 .AND.  &
    andf(ib(4),two( 2)) == 0 .AND. andf(ib(4),two( 3)) == 0 .AND.  &
    andf(ib(4),two( 4)) == 0 ) GO TO 1990
ib(mn) = orf(ib(mn),two(19))
1990 CALL conmsg (ifpna2,2,0)

CALL sswtch (27,l27)
IF (l27 == 0) GO TO 2060
CALL page1
line = line + 8
WRITE  (nout,2000)
2000 FORMAT ('0DIAG 27 DUMP OF IFP TABLES AFTER IFP PROCESSING',     /,  &
    1H0,6X,6HIFX1BD,9X,6HIFX2BD,7X,6HIFX3BD,2X,6HIFX4BD,3X,  &
    6HIFX5BD,6X,6HIFX6BD                            ,/,  &
    1H ,5X,8(1H-),2X,17(1H-),2X,6(1H-),2X,6(1H-),2X,8(1H-),  &
    2X,12(1H-)                                         ,/,  &
    1H ,1X,3H(a),3X,3H(b),5X,3H(c),3X,3H(d),5X,3H(e),2X,3H(n),  &
    5X,3H(f),   3H(g),3X,3H(h),1X,3H(i),3X,3H(j),5X,3H(k),  &
    4X,3H(l),3X,3H(m),4X,3H(o),3X,4HFLAG               ,/ 1H0 )
2010 FORMAT (1H ,i4,1X,2A4,i4,1X,1H(,2A4,1H),i3,i5,i4,  &
    i4,i4,i6,i3,i7,i8,16X, i1,i1,a1,i1,4X,i2)
DO  j = 1,ncdsmx
  id = t3(1,j)
  IF (id <= 0) GO TO 2020
  lf = fnm(1,id)
  lm = fnm(2,id)
  GO TO 2030
  2020 CONTINUE
  lf = BLANK
  lm = BLANK
  2030 CONTINUE
  n  = j
  k  = n/90 + MIN0(1,MOD(n,90))
  n  = n - 90*(k-1)
  kx = n/30 + MIN0(1,MOD(n,30))
  n  = n - 30*(kx-1)
  ky = n/6  + MIN0(1,MOD(n, 6))
  l  = n - 6*(ky-1)
  iflag = 0
  IF (eject(1) == 0) GO TO 2040
  WRITE (nout,2000)
  line = line + 8
  2040 CONTINUE
  line = line + 1
  WRITE (nout,2010) j,t1(1,j),t1(2,j), t3(1,j),lf,lm,t3(2,j),  &
      t4(1,j),t4(2,j), t5(1,j),t5(2,j),  &
      t6(1,j),t6(2,j), t7(1,j),t7(2,j),  &
      k,kx,ooo(ky),l,iflag
END DO

2060 RETURN
END SUBROUTINE ifp
