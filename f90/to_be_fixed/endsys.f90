SUBROUTINE endsys (jobseg,jobend)
     
!     ENDSYS SAVES VARIOUS EXEC TABLES ON A SCRATCH FILE
 
!     LAST REVISED  5/91 BY G.CHAN/UNISYS  FOR SUPERLINK OPERATION
!          IF SPERLK = 0, WE ARE IN NASTRAN MULTI-LINK COMPUTATION
!          IF SPERLK = NON-ZERO, WE ARE IN NASTRAN SUPERLINK
!          SPERLK IS THE 95TH WORD OF /SYSTEM/
 
 
 INTEGER, INTENT(IN)                      :: jobseg
 INTEGER, INTENT(IN OUT)                  :: jobend
 EXTERNAL        lshift,rshift,andf,orf,link
 LOGICAL :: bitpas
 INTEGER :: andf,fist,SAVE,scrn1,scrn2,thcrmk,pool,sperlk,  &
     nopref(2),rshift,buf,msgbuf(8),bcdnum(10),units,  &
     tens,orf,unitab(75),fcb(75),databf,msg(2),NAME(2),  &
     FILE,filex,lnknum(15),comm,xf1at,prefac
 CHARACTER (LEN=7) :: fortxx
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /BLANK / iblkcm(58),prefac(2)
 COMMON /xpfist/ npfist
 COMMON /xfist / fist(2)
 COMMON /msgx  / itab1(1)
 COMMON /stime / itab2(2)
 COMMON /stapid/ itab3(1)
 COMMON /xdpl  / itab4(3)
 COMMON /xxfiat/ itab5(1)
 COMMON /xfiat / itab6(4)
 COMMON /xvps  / itab7(2)
 COMMON /xceitb/ itab8(2)
 COMMON /ginox / itab9(170)
 COMMON /system/ itab10(22),lsystm,icfiat,jtab10(11),nprus,  &
     ktab10(35),bitpas,ltab10(18),lpch,ldict,mtab10(2), sperlk,ntab10(5)
 COMMON /output/ itab11(1)
 COMMON /ntime / itab13(1)
 COMMON /xlink / itab14(1)
 COMMON /sofcom/ itab15(1)
 COMMON /bitpos/ bt(32,2)
 COMMON /oscent/ inoscr(2)
!ZZ   COMMON /ZZENDS/ DATABF(1)
 COMMON /zzzzzz/ databf(1)
 COMMON /sem   / mask,thcrmk,imask,links(15)
 COMMON /l15 l8/ l15,l8,l13
 COMMON /xsfa1 / dummy(1902),comm(20),xf1at(1100)
!                           1902 = 401+1501
 
 EQUIVALENCE    (itab10( 1),isybuf),  (itab10(2),nout     ),  &
     (itab9 ( 2),filex ),  (itab9(12),unitab(1)), (itab9(170),fcb(1))
 DATA msgbuf(1)/ 4HLINK /
 DATA msgbuf(3)/ 4H     /
 DATA msgbuf(5)/ 4H---- /
 DATA msgbuf(6)/ 4H---- /
 DATA msgbuf(7)/ 4H---- /
 DATA msgbuf(8)/ 4H---- /
 DATA scrn1    , scrn2  /4HSCRA,4HTCH0/, SAVE/4HSAVE/, pool     / 4HPOOL /,  &
     nopref   / 4HNOT , 4HPREF/
DATA msg      / 4HBEGN, 4HEND /
DATA bcdnum   / 1H0, 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9 /
DATA lnknum   / 4H 1  , 4H 2  , 4H 3  , 4H 4  , 4H 5  ,  &
    4H 6  , 4H 7  , 4H 8  , 4H 9  , 4H10  ,  &
    4H11  , 4H12  , 4H13  , 4H14  , 4H15  /
DATA NAME     / 4HENDS,4HYS   /


!     PUNCH RESTART DICTIONAY
!     LDICT MAY NOT BE A SYSTEM PUNCH FILE, PUNCH THE CARDS OUT FIRST
!     BEFORE THE RESTART DICTIONARY CARDS GET LOST

IF (mach >= 5 .OR. ldict == lpch) GO TO 8
ENDFILE ldict
REWIND  ldict
5 READ   (ldict,6,ERR=7,END=7) (databf(j),j=1,20)
6 FORMAT (20A4)
WRITE  (lpch,6) (databf(j),j=1,20)
GO TO  5
7 REWIND ldict

8 msgbuf(2) = 0
j = 0
DO  i = 1,15
  IF (jobend == links(i)) msgbuf(2) = lnknum(i)
  IF (jobseg == links(i)) j = i
END DO
IF (msgbuf(2) /= 0) GO TO 15
WRITE  (nout,12) sfm,jobend
12 FORMAT (a25,', ILLEGAL LINK NUMBER ',a4,' ENCOUNTERED BY ENDSYS')
CALL mesage (-61,0,0)
15 msgbuf(4) = msg(2)

IF (sperlk == 0) GO TO 30

!     SIMPLIFIED OPERATION IF SUPERLINK (USED IN UNIX VERSION)

sperlk = j
itab10(22) = jobseg
!     PREFAC(1)  = NOPREF(1)
!     PREFAC(2)  = NOPREF(2)
CALL conmsg (msgbuf   ,4,0)
CALL conmsg (msgbuf(5),4,0)
DO  j = 2,11
  itab9(j) = 0
END DO
DO  j = 87,161
  IF (itab9(j) == 0) CYCLE
  i = j - 86
  WRITE  (nout,23) sfm,i,jobend
23 FORMAT (a25,', LOGICAL UNIT',i5,' WAS NOT CLOSED AT END OF ',a4)
!     ITAB9(J) = 0
CALL mesage (-37,0,0)
END DO
GO TO 400

!     SEARCH FIAT FOR A SAVE FILE -- FILE MUST SATISFY THE FOLLOWING
!     (1) FILE MUST BE SCRATCHX OR TRAILERS=0 OR EXPIRED*
!     (2) IF (1) IS TRUE, NO UNEXPIRED SECONDARY ALLOCATIONS WITH
!     NON-ZERO TRAILERS MAY EXIST. (ALSO FILE MUST NOT BE PURGED)
!     AN EXPIRED FILE HAS AN LTU LESS THAN THE CURRENT OSCAR POSITION.

30 FILE = SAVE
lmt  = itab6(3)*icfiat + 3
next = lshift(inoscr(2),16)
ifound  = 0
fist(2) = npfist + 1
fist(2*npfist+3) = SAVE

k = andf(thcrmk,scrn2)
loop50:  DO  i = 4,lmt,icfiat
  IF (itab6(i+1) == scrn1 .AND. andf(thcrmk,itab6(i+2)) == k) GO TO 35
  IF (itab6(i+3) /= 0 .OR. itab6(i+4) /= 0 .OR. itab6(i+5) /= 0) GO TO 32
  IF (icfiat == 11 .AND.  (itab6(i+8) /= 0 .OR. itab6(i+9) /= 0 .OR.  &
      itab6(i+10) /= 0)) GO TO 32
  GO TO 35
  32 ltu = andf(itab6(i),1073676288)
!                         1073676288 = 2**30 - 2**16 = 3FFF0000 HEX
!                                    = 0 SIGN BIT + LEFT 14 BITS OF 1's
  IF (ltu >= next .OR. ltu == 0) CYCLE loop50
  35 iucb = andf(itab6(i),32767)
!                          32767 = 2**15 - 1 = RIGHT 15 BITS OF 1's
  IF (iucb == 32767) CYCLE loop50
  DO  j = 4,lmt,icfiat
    IF (andf(itab6(j),32767) /= iucb) CYCLE
    IF (i == j) CYCLE
    ltu = andf(itab6(j),1073676288)
    IF (ltu < next .AND. ltu /= 0) CYCLE
    IF (itab6(j+3) /= 0 .OR. itab6(j+4) /= 0 .OR. itab6(j+5) /= 0)  &
        CYCLE loop50
    IF (icfiat == 11 .AND. (itab6(j+8) /= 0 .OR. itab6(j+9) /= 0 .OR.  &
        itab6(j+10) /= 0)) CYCLE loop50
  END DO
  IF (ifound == 0) ifound = i
  
!     FLUSH FILE IN CASE DATA EXISTS ON FILE
!     THIS WILL FREE UP SECONDARIES ON 360 AND DISK ON CDC AND UNIVAC
  
  IF (itab6(i+3) /= 0 .OR. itab6(i+4) /= 0 .OR. itab6(i+5) /= 0) GO TO 45
  IF (icfiat == 11 .AND. (itab6(i+8) /= 0 .OR. itab6(i+9) /= 0 .OR.  &
      itab6(i+10) /= 0)) GO TO 45
  CYCLE loop50
  45 fist(2*npfist+4) = i - 1
  CALL OPEN (*360,SAVE,databf,1)
  CALL CLOSE (SAVE,1)
END DO loop50

IF (ifound == 0) CALL mesage (-39,0,0)
i = -2
IF (itab11(1)+itab11(-i) == i) icfiat = icfiat + i

!     GOOD NEWS - WE FOUND A FILE FOR SAVE PURPOSES.
!     SAVE POINTER TO FILE IN BLANK COMMON.

i = ifound
iblkcm(1) = itab6(i)

!     SAVE UNIT = 2 FOR ALL MACHINES, IBM INCLUDED
!     (IBM USED 51 BEFORE)

iunitu = 2

!     FCB ARREY OF 75 WORDS IS NOT USED BY VAX AND UNIX

REWIND iunitu
IF (mach < 5) WRITE (iunitu) itab6(i),isybuf,fcb
IF (mach >= 5) WRITE (iunitu) itab6(i),isybuf
REWIND iunitu
fist(2*npfist+4) = i - 1

!     SET PREFAC FLAG SO LINK 1 IS RE-ENTRANT

!     PREFAC(1) = NOPREF(1)
!     PREFAC(2) = NOPREF(2)

!     SAVE THE NEXT LINK NO. IN THE 22ND WORD OF /SYSTEM/

itab10(22) = jobseg

!     WRITE EXEC TABLES ON THE FILE JUST FOUND.

CALL OPEN  (*360,SAVE,databf,1)
ltab10(7)  = 0
CALL WRITE (SAVE,itab10,lsystm,1)
CALL WRITE (SAVE,itab1,itab1(1)*4+2,1)
CALL WRITE (SAVE,itab2,1,1)
CALL WRITE (SAVE,itab3,6,1)
CALL WRITE (SAVE,itab4,itab4(3)*3+3,1)
CALL WRITE (SAVE,itab5,npfist,1)
CALL WRITE (SAVE,itab6,itab6(3)*icfiat+3,1)
CALL WRITE (SAVE,itab7,itab7(2),1)
CALL WRITE (SAVE,itab8,itab8(2),1)
CALL WRITE (SAVE,itab9(12),75,1)
CALL WRITE (SAVE,itab11,224,1)
CALL WRITE (SAVE,itab13,itab13(1)+1,1)
CALL WRITE (SAVE,itab14,itab14(1)+2,1)
CALL WRITE (SAVE,itab15,27,1)
CALL WRITE (SAVE,bt,64,1)
CALL CLOSE (SAVE,1)

!     FLUSH ANY QUEUED SYSTEM OUTPUT.
!     LOAD NEXT LINK NO. INTO UNIT 97, AND TERMINATE PRESENT LINK.

kk = itab10(2)
WRITE  (kk,55)
55 FORMAT (//)
CALL conmsg (msgbuf   ,4,0)
CALL conmsg (msgbuf(5),4,0)
IF (mach == 4) GO TO 67
IF (itab10(7) < 0) ENDFILE 52

!     IF IBM NEW LOGIC OF LINK SWITCHING VIA FILE 97 IS NOT AVAILBLE,
!     WE STILL NEED THE NEXT 3 LINES FOR DEAR OLD IBM

IF (mach /= 2) GO TO 60
!     CALL SEARCH (JOBSEG,SYSLB2,NOTUSE)
CALL search (jobseg)
GO TO 400

60 i = khrfn3(msgbuf(3),jobseg,2,1)
IF (mach == 9 .OR. mach == 12) GO TO 61
OPEN (UNIT=97,ACCESS='SEQUENTIAL',STATUS='NEW',ERR=64)
GO TO 62
61 CALL flunam (97,fortxx)
OPEN (UNIT=97,ACCESS='SEQUENTIAL',STATUS='NEW',ERR=64,FILE=fortxx)
62 WRITE  (97,63) i
63 FORMAT ('NAST',a2)
CLOSE (UNIT=97)
CALL EXIT
!SUN  CALL EXIT (0)
64 WRITE  (nout,65)
65 FORMAT ('0*** SYSTEM ERROR, CAN NOT OPEN FORTRAN UNIT 97 FOR ',  &
    'LINK SWITCH')
CALL mesage (-37,0,NAME)

!     DETERMINE LINK NUMBER FOR 6600

67 i = andf(4095,rshift(jobseg,36))
i1 = i/64
i2 = i - i1*64
i  = 10*i1 + i2 - 297
i76= 76
CALL link (i,itab10(i76),0)
GO TO 350


ENTRY bgnsys
!     ============

nprus      = 0
bitpas     = .true.
msgbuf(4)  = msg(1)
!     PREFAC(1)  = 0
!     PREFAC(2)  = 0
IF (sperlk == 0) GO TO 70

!     SIMPLEFIED OPERATION IF SUPERLINK (USED IN UNIX VERSION)

IF (sperlk < 1 .OR. sperlk > 15) GO TO 225
itab10(22) = links(sperlk)
msgbuf(2)  = lnknum(sperlk)
jobsxx     = itab10(22)
GO TO 228

!     BGNSYS RESTORES THE EXEC TABLES SAVED BY ENDSYS
!     THEN REPOSITIONS THE OSCAR TO THE ENTRY FOR THE MODULE
!     IN THE CURRENT LINK.

70 iunitu = 2
IF (mach < 5) READ (iunitu) itab6(4),isybuf,fcb
IF (mach >= 5) READ (iunitu) itab6(4),isybuf
fist(2) = npfist + 1
fist(2*npfist+3) = SAVE
fist(2*npfist+4) = 3
j = 5000
CALL OPEN (*360,SAVE,databf(j),0)
CALL READ (*340,*80,SAVE,itab10,900,1,flg)
80 CALL READ (*340,*90,SAVE,itab1,900,1,flg)
GO TO 350
90 CALL READ (*340,*100,SAVE,itab2,900,1,flg)
GO TO 350
100 CALL READ (*340,*110,SAVE,itab3,900,1,flg)
GO TO 350
110 CALL READ (*340,*120,SAVE,itab4,900,1,flg)
GO TO 350
120 CALL READ (*340,*130,SAVE,itab5,900,1,flg)
GO TO 350
130 CALL READ (*340,*140,SAVE,itab6,900,1,flg)
GO TO 350
140 CALL READ (*340,*150,SAVE,itab7,900,1,flg)
GO TO 350
150 CALL READ (*340,*160,SAVE,itab8,900,1,flg)
GO TO 350
160 CALL READ (*340,*170,SAVE,itab9(12),900,1,flg)
GO TO 350
170 CALL READ (*340,*190,SAVE,itab11,900,1,flg)
GO TO 350
190 CALL READ (*340,*210,SAVE,itab13,900,1,flg)
GO TO 350
210 CALL READ (*340,*220,SAVE,itab14,900,1,flg)
GO TO 350
220 CALL READ (*340,*221,SAVE,itab15,900,1,flg)
GO TO 350
221 CALL READ (*340,*222,SAVE,bt,900,1,flg)
GO TO 350
222 CALL CLOSE (SAVE,1)

!     RETRIEVE THE CURRENT LINK NO. FROM THE 22ND WORD OF /SYSTEM/

jobsxx = itab10(22)
DO  i = 1,15
  IF (jobsxx /= links(i)) CYCLE
  msgbuf(2) = lnknum(i)
  GO TO 228
END DO
225 WRITE  (nout,226) sfm,jobsxx,sperlk
226 FORMAT (a25,', ILLEGAL LINK NUMBER ',a4,' ENCOUNTERED BY BGNSYS.',  &
    4X,'SPERLK=',i14)
CALL mesage (-61,0,0)

228 CALL pressw (jobsxx,i)
CALL conmsg (msgbuf,4,0)
CALL sswtch (15,l15)
CALL sswtch ( 8,l 8)
CALL sswtch (13,l13)
IF (mach /= 3) GO TO 320

IF (itab10(7) >= 0) GO TO 238
232 READ (52,234,END=236) i
234 FORMAT (a1)
GO TO 232
236 BACKSPACE 52
238 CONTINUE

!     REPOSITION DRUM FILES OFF LOAD POINT (1108 ONLY)

CALL defcor
CALL contin

!     TAPE-FLAG IS THE 45TH WORD OF /SYSTEM/
!     IF THE 7TH BIT (COUNTING FROM RIGHT TO LEFT) OF TAPE-FLAG IS NOT
!     ON (=1), AND PLT2 HAS NOT BEEN EXTERNALLY ASSIGNED AS A MAGNETIC
!     TAPE, SET PLT2 IS TO DISK. SIMILARILY,
!     IF THE 6TH BIT IS NOT SET, AND PLT1 IS NOT TAPE ASSIGNED, SET PLT1
!     TO DISK

i45   = 45
istat = andf(itab10(i45),64)
jstat = andf(itab10(i45),32)

DO  i = 1,75
  
!     CALL FACIL TO DETERMINE IF UNIT IS TAPE
  
  tens  = i/10
  units = i - 10*tens
  nbcd  = bcdnum(units+1)
  IF (tens == 0) GO TO 295
  maskk = 255
  maskk = lshift(maskk,27)
  nbcd  = orf(andf(bcdnum(tens+1),maskk),rshift(nbcd,9))
  295 CALL facil (nbcd,j)
  
!     DECODE UNITAB ENTRY
  
  nblock = andf(rshift(unitab(i),12),262143)
  nlr = andf(unitab(i),4095)
  IF (j == 7 .OR. j == 9) GO TO 298
  
!     POSITION DRUM UNIT NOW OFF LOAD POINT
  
  IF (nblock+nlr == 1) CYCLE
  CALL ntran (i,10,22)
  nosect = nblock*itab9(164)
  IF (nlr == 0) nosect = nosect - itab9(164)
  IF (i == 13 .AND. istat /= 0) nosect = unitab(13)
  IF (i == 12 .AND. jstat /= 0) nosect = unitab(12)
  CALL ntran (i,6,nosect)
  
!     RESET FCB ENTRY
!     COMMENTS FROM G.CHAN/UNISYS   11/90
!     FCB ARRAY OF 75 WORDS IS USED ONLY BY UNIVAC AND IBM. IT BEGINS
!     AT THE 170TH WORD OF /GINOX/
  
  298 IF (nlr /= 0) nblock = nblock + 1
  fcb(i) = nblock
END DO

320 IF (sperlk /= 0) GO TO 330

!     DEFINE OPEN CORE FOR VAX AND UNIX

IF (mach >= 5) CALL defcor

!     REPOSITION POOL TO OSCAR ENTRY TO BE EXECUTED.

330 buf  = korsz(databf) - itab10(1)
FILE = pool
CALL OPEN (*360,pool,databf(buf),2)
CALL bckrec (pool)
IF (sperlk == 0) GO TO 400
DO  j = 1,60
  iblkcm(j)= 0
END DO
!     DO 334 J = 1,1902
! 334 DUMMY(J) = 0
!     COMM( 1) = 0
!     COMM( 3) = 0
comm( 8) = 0
!     COMM( 9) = 0
!     COMM(12) = 0
!     COMM(15) = 0
!     COMM(18) = 0
DO  j = 1,1100
  xf1at(j) = 0
END DO
GO TO 400

340 CONTINUE
350 CALL mesage (-37,0,NAME)
360 CALL mesage (-1,FILE,NAME)

400 RETURN
END SUBROUTINE endsys
