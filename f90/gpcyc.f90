SUBROUTINE gpcyc
     
!     GPCYC IS THE GEOMETRY PROCESSOR FOR CYCLIC PROBLEM
 
!     INPUT DATA BLOCKS - GEOM4,EQEXIN,USET
 
!     OUTPUT DATA BLOCKS - CYCD
 
!     PARAMETERS  CTYPE - INPUT,BCD -
!                 NOGO  - OUTPUT--+1 UNLESS ERROR--THEN-1
 
!     SCRATCH FILES (2)
!     DEFINITION OF VARIABLES
!     NZ       OPEN CORE LENGTH
!     NX       ORIGINAL OPEN CORE
!     NENT     NUMBER OF ENTRIES IN EQEXIN
!     ITYP     PROBLEM TYPE (ROT=0 ,OTHERWISE=1)
!     LCYJ     LENGTH OF CJOIN CARDS
!     ISID1    POINTER TO START OF SIDE 1 CZRDS
!     ISID2    POINTER TO START OF SIDE 2 CZRDS
 
 EXTERNAL        andf
 INTEGER :: geom4,eqexin,uset,cycd,ctype,sysbuf,FILE,NAME(2),  &
     scr1,scr2,REC,cyl,sph,rot,cyjoin(2),ib(5),andf, ibb(4),mcb(7),blk
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /two   / itwo(32)
 COMMON /bitpos/ isk(6),iua
 COMMON /system/ ksystm(65)
 COMMON /zzzzzz/ iz(1)
 COMMON /BLANK / ctype(2),nogo
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(2),nout)
 DATA    geom4 , eqexin,uset,cycd,scr1,scr2,NAME         /  &
     101   , 102   ,103 ,201 ,301 ,302 ,4HGPCY,4HC   /
 DATA    rot   / 4HROT /,REC,cyl,sph  / 1HR,1HC,1HS      /
 DATA    cyjoin/ 5210,52  /
 DATA    nocy  , nosid1, isid1, iblen, icm,  isam, nocnt, nopar /  &
     4024  , 4025,   4026,  4027,  4028, 4029, 4030,  4032  /
 DATA    noeq  , ncord    / 4037  , 4039     /
 DATA    mcb   / 7*0/, blk/ 1H  /
 
 
 nz   = korsz(iz)
 nogo = 1
 
!     IBUF1 IS PRELOC BUFFER
 
 ibuf1 = nz - sysbuf
 ibuf2 = ibuf1 - sysbuf
 nz = ibuf2 - 1
 nx = nz
 IF (nz <= 0) CALL mesage (-8,0,NAME)
 
!     PUT  SECOND RECORD OF EQEXIN ITO CORE
 
 FILE = eqexin
 CALL gopen (eqexin,iz(ibuf1),0)
 CALL fwdrec (*560,eqexin)
 CALL READ (*560,*10,eqexin,iz,nz,0,iflag)
 CALL mesage (-8,0,NAME)
 10 CALL CLOSE (eqexin,1)
 nent = iflag/2
 
!     DECIDE ON TYPE
 
 ityp = 1
 IF (ctype(1) == rot) ityp = 0
 
!     FIND  CYJOIN CARDS ON GEOM4
 
 FILE = geom4
 CALL preloc (*540,iz(ibuf1),geom4)
 CALL locate (*580,iz(ibuf1),cyjoin,idum)
 nz = nz - iflag
 k  = iflag + 1
 CALL READ (*560,*20,geom4,iz(k),nz,0,lcyj)
 CALL mesage (-8,0,NAME)
 20 CALL CLOSE (geom4,1)
 lcyj = lcyj + k - 1
 IF (iz(k) == 1) GO TO 40
 WRITE  (nout,590) ufm,nosid1
 WRITE  (nout,30)
 30 FORMAT ('0NO SIDE 1 DATA FOUND.')
 GO TO 620
 31 WRITE  (nout,590) ufm,nosid1
 WRITE  (nout,32)
 32 FORMAT ('0NO SIDE 2 DATA FOUND.')
 GO TO 620
 
!     FIND SIDE TWO DATA
 
 40 l = k
 50 IF (l > lcyj) GO TO  31
 IF (iz(l) == -1) GO TO 70
 60 l = l + 1
 GO TO 50
 
!     END OF CARD FOUND
 
 70 IF (l+1  > lcyj) GO TO 31
 IF (iz(l+1) == 2) GO TO 90
 IF (ityp    == 1) GO TO 60
 WRITE  (nout,590) ufm,isid1
 WRITE  (nout,80)
 80 FORMAT ('0TOO MANY SIDE 1 CARDS.')
 GO TO 620
 
!     FOUND SIDE TWO LIST
 
 90 isid2 = l + 1
 IF (ityp /= 0) GO TO 370
 
!     CHECK LENGTH OF SIDE TWO LIST
 
 ns1 = isid2 - k - 4
 ns2 = lcyj - isid2 - 3
 IF (ns1 == ns2) GO TO 110
 WRITE  (nout,590) ufm,iblen
 WRITE  (nout,100)
 100 FORMAT ('0NUMBER OF ENTRIES IN SIDE 1 NOT EQUAL TO NUMBER IN ', 'SIDE 2')
 nogo = -1
 GO TO 620
 
!     BUILD 5 WORDS FOR EACH PAIR
 
 
!     FIVE WORD ENTRY FOR EACH PAIR APPEARS AS FOLLOWS
 
!     1        CODE(1 = GRID   2 = SCALAR)
!     2        INTERNAL INDEX (SIL)      SIDE 1
!     3        GRID ID (EXTERNAL)       SIDE 1
!     4        INTERNAL INDEX (SIL)      SIDE 2
!     5        GRID ID (EXTERNAL)       SIDE 2
 
 110 CALL gopen (scr1,iz(ibuf1),1)
 l = isid2 + 3
 k = k + 3
 DO  i = 1,ns1
   IF (iz(k) /= iz(l)) GO TO  130
   WRITE  (nout,590) ufm,isam
   WRITE  (nout,120) iz(k)
   120 FORMAT ('0GRID POINT',i10,' APPEARS IN BOTH SIDE LISTS.')
   GO TO 620
   130 CONTINUE
   ip = iz(k)
   CALL bisloc (*610,ip,iz(1),2,nent,m)
   ix1 = iz(m+1)/10
   ic1 = iz(m+1) - ix1*10
   ip  = iz(l)
   CALL bisloc (*610,ip,iz(1),2,nent,m)
   ix2 = iz(m+1)/10
   ic2 = iz(m+1) - ix2*10
   IF (ic1 == ic2) GO TO 150
   WRITE  (nout,590) ufm,icm
   WRITE  (nout,140) iz(k),iz(l)
   140 FORMAT ('0THE CODE FOR GRID POINT',i10,' DOES NOT MATCH THE CODE',  &
       ' FOR GRID POINT',i10)
   GO TO 620
   150 ib(1) = ic1
   ib(2) = ix1
   ib(3) = iz(k)
   ib(4) = ix2
   ib(5) = iz(l)
   CALL WRITE (scr1,ib,5,0)
   k = k + 1
   l = l + 1
 END DO
 170 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,1)
 
!     SET UP USET
 
 CALL gopen (uset,iz(ibuf1),0)
 FILE = uset
 nz   = nx
 CALL READ (*560,*190,uset,iz,nz,0,luset)
 CALL mesage (-8,0,NAME)
 190 CALL CLOSE (uset,1)
 
!     SET UP REDUCED USET TABLE
 
 k = 0
 m = itwo(iua)
 DO  i = 1,luset
   IF (andf(iz(i),m) == 0.0) THEN
     GO TO   200
   ELSE
     GO TO   210
   END IF
   200 iz(i) = 0
   CYCLE
   210 k = k + 1
   iz(i) = -k
 END DO
 lua = k
 
!     FORM SILA VALUES
 
 FILE = scr1
 CALL gopen (scr1,iz(ibuf1),0)
 CALL gopen (scr2,iz(ibuf2),1)
 IF (ityp /= 0) GO TO  410
 230 CALL READ (*560,*300,scr1,ib(1),5,0,iflag)
 np = 1
 IF = 0
 IF (ib(1) == 1) np = 6
 k = 0
 240 l = ib(2) + k
 m = ib(4) + k
 
!     IF NEITHER IGNORE
 
 IF (iz(l) == 0 .AND. iz(m) == 0) GO TO 280
 IF (iz(l) < 0 .AND. iz(m) < 0) GO TO 270
 WRITE  (nout,250) uwm,nocnt
 250 FORMAT (a25,i5)
 m = k + 1
 WRITE  (nout,260) m,ib(3),ib(5)
 260 FORMAT ('0COMPONENT',i4,' OF GRID POINTS',i10,5H AND ,i10,  &
     ' CANNOT BE CONNECTED.')
 GO TO 280
 270 IF = IF + 1
 ibb(1) = IABS(iz(l))
 ibb(2) = ib(3)
 ibb(3) = IABS(iz(m))
 ibb(4) = ib(5)
 CALL WRITE (scr2,ibb,4,0)
 280 k = k + 1
 IF (k /= np) GO TO 240
 IF (IF /= 0) GO TO 230
 WRITE  (nout,250) uwm,nopar
 WRITE  (nout,290) ib(3),ib(5)
 290 FORMAT ('0NO COMPONENTS OF GRID POINTS',i10,5H AND ,i10,  &
     ' WERE CONNECTED.')
 GO TO 230
 
!     CLOSE UP
 
 300 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr1,1)
 CALL CLOSE (scr2,1)
 
!     BUILD CYCD
 
 DO  i = 1,lua
   iz(i) = 0
 END DO
 FILE = scr2
 CALL gopen (scr2,iz(ibuf1),0)
 IF (ityp /= 0) GO TO 520
 320 CALL READ (*560,*360,scr2,ibb,4,0,iflag)
 k = ibb(1)
 m = ibb(3)
 IF (iz(k) == 0) GO TO 340
 WRITE  (nout,590) ufm,noeq
 WRITE  (nout,330) ibb(2)
 330 FORMAT ('0GRID POINT',i10,' IS LISTED MORE THAN ONCE.')
 nogo = -1
 340 IF (iz(m) == 0) GO TO 350
 WRITE (nout,590) ufm,noeq
 WRITE (nout,330) ibb(4)
 nogo  = -1
 350 iz(k) =  m
 iz(m) = -k
 GO TO  320
 
!     END OF PAIRS
 
 360 CALL CLOSE (scr2,1)
 CALL gopen (cycd,iz(ibuf1),1)
 CALL WRITE (cycd,iz(1),lua,1)
 CALL CLOSE (cycd,1)
 mcb(1) = cycd
 mcb(2) = ityp + 1
 mcb(3) = lua
 CALL wrttrl (mcb)
 IF (nogo /= -1) RETURN
 GO TO 620
 
!     1. DIHEDRAL TYPE
 
!     BUILD FIVE WORD LIST
 
 
!     FIVE WORD ENTRY FOR EACH POINT IN SIDE 1 OR SIDE TWO LOOKS AS
!         FOLLOWS
!     1        SIDE (1,2)
!     2        COORD SYS (R = 1,C = 1,S = 2,BLANK = 0)
!     3        CODE ( 1 = GRID   2 = SCALAR)
!     4        INTERNAL INDEX (SIL)
!     5        GRID ID (EXTERNAL)
 
 370 l = k
 CALL gopen (scr1,iz(ibuf1),1)
 380 icid = iz(l+1)
 isid = iz(l  )
 IF (icid == REC) icid = 1
 IF (icid == cyl) icid = 1
 IF (icid == sph) icid = 2
 IF (icid == blk) icid = 0
 l = l + 3
 390 IF (iz(l) == -1) GO TO 400
 ip = iz (l)
 CALL bisloc (*610,ip,iz(1),2,nent,m)
 ib(1) = isid
 ib(2) = icid
 ib(4) = iz(m+1)/10
 ib(3) = iz(m+1) - ib(4)*10
 ib(5) = ip
 CALL WRITE (scr1,ib,5,0)
 l = l + 1
 GO TO 390
 
!     END OF LIST
 
 400 IF (l >= lcyj) GO TO 170
 l = l + 1
 GO TO 380
 
!     END OF CYJOIN LISTS
 
 
!     PRODUCE CYCD CODES
 
 410 CALL READ (*560,*300,scr1,ib(1),5,0,iflag)
 np = 1
 IF (ib(3) == 1) np = 6
 IF = 0
 k  = 0
 IF (ib(1) == 2) ib(1) = ib(1) + 1
 420 l = ib(4) + k
 IF (iz(l) == 0) GO TO 500
 
!     POINT IS IN  A SET
 
 ibb(2) = IABS(iz(l))
 ibb(3) = ib(5)
 IF (ib(3) == 2) GO TO 480
 IF (ib(2) == 1) GO TO 440
 IF (ib(2) == 2) GO TO 460
 
!     COORD SYS = 0
 
 WRITE  (nout,590) ufm,ncord
 WRITE  (nout,430) ibb(3)
 430 FORMAT ('0NO COORDINATE SYSTEM DEFINED FOR GRID POINT',i10)
 nogo = -1
 GO TO 480
 
!     RECTANGULAR OR CYL
 
 440 IF (MOD(k+1,2) == 1) GO TO 480
 450 m = 1
 GO TO 490
 
!     SPH
 
 460 IF (k < 2 .OR. k == 5 .OR. np < 3 .OR. np == 6) GO TO 480
 GO TO 450
 
!     EVEN
 
 480 m = 0
 490 ibb(1) = ib(1) + m
 IF = IF + 1
 CALL WRITE (scr2,ibb,3,0)
 500 k = k + 1
 IF (k /= np) GO TO 420
 IF (IF /= 0) GO TO 410
 WRITE  (nout,250) uwm,nopar
 WRITE  (nout,510) ib(5)
 510 FORMAT ('0NO COMPONENTS OF GRID POINT',i10,' WERE IN THE A SET')
 GO TO 410
 
!     BUILD CYCD FOR DIH
 
 520 CALL READ (*540,*360,scr2,ibb,3,0,iflag)
 k = ibb(2)
 IF (iz(k) == 0) GO TO 530
 WRITE (nout,590) ufm,noeq
 WRITE (nout,330) ibb(3)
 nogo  = -1
 530 iz(k) = ibb(1)
 GO TO 520
 
!     ERROR MESSAGES
 
 540 ip1 = -1
 550 CALL mesage (ip1,FILE,NAME)
 RETURN
 560 ip1 = -2
 GO TO 550
 580 WRITE  (nout,590) ufm,nocy
 590 FORMAT (a23,i5)
 WRITE  (nout,600)
 600 FORMAT ('0NO CYJOIN CARDS WERE SUPPLIED.')
 GO TO 620
 610 CALL mesage (-30,2,ip)
 620 CALL mesage (-61,0,NAME)
 RETURN
END SUBROUTINE gpcyc
