SUBROUTINE tmtsio (*,debug1)
     
!     TMTSIO TIME TESTS GINO AND THE PACK ROUTINES
 
!     COMMENT FORM G.CHAN/UNISYS   5/91
!     BASICALLY THIS ROUTINE IS SAME AS TIMTS1.
 
 
 INTEGER, INTENT(IN OUT)                  :: debug1
 INTEGER :: sysbuf, output, files(2), f1 , f2, buf1  , buf2  ,  &
            end   , mcb(7), eol     , eor    , type  , typin1,  &
            typou1, typou2, ablk(15), bblk(15),zero  , isubr(2)
REAL :: x(1)  , z(1)  , t(23)
DOUBLE  PRECISION       zd    , xd
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /xmssg / ufm   , uwm   , uim     , sfm
COMMON /ntime / nitems, tgino , tbldpk  , tintpk , tpack ,  &
    tunpak, tgetst, tputst  , ttlrsp, ttlrdp, ttlcsp  , ttlcdp ,  &
    tllrsp, tllrdp, tllcsp  , tllcdp , tgetsb  &
    , rgino , rbldpk  , rintpk , rpck   , runpak, rgetst  , rputst
COMMON /system/ sysbuf, output, idum(52), iprec  , jdum(21), isy77
COMMON /ginox / g(86) , ig(75), nbuff3  , pu(226), iwr(75)
COMMON /zblpkx/ zd(2),  iz
COMMON /zntpkx/ xd(2),  ix    , eol     , eor
COMMON /packx / typin1, typou1, i1 , j1 , incr1
COMMON /unpakx/ typou2, i2, j2, incr2
COMMON /zzzzzz/ a(1)
EQUIVALENCE     (zd(1),z(1)), (xd(1),x(1)), (t(1),tgino)
DATA    files / 301, 304/, zero  / 0    /
DATA    i1000 / 1000    /, i1001 / 1001 /
DATA    isubr / 4HTMTS   , 4HIO  /


!     CHECK KORSZ AND DUMMY SUBROUTINES HERE.
!     IF NASTRAN SUBROUTINES WERE NOT COMPILED WITH STATIC OPTION, BUF1
!     COULD BE NEGATIVE HERE.
!     CALL DUMMY NEXT TO SEE WHETHER THE RIGHT DUMMY ROUTINE IS SET UP
!     FOR THIS MACHINE

kerr = 1
buf1 = korsz(a)
IF (buf1 <= 0) GO TO 930
IF (debug1 > 0) WRITE (output,10)
10 FORMAT (' -LINK1 DEBUG- TMTSIO CALLINS DUMMY NEXT')
!WKBD CALL DUMMY

!     NOTE - ISY77 (WHICH IS BULKDATA OPTION) AND TGINO DETERMINE TO
!            SKIP TMTSIO AND TMTSLP OR NOT. DIAG 35 CAN NOT BE USED AT
!            THIS POINT SINCE THE DIAG CARD HAS NOT BEEN READ YET.

IF (tgino > 0. .AND. isy77 /= -3) RETURN 1

!     INITIALIZE

CALL page1
WRITE  (output,20)
20 FORMAT ('0*** USER INFORMATION MESSAGE 225, GINO TIME CONSTANTS ',  &
    'ARE BEING COMPUTED', /5X,  &
    '(SEE NASINFO FILE FOR ELIMINATION OF THESE COMPUTATIONS)')
IF (tgino > 0.) WRITE (output,30) t
30 FORMAT ('0*** EXISTING TIME CONSTANTS IN /NTIME/ -', /,2(/5X,9F8.3))
n = 50
m = n
TYPE = iprec
!     NITEMS = 23

f1   = files(1)
f2   = files(2)
buf1 = buf1 - sysbuf
buf2 = buf1 - sysbuf
END  = n*m
IF (END >= buf1-1) CALL mesage (-8,0,isubr)
DO  i = 1,END
a(i) = i
END DO
n10 = n*10
m10 = m/10
IF (m10 <= 0) m10 = 1
fn  = n
fm  = m

!     WRITE TEST

IF (debug1 > 0) WRITE (output,50) nbuff3,ig
50 FORMAT (' -LINK1 DEBUG- OPEN OUTPUT FILE NEXT FOR WRITE. NBUFF3 ='  &
    ,       i5, /5X,'GINO BUFADD 75 WORDS =', /,(2X,11I7))
CALL OPEN (*900,f1,a(buf1),1)
IF (debug1 <= 0) GO TO 60
WRITE  (output,53) nbuff3,ig
53 FORMAT (' -LINK1 DEBUG- FILE OPEN OK. NBUFF3 =',i5, /5X,  &
    'GINO BUFADD 75 WORDS =', /,(2X,11I7))
WRITE  (output,55) iwr(41)
55 FORMAT (5X,'RWFLG(41) =',i7, //, ' -LINK1 DEBUG- CALLING SECOND NEXT')
60 CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL WRITE (f1,a,m,1)
END DO
CALL cputim (t2,t2,1)
IF (debug1 > 0) WRITE (output,80)
80 FORMAT (' -LINK1 DEBUG- CLOSE FILE NEXT')
CALL CLOSE  (f1,1)
IF (debug1 > 0) WRITE (output,90)
90 FORMAT (' -LINK1 DEBUG- OPEN ANOTHER OUTPUT FILE NEXT FOR WRITE')
CALL OPEN   (*900,f2,a(buf2),1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL WRITE (f2,a,m10,1)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 610 TO iret
GO TO 600

!     READ TEST

110 IF (debug1 > 0) WRITE (output,120)
120 FORMAT (' -LINK1 DEBUG- OPEN INPUT FILE NEXT FOR READ')
CALL OPEN (*900,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL READ (*910,*920,f1,a(i1000),m,1,flag)
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,2)
CALL OPEN   (*900,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL READ (*910,*920,f2,a(i1000),m10,1,flag)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,2)
ASSIGN 620 TO iret
GO TO 600

!     BACKWARD READ TEST

150 CONTINUE
CALL OPEN (*900,f1,a(buf1),2)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL bckrec (f1)
  CALL READ   (*910,*920,f1,a(i1000),m,1,flag)
  CALL bckrec (f1)
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
CALL OPEN   (*900,f2,a(buf2),2)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL bckrec (f2)
  CALL READ   (*910,*920,f2,a(i1000),m10,1,flag)
  CALL bckrec (f2)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 630 TO iret
GO TO 600

!     BLDPK TEST

180 CONTINUE
CALL OPEN   (*900,f1,a(buf1),1)
CALL makmcb (mcb,f1,m,2,TYPE)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL bldpk (TYPE,TYPE,f1,0,0)
  DO  j = 1,m
    z(1) = 1.0
    iz   = j
    CALL zblpki
  END DO
  CALL bldpkn (f1,0,mcb)
END DO
CALL cputim (t2,t2,1)
CALL wrttrl (mcb)
CALL CLOSE  (f1,1)
CALL makmcb (mcb,f2,m10,2,TYPE)
CALL OPEN   (*900,f2,a(buf2),1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL bldpk (TYPE,TYPE,f2,0,0)
  DO  j = 1,m10
    z(1) = 2.0
    iz   = j
    CALL zblpki
  END DO
  CALL bldpkn (f2,0,mcb)
END DO
CALL cputim (t4,t4,1)
CALL wrttrl (mcb)
CALL CLOSE  (f2,1)
ASSIGN 640 TO iret
GO TO 600

!     INTPK TEST

230 CONTINUE
CALL OPEN   (*900,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL intpk  (*910,f1,0,TYPE,0)
  DO  j = 1,m
    CALL zntpki
    IF (ix  /= j) GO TO 800
    IF (eol == 0) CYCLE
    IF (ix  /= m) GO TO 800
  END DO
  IF (eol == 0) GO TO 800
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
CALL OPEN   (*900,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL intpk (*910,f2,0,TYPE,0)
  DO  j = 1,m10
    CALL zntpki
    IF (ix  /= j) GO TO 800
    IF (eol == 0) CYCLE
    IF (ix  /= m10) GO TO 800
  END DO
  IF (eol == 0) GO TO 800
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 650 TO iret
GO TO 600

!     PACK TEST

280 CONTINUE
CALL makmcb (mcb,f1,m,2,TYPE)
typin1 = TYPE
typou1 = TYPE
i1     = 1
j1     = m
incr1  = 1
mx     = m*TYPE
DO  i = 1,mx
  a(i+1000) = i
END DO
CALL OPEN (*900,f1,a(buf1),1)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL pack (a(i1001),f1,mcb)
END DO
CALL cputim (t2,t2,1)
CALL wrttrl (mcb)
CALL CLOSE  (f1,1)
CALL makmcb (mcb,f2,m10,2,TYPE)
j1 = m10
CALL OPEN (*900,f2,a(buf2),1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL pack (a(i1001),f2,mcb)
END DO
CALL cputim (t4,t4,1)
CALL wrttrl (mcb)
CALL CLOSE  (f2,1)
ASSIGN 660 TO iret
GO TO 600

!     UNPACK TEST

320 CONTINUE
typou2 = TYPE
i2     = 1
j2     = m
incr2  = 1
CALL OPEN (*900,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL unpack (*910,f1,a(i1001))
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
j2 = m10
CALL OPEN (*900,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL unpack (*910,f2,a(i1001))
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,2)
ASSIGN 670 TO iret
GO TO 600
350 CONTINUE

!     PUTSTR TEST

kerr    = 2
ablk(1) = f1
ablk(2) = TYPE
ablk(3) = 1
CALL gopen (f1,a(buf1),1)
nwds = TYPE
IF (TYPE == 3) nwds = 2
CALL cputim (t1,t1,1)
DO  i = 1,n
  ablk(4) = 0
  ablk(8) = -1
  DO  j = 1,10
    nbrstr  = m10
    360 CALL putstr (ablk)
    IF (nbrstr == 0) GO TO 930
    ablk(7) = MIN0(ablk(6),nbrstr)
    ablk(4) = ablk(4) + ablk(7) + 4
    mm      = ablk(7)*nwds
    DO  k = 1,mm
      x(1) = a(k)
    END DO
    IF (ablk(7) == nbrstr) GO TO 380
    CALL endput (ablk)
    nbrstr = nbrstr - ablk(7)
    GO TO 360
    380 IF (j == 10) ablk(8) = 1
    CALL endput (ablk)
  END DO
END DO
CALL cputim (t2,t2,1)
CALL CLOSE (f1,1)
m100 = MAX0(m10/10,1)
CALL gopen (f2,a(buf2),1)
kerr    = 3
bblk(1) = f2
bblk(2) = TYPE
bblk(3) = 1
CALL cputim (t3,t3,1)
DO  i = 1,n10
  bblk(4) = 0
  bblk(8) =-1
  DO  j = 1,10
    nbrstr = m100
    410 CALL putstr (bblk)
    IF (nbrstr == 0) GO TO 930
    bblk(7) = MIN0(bblk(6),nbrstr)
    bblk(4) = bblk(4) + bblk(7) + 4
    mm = bblk(7)*nwds
    DO  k = 1,mm
      x(1) = a(k)
    END DO
    IF (bblk(7) == nbrstr) GO TO 430
    nbrstr = nbrstr - bblk(7)
    GO TO 410
    430 IF (j == 10) bblk(8) = 1
    CALL endput (bblk)
  END DO
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 680 TO iret
GO TO 600

!     GETSTR TEST (GET STRING FORWARD)

460 CONTINUE
CALL gopen (f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  ablk(8) = -1
  470 CALL getstr (*490,ablk)
  mm = ablk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endget (ablk)
  GO TO 470
  490 CONTINUE
END DO
CALL cputim (t2,t2,1)
!     CALL CLOSE  (F1,1)
CALL gopen  (f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  bblk(8) = -1
  500 CALL getstr (*520,bblk)
  mm = bblk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endget (bblk)
  GO TO 500
  520 CONTINUE
END DO
CALL cputim (t4,t4,1)
!     CALL CLOSE  (F2,1)
ASSIGN 690 TO iret
GO TO 600

!     GETSTB TEST, (GET BACKWARD STRING)
!     F1 AND F2 FILES ARE STILL OPENED, AND POSITIONED AT THE END

530 CONTINUE
!     CALL GOPEN (F1,A(BUF1),0)
!     CALL REWIND (F1)
!     CALL SKPFIL (F1,N+1)
CALL cputim (t1,t1,1)
DO  i = 1,n
  ablk(8) = -1
  540 CALL getstb (*560,ablk)
  mm = ablk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endgtb (ablk)
  GO TO 540
  560 CONTINUE
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
!     CALL GOPEN  (F2,A(BUF2),0)
!     CALL REWIND (F2)
!     CALL SKPFIL (F2,N10+1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  bblk(8) = -1
  570 CALL getstb (*590,bblk)
  mm = bblk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endgtb (bblk)
  GO TO 570
  590 CONTINUE
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 700 TO iret

!     INTERNAL ROUTINE TO STORE TIMING DATA IN /NTIME/ COMMON BLOCK

600 CONTINUE
time1  = t2 - t1
time2  = t4 - t3
tprrec = 1.0E6*(time2 - time1)/(9.0*fn)
tprwrd = (1.0E6*time1 - fn*tprrec)/(fn*fm)
GO TO iret, (610,620,630,640,650,660,670,680,690,700)
610 tgino  = tprwrd
rgino = tprrec
GO TO 110
620 tgino  = tgino + tprwrd
rgino  = rgino + tprrec
GO TO 150
630 tgino  = tgino + tprwrd
tgino  = tgino/3.0
rgino  = rgino + tprrec
rgino  = rgino/3.0
GO TO 180
640 tbldpk = tprwrd
rbldpk = tprrec
GO TO 230
650 tintpk = tprwrd
rintpk = tprrec
GO TO 280
660 tpack  = tprwrd
rpack  = tprrec
GO TO 320
670 tunpak = tprwrd
runpak = tprrec
GO TO 350
680 tputst = tprwrd
rputst = tprrec
GO TO 460
690 tgetst = tprwrd
rgetst = tprrec
GO TO 530
700 tgetsb = tprwrd
IF (debug1 > 0) WRITE (output,710)
710 FORMAT (' -LINK1 DEBUG- TMTSIO FINISHED')
RETURN

!     INTERNAL ROUTINE CALLED FOR AN ABORT IN THE INTPK TEST

800 WRITE  (output,810) sfm
810 FORMAT (a25,' 2197, ABORT CALLED DURING TIME TEST OF INTPK')

!     ABNORMAL RETURNS FROM GINO - ALL FATAL ERRORS

900 CONTINUE
910 CONTINUE
920 CALL mesage (-61,0,0)
930 WRITE  (output,940) kerr
940 FORMAT ('0*** TMTSIO FATAL ERROR',i7)
GO TO 920

END SUBROUTINE tmtsio
