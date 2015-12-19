SUBROUTINE frbk2 (v1,v2,v3,vb)
     
!     LAST REVISED BY G.CHAN/UNISYS  11/1991
!     . ELIMINATE UN-NECCESSARY REWIND AND SKIP AFTER FIRST CALL TO THIS
!       ROUTINE (NASTRAN ORIGINAL METHOD)
!     . ADDITION OF A NEW BACKWARD-FORWARD SUBSTITUTION METHOD WHICH IS
!       MORE EFFICIENT, AND IS ALREADY GOOD FOR VECTORIZATION
 
!DB   LOGICAL          DEBUG
 
 DOUBLE PRECISION, INTENT(IN)             :: v1(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v2(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v3(1)
 DOUBLE PRECISION, INTENT(IN)             :: vb(1)
 INTEGER :: base    ,buf(6) ,ljj(2)  ,iblk(15)
 DOUBLE PRECISION :: xl(1)  ,xljj , zero    ,v3j    ,sum
 COMMON  /opinv / mcblt(7),mcbsma(7)
 COMMON  /system/ ksystm  ,io
 COMMON  /feerxx/ dumm(18),nzvb
 COMMON  /zzzzzz/ iz(1)
 EQUIVALENCE      (xl(1),iz(1))
 EQUIVALENCE      (ljj(1) ,xljj)  ,(l16,dumm(6))
DATA     buf   / 4HFRBK  ,4H2    ,2*4HBEGN ,4HEND  ,4HBGIN /
DATA     zero  / 0.0D+0       /
!DB   DATA     DEBUG , ITER    ,MAX /  .FALSE.   ,0      ,3      /

!DB   IF (.NOT.DEBUG) GO TO 20
!     ITER = ITER + 1
!     IF (ITER .GT. MAX) DEBUG = .FALSE.
!     WRITE  (IO,10) NZVB,ITER
!  10 FORMAT ('  .... IN FRBK2.  NZVB =',I8,',  ITER =',I3)
!  20 CONTINUE
nrow = mcblt(2)
DO  i = 1,nrow
  v2(i) = v1(i)
END DO

!     SELECTION OF ORIGINAL OR NEW FBS METHOD

j = nrow
IF (mcblt(7) < 0) GO TO 200

!     NASTRAN ORIGIANL METHOD

iblk( 1) = mcblt(1)
iblk( 9) = 1
iblk(10) = 1

!     BACKWARD SUBSTITUTION

IF (buf(3) == buf(5)) GO TO 40
!     BUF(3) = BUF(4)
IF (l16 /= 0) CALL conmsg (buf,3,0)

!     REWIND AND SKIP TO COLUMN N

CALL REWIND (mcblt)
CALL skprec (mcblt,nrow+1)
GO TO 50

!     ALREADY AT END, NO SKIP NEEDED

40 buf(3) = buf(6)
IF (l16 /= 0) CALL conmsg (buf,3,0)

50 iblk(8) = -1
60 CALL getstb (*100,iblk(1))
ntms = iblk(6)
ji   = iblk(5)
ik   = iblk(4)
IF (ik-ntms+1 /= j) GO TO 70
ntms = ntms - 1
xljj = xl(ji-ntms)
IF (ntms == 0) GO TO 90
70 sum  = zero
DO  ii = 1,ntms
  sum  = sum + xl(ji)*v2(ik)
  ji   = ji - 1
  ik   = ik - 1
END DO
v2(j)= v2(j) + sum
90 CALL endgtb (iblk(1))
GO TO 60
100 v2(j)= v2(j)/xljj
IF (j == 1) GO TO 110
j    = j - 1
GO TO 50
110 CALL frmltd (mcbsma(1),v2(1),v3(1),vb(1))

!     FORWARD SWEEP DIRECTLY ON V3

DO  j = 1,nrow
  iblk(8) = -1
  120 CALL getstr (*160,iblk(1))
  ji   = iblk(5)
  ntms = iblk(6)
  ik   = iblk(4)
  IF (ik /= j) GO TO 130
  ntms = ntms - 1
  v3(j)= v3(j)/xl(ji)
  ji   = ji + 1
  ik   = ik + 1
  130 IF (ntms == 0) GO TO 150
  v3j  = v3(j)
  IF (v3j == zero) GO TO 150
  DO  ii = 1,ntms
    v3(ik) = v3(ik) + xl(ji)*v3j
    ik   = ik + 1
    ji   = ji + 1
  END DO
  150 CALL endget (iblk(1))
  GO TO 120
  160 CONTINUE
END DO
GO TO 500

!     NEW METHOD

!     MATRIX MCBLT HAS BEEN RE-WRITTEN TO MCBLTX BY UNPSCR/FEER3. NO
!     STRING OPERATIONS HERE.

200 IF (buf(3) == buf(5)) buf(3) = buf(6)
IF (l16 /= 0) CALL conmsg (buf,3,0)
mcbltx = -mcblt(7)
IF (MOD(mcblt(4),10) /= 3) GO TO 440
nrec = 0
CALL REWIND (mcbltx)
CALL fwdrec (*400,mcbltx)
nwds = mcblt(5)

!     IZ(1)                                                      GINO
!      / V1   V2    V3           VB (OPEN CORE LENGTH = NZVB)   BUFFERS
!     +-----+-----+-----+-----+-------------------------------+---------
!                          OPEN  CORE


!     BACKWARD SUBSTITUTION


ll2  = 0
base = 1
ifb  = -350
DO  ik = 1,nrow
  IF (base < ll2) GO TO 240
  nrec = nrec + 1
!DB   IF (DEBUG) WRITE (IO,210) NREC,IFB
! 210 FORMAT ('  ...READING RECORD',I5,'.   IFB =',I5)
  CALL READ (*400,*220,mcbltx,vb,nzvb,1,ll)
  CALL mesage (-8,0,nam)
  220 ll2  = ll/nwds
!DB   LL3  = LL2/30
!     LL4  = LL2 - LL3
!     IF (DEBUG) WRITE (IO,230) LL,NREC,LL2
! 230 FORMAT (1X,I10,' WORDS READ FROM RECORD NO.',I5,'.   LL2 =',I10)
  base = 1
  240 xljj = vb(base)
  ii   = ljj(1)
  jj   = ljj(2)
!DB   IF (DEBUG .AND. (BASE.LT.LL3 .OR. BASE.GT.LL4))
!    1    WRITE (IO,250) J,BASE,II,JJ,IFB
! 250 FORMAT (11X,'J,BASE,II,JJ,IFB =',5I8)
  IF (ii /= j) GO TO 420
  ntms = jj - ii + 1
  ib   = base + 2
  ie   = base + ntms
  base = ie + 1
  IF (ntms <= 1) GO TO 270
  sum  = zero
  DO  i = ib,ie
    ii   = ii + 1
    sum  = sum + vb(i)*v2(ii)
  END DO
  v2(j)= v2(j) + sum
  270 v2(j)= v2(j)/vb(ib-1)
  j    = j - 1
END DO
CALL frmltd (mcbsma(1),v2(1),v3(1),vb(1))

!     FORWARD SWEEP DIRECTLY ON V3

IF (nrow == 1) GO TO 500
nrec = 0
ll2  = 0
base = 1
ifb  = +390
DO  j = 1,nrow
  IF (base < ll2) GO TO 300
  nrec = nrec + 1
!DB   IF (DEBUG) WRITE (IO,210) NREC,IFB
  CALL READ (*400,*290,mcbltx,vb,nzvb,1,ll)
  CALL mesage (-8,0,nam)
  290 ll2  = ll/nwds
!DB   LL3  = LL2/30
!     LL4  = LL2 - LL3
!     IF (DEBUG) WRITE (IO,230) LL,NREC,LL2
  base = 1
  300 xljj = vb(base)
  ii   = ljj(1)
  jj   = ljj(2)
!DB   IF (DEBUG .AND. (BASE.LT.LL3 .OR. BASE.GT.LL4))
!    1    WRITE (IO,250) J,BASE,II,JJ,IFB
  IF (ii /= j) GO TO 420
  ntms = jj - ii + 1
  v3(j)= v3(j)/vb(base+1)
  IF (ntms <= 1) GO TO 320
  v3j  = v3(j)
  IF (v3j == zero) GO TO 320
  ib   = base + 2
  ie   = base + ntms
  DO  i = ib,ie
    ii   = ii + 1
    v3(ii) = v3(ii) + vb(i)*v3j
  END DO
  320 base = base + ntms + 1
END DO
GO TO 500

400 i = mcblt(4)/10
WRITE  (io,410) nrec,j,i,ifb
410 FORMAT ('0*** TRY TO READ RECORD',i5,'.  J,MCBLT(4),IFB =',i7,2I5)
CALL mesage (-3,mcbltx,nam)
420 WRITE  (io,430) j,ii,ifb
430 FORMAT ('0*** ROW MISMATCH.  J,II,(IFB =',i7,i12,3H  (,i4)
GO TO 460
440 j = MOD(mcblt(4),10)
WRITE  (io,450) j
450 FORMAT ('0*** MCBLT MATRIX IN WRONG FORM.   UNPSCR FLAG =',2I3)
460 CALL mesage (-37,0,buf(1))

500 buf(3) = buf(5)
IF (l16 /= 0) CALL conmsg (buf,3,0)
RETURN
END SUBROUTINE frbk2
