SUBROUTINE frsw2 (v1,v2,v3,vb)
     
!     LAST REVISED  11/91, BY G.CHAN/UNISYS
!     ADDITION OF A NEW FORWARD-BACKWARD SUBSTITUTION METHOD, WHICH IS
!     MORE EFFICIENT, AND IS ALREADY GOOD FOR VECTORIZATION.
 
!DB   LOGICAL          DEBUG
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: v1(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v2(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v3(1)
 DOUBLE PRECISION, INTENT(IN)             :: vb(1)
 INTEGER :: nam(6)  ,ljj(2) ,iblk(15),base
 DOUBLE PRECISION :: xl(1) ,xljj  , v3j     ,zero   ,sum
 COMMON  /opinv / mcblt(7),mcbsma(7)
 COMMON  /system/ ksystm  ,io
 COMMON  /feerxx/ dumm(18),nzvb
 COMMON  /zzzzzz/ iz(1)
 EQUIVALENCE      (xl(1),iz(1))
 EQUIVALENCE      (ljj(1) ,xljj)  ,(l16,dumm(6))
DATA     nam   / 4HFRSW  ,4H2    ,2*4HBEGN,4HEND ,4HBGIN /
DATA     zero  / 0.0D+0        /
!DB   DATA     DEBUG , ITER    ,MAX  / .FALSE.  ,0     ,3      /

!DB   IF (.NOT.DEBUG) GO TO 20
!     ITER = ITER + 1
!     IF (ITER .GT. MAX) DEBUG = .FALSE.
!     WRITE  (IO,10) NZVB,ITER
!  10 FORMAT ('  .... IN FRSW2.  NZVB =',I8,',   ITER =',I3)
!  20 CONTINUE
nrow = mcblt(2)
CALL frmltd (mcbsma(1),v1(1),v3(1),vb(1))
IF (mcblt(7) < 0) GO TO 200

!     NASTRAN ORIGINAL METHOD

iblk( 1) = mcblt(1)
iblk( 9) = 1
iblk(10) = 1
CALL REWIND (mcblt)
CALL skprec (mcblt,1)

!     FORWARD SWEEP DIRECTLY ON V3

DO  j = 1,nrow
  iblk(8) = -1
  30 CALL getstr (*70,iblk(1))
  ji   = iblk(5)
  ntms = iblk(6)
  ik   = iblk(4)
  IF (ik /= j) GO TO 40
  ntms = ntms - 1
  xljj = xl(ji)
  ji   = ji + 1
  ik   = ik + 1
  40 IF (ntms == 0) GO TO 60
  v3j  = v3(j)
  IF (v3j == zero) GO TO 60
  DO  ii = 1,ntms
    v3(ik) = v3(ik) + xl(ji)*v3j
    ik   = ik + 1
    ji   = ji + 1
  END DO
  60 CALL endget (iblk(1))
  GO TO 30
  70 v3(j)= v3(j)/xljj
END DO

!     BACKWARD SUBSTITUTION OMIT DIAGONAL

IF (nrow == 1) GO TO 500
j    = nrow
90 iblk(8) = -1
100 CALL getstb (*130,iblk(1))
ntms = iblk(6)
ji   = iblk(5)
ik   = iblk(4)
IF (ik-ntms+1 == j) ntms = ntms - 1
IF (ntms == 0) GO TO 120
sum  = zero
DO  ii = 1,ntms
  sum  = sum + xl(ji)*v3(ik)
  ji   = ji - 1
  ik   = ik - 1
END DO
v3(j)= v3(j) + sum
120 CALL endgtb (iblk(1))
GO TO 100
130 IF (j == 1) GO TO 500
j    = j - 1
GO TO 90

!     NEW METHOD

!     THE MCBLT MATRIX HAS BEEN RE-WRITTEN FORWARD FIRST THAN BACKWARD
!     BY UNPSCR IN FEER3. NO STRING OPERATION HERE

200 IF (nam(3) == nam(5)) nam(3) = nam(6)
IF (l16 /= 0) CALL conmsg (nam,3,0)
mcbltx =-mcblt(7)
IF (MOD(mcblt(4),10) /= 2) GO TO 440
nrec = 0
CALL REWIND (mcbltx)
CALL fwdrec (*400,mcbltx)
nwds = mcblt(5)

!     IZ(1)                                                     GINO
!      / V1   V2    V3          VB (OPEN CORE LENGTH = NZVB)   BUFFERS
!     +-----+-----+-----+-----+-------------------------------+--------
!                         OPEN  CORE

!     FORWARD SWEEP DIRECTLY ON V3

ll2  = 0
base = 1
ifb  = +450
DO  j = 1,nrow
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
! 230 FORMAT (5X,I10,' WORDS READ FROM RECORD',I5,'.   LL2 =',I8)
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
  v3j  = v3(j)
  IF (v3j == zero) GO TO 270
  DO  i = ib,ie
    ii   = ii + 1
    v3(ii) = v3(ii) + vb(i)*v3j
  END DO
  270 v3(j)= v3(j)/vb(ib-1)
END DO

!     BACKWARD SUBSTITUTION OMIT DIAGONAL

IF (nrow == 1) GO TO 500
nrec = 0
ll2  = 0
base = 1
j    = nrow
ifb  = -490
DO  jx = 1,nrow
  IF (base < ll2) GO TO 290
  nrec = nrec + 1
!DB   IF (DEBUG) WRITE (IO,210) NREC,IFB
  CALL READ (*400,*280,mcbltx,vb,nzvb,1,ll)
  CALL mesage (-8,0,nam)
  280 ll2  = ll/nwds
!DB   LL3  = LL2/30
!     LL4  = LL2 - LL3
!     IF (DEBUG) WRITE (IO,230) LL,NREC,LL2
  base = 1
  290 xljj = vb(base)
  ii   = ljj(1)
  jj   = ljj(2)
!DB   IF (DEBUG .AND. (BASE.LT.LL3 .OR. BASE.GT.LL4))
!    1    WRITE (IO,250) J,BASE,II,JJ,IFB
  IF (ii /= j) GO TO 420
  ntms = jj - ii + 1
  ib   = base + 2
  ie   = base + ntms
  base = ie + 1
  IF (ntms <= 1) GO TO 310
  sum  = zero
  DO  i = ib,ie
    ii   = ii + 1
    sum  = sum + vb(i)*v3(ii)
  END DO
  v3(j)= v3(j) + sum
  310 j    = j - 1
END DO
GO TO 500

!     ERROR

400 i = mcblt(4)/10
WRITE  (io,410) nrec,j,i,ifb
410 FORMAT ('0*** TRY TO READ RECORD',i5,'.  J,MCBLT(4),IFB =',i7,2I5)
CALL mesage (-2,mcbltx,nam)
420 WRITE  (io,430) ifb,ii,j
430 FORMAT ('0*** ERROR.   IFB),II,J =',i5,1H),2I8)
GO TO 460
440 j = MOD(mcblt(4),10)
WRITE  (io,450) j
450 FORMAT ('0*** MCBLT MATRIX IN WRONG FORM.  UNPSCR FLAG =',i3)
460 CALL mesage (-37,0,nam)

500 nam(3) = nam(5)
IF (l16 /= 0) CALL conmsg (nam,3,0)
RETURN
END SUBROUTINE frsw2
