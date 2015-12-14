SUBROUTINE fnxtvc (v1,v2,v3,v4,v5,zb,ifn)
     
!     FNXTVC OBTAINS THE REDUCED TRIDIAGONAL MATRIX B WHERE FRBK2
!     PERFORMS THE OPERATIONAL INVERSE.   (DOUBLE PREC VERSION)
 
!           T   -
!      B = V  * A  * V
 
!     V1  = SPACE FOR THE PREVIOUS CURRENT TRIAL VECTOR. INITALLY NULL
!     V2  = SPACE FOR THE CURRENT TRIAL VECTOR. INITIALLY A PSEUDO-
!           RANDOM START VECTOR
!     V3,V4,V5 = WORKING SPACES FOR THREE VECTORS
!     IFN = NO. OF TRIAL VECOTRS EXTRACTED. INITIALLY ZERO.
!     SEE FEER FOR DEFINITIONS OF OTHER PARAMETERS. ALSO PROGRAMMER'S
!           MANUAL PP. 4.48-19G THRU I
 
!     REAL*16, MARKED BY 'CQ', WAS TRIED FOR IMPROVED ACCURACY. BUT THE
!     REAL*16  OPERATIONS ON VAX WERE 10 TIMES SLOWER THAN REAL*8
!     (NUMERIC ACCURACY IS VERY IMPORTANT IN THIS SUBROUTINE)
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: v1(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v2(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v3(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v4(1)
 DOUBLE PRECISION, INTENT(IN)             :: v5(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zb(1)
 INTEGER, INTENT(OUT)                     :: ifn
 INTEGER :: sysbuf    ,cndflg   ,sr5fle   ,NAME(5)  , vcdot
 DOUBLE PRECISION :: lmbda    ,lambda   ,b(2)     , zero
!Q    REAL*16            D         ,DB       ,DSQ      ,SD       ,
 DOUBLE PRECISION :: d         ,db       ,dsq      ,sd       ,  &
     aii       ,dbi      ,depx     ,depx2    ,  &
     sdmax     ,dtmp     ,opdepx   ,omdepx
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON   /xmssg /  ufm       ,uwm
 COMMON   /feercx/  ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7),  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle   ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle   ,  &
     dmpfle    ,nord     ,xlmbda   ,neig     ,  &
     mord      ,ibk      ,critf    ,northo   , iflrva    ,iflrvc
 COMMON   /feerxx/  lambda    ,cndflg   ,iter     ,timed    ,  &
     l16       ,ioptf    ,epx      ,errc     ,  &
     ind       ,lmbda    ,ifset    ,nzero    ,  &
     nonul     ,idiag    ,mrank    ,istart
 COMMON   /system/  ksystm(65)
 COMMON   /opinv /  mcblt(7)  ,mcbsma(7),mcbvec(7),mcbrm(7)
 COMMON   /unpakx/  iprc      ,ii       ,nn       ,incr
 COMMON   /packx /  itp1      ,itp2     ,iip      ,nnp      , incrp
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw
 EQUIVALENCE        (ksystm(1),sysbuf)  ,(ksystm(2),io)
DATA      NAME  /  4HFNXT    ,4HVC     ,2*4HBEGN ,4HEND    /
DATA      vcdot ,  zero /     4HVC.    ,0.0D+0             /

!     SR5FLE CONTAINS THE REDUCED TRIDIAGONAL ELEMENTS

!     SR6FLE CONTAINS THE G VECTORS
!     SR7FLE CONTAINS THE ORTHOGONAL  VECTORS
!     SR8FLE CONTAINS THE CONDITIONED MAA MATRIX

IF (mcblt(7) < 0) NAME(2) = vcdot
NAME(3) = NAME(4)
CALL conmsg (NAME,3,0)
iter  = iter + 1
iprc  = 2
incr  = 1
incrp = incr
itp1  = iprc
itp2  = iprc
ifg   = mcbrm(1)
ifv   = mcbvec(1)
depx  = epx
depx2 = depx**2
opdepx= 1.0D+0 + depx
omdepx= 1.0D+0 - depx
!Q    OPDEPX= 1.0Q+0 + DEPX
!Q    OMDEPX= 1.0Q+0 - DEPX
d     = zero
nord1 = nord - 1

!     NORMALIZE START VECTOR

dsq = zero
IF (ioptf == 1) GO TO 20
CALL frmltd (mcbsma(1),v2(1),v3(1),v5(1))
DO  i = 1,nord
  dsq = dsq + v2(i)*v3(i)
END DO
GO TO 40
20 DO  i = 1,nord
  dsq = dsq + v2(i)*v2(i)
END DO
40 dsq = 1.0D+0/DSQRT(dsq)
!Q 40 DSQ = 1.0D+0/QSQRT(DSQ)
DO  i = 1,nord
  v2(i) = v2(i)*dsq
END DO
IF (northo == 0) GO TO 200

!     ORTHOGONALIZE WITH PREVIOUS VECTORS

DO  i = 1,nord
  v3(i) = v2(i)
END DO
70 DO  ix = 1,14
  nonul = nonul + 1
  CALL gopen (ifv,zb(1),rdrew)
  IF (ioptf == 0) CALL frmltd (mcbsma(1),v2(1),v3(1),v5(1))
  sdmax = zero
  DO  iy = 1,northo
    ii = 1
    nn = nord
    sd = zero
    CALL unpack (*90,ifv,v5(1))
    DO  i = 1,nord
      sd = sd + v3(i)*v5(i)
    END DO
    90 IF (DABS(sd) > sdmax) sdmax = DABS(sd)
!Q 90 IF (QABS(SD) .GT. SDMAX) SDMAX = QABS(SD)
    DO  i = 1,nord
      v2(i) = v2(i) - sd*v5(i)
    END DO
  END DO
  CALL CLOSE (ifv,eofnrw)
  dsq = zero
  IF (ioptf == 1) GO TO 130
  CALL frmltd (mcbsma(1),v2(1),v3(1),v5(1))
  DO  i = 1,nord1
    dsq = dsq + v2(i)*v3(i)
  END DO
  GO TO 150
  130 DO  i = 1,nord1
    dsq = dsq + v2(i)*v2(i)
  END DO
  
! 150 IF (DSQ .LT. DEPX2) GO TO 500
  
!     COMMENTS FORM G.CHAN/UNISYS ABOUT DSQ AND DEPX2 ABOVE,   1/92
  
!     DEPX2 IS SQUARE OF EPX. ORIGINALLY SINCE DAY 1, EPX (FOR VAX AND
!     IBM) IS 10.**-14 AND THEREFORE DEPX2 = 10.**-28. (10.**-24 FOR
!     THE 60/64 BIT MACHINES, USING S.P. COMPUTATION)
!     (EPX WAS SET TO 10.**-10 FOR ALL MACHINES, S.P. AND D.P., 1/92)
  
!     NOTICE THAT DSQ IS THE DIFFERENCE OF TWO CLOSE NUMERIC NUMBERS.
!     THE FINAL VAULES OF DSQ AND THE PRODUCT OF V2*V2 OR V2*V3 APPROACH
!     ONE ANOTHER, AND DEFFER ONLY IN SIGN. THEREFORE, THE NUMBER OF
!     DIGITS (MANTISSA) AS WELL AS THE EXPONENT ARE IMPORTANT HERE
!     (PREVIOUSLY, DO LOOPS 120 AND 140 COVERED 1 THRU NORD)
  
!     MOST OF THE 32 BIT MACHINES HOLD 15 DIGIT IN D.P. WORD, AND SAME
!     FOR THE 64 BIT MACHINES USING S.P. WORD. THEREFORE, CHECKING DSQ
!     DOWN TO 10.**-28 (OR 10.**-24) IS BEYOND THE HARDWARE LIMITS.
!     THIS MAY EXPLAIN SOME TIMES THE RIGID BODY MODES (FREQUENCY = 0.0)
!     GO TO NEGATIVE; IN SOME INSTANCES REACHING -1.E+5 RANGE
  
!     NEXT 7 LINES TRY TO SOLVE THE ABOVE DILEMMA
  
  150 d = v3(nord)
  IF (ioptf == 1) d = v2(nord)
  d = v2(nord)*d
  dtmp = dsq
  dsq  = dsq + d
  IF (dsq < depx2) GO TO 500
  dtmp = DABS(d/dtmp)
!Q    DTMP = QABS(D/DTMP)
  IF (dtmp > omdepx .AND. dtmp < opdepx) GO TO 500
  d = zero
  
  dsq = DSQRT(dsq)
!Q    DSQ = QSQRT(DSQ)
  IF (l16 /= 0) WRITE (io,620) ix,sdmax,dsq
  dsq = 1.0D+0/dsq
  DO  i = 1,nord
    v2(i) = v2(i)*dsq
    v3(i) = v2(i)
  END DO
  IF (sdmax < depx) GO TO 200
END DO

GO TO 500
200 IF (ifn /= 0) GO TO 300

!     SWEEP START VECTOR FOR ZERO ROOTS

dsq = zero
IF (ioptf == 1) GO TO 220
CALL frsw2 (v2(1),v4(1),v3(1),v5(1))
CALL frmltd (mcbsma(1),v3(1),v4(1),v5(1))
DO  i = 1,nord
  dsq = dsq + v3(i)*v4(i)
END DO
GO TO 240
220 CALL frbk2 (v2(1),v4(1),v3(1),v5(1))
DO  i = 1,nord
  dsq = dsq + v3(i)*v3(i)
END DO
240 dsq = 1.0D+0/DSQRT(dsq)
!Q240 DSQ = 1.0D+0/QSQRT(DSQ)
DO  i = 1,nord
  v2(i) = v3(i)*dsq
END DO
GO TO 320

!     CALCULATE OFF DIAGONAL TERM OF B

300 d = zero
DO  i = 1,nord
  d = d + v2(i)*v4(i)
END DO

!     COMMENTS FROM G.CHAN/UNISYS 1/92
!     WHAT HAPPENS IF D IS NEGATIVE HERE? NEXT LINE WOULD BE ALWAY TRUE.

IF (d < depx*DABS(aii)) GO TO 500
!Q    IF (D .LT. DEPX*QABS(AII)) GO TO 500
320 CALL gopen (ifg,zb(1),wrt)
iip = 1
nnp = nord
IF (ioptf == 1) GO TO 330
CALL frsw2 (v2(1),v4(1),v3(1),v5(1))
CALL frmltd (mcbsma(1),v3(1),v4(1),v5(1))
CALL pack (v2(1),ifg,mcbrm(1))
GO TO 350
330 CALL frbk2 (v2(1),v4(1),v3(1),v5(1))
CALL pack (v4(1),ifg,mcbrm(1))
DO  i = 1,nord
  v4(i) = v3(i)
END DO
350 CALL CLOSE (ifg,norew)

!     CALCULATE DIAGONAL TERM OF B

aii = zero
DO  i = 1,nord
  aii = aii + v2(i)*v4(i)
END DO
IF (d == zero) GO TO 420
DO  i = 1,nord
  v3(i) = v3(i) - aii*v2(i) - d*v1(i)
END DO
GO TO 440
420 DO  i = 1,nord
  v3(i) = v3(i) - aii*v2(i)
END DO
440 db = zero
IF (ioptf == 1) GO TO 460
CALL frmltd (mcbsma(1),v3(1),v4(1),v5(1))
DO  i = 1,nord
  db = db + v3(i)*v4(i)
END DO
GO TO 480
460 DO  i = 1,nord
  db = db + v3(i)*v3(i)
END DO
480 db = DSQRT(db)
!Q480 DB = QSQRT(DB)
errc = SNGL(db)
b(1) = aii
b(2) = d
CALL WRITE (sr5fle,b(1),4,1)
CALL gopen (ifv,zb(1),wrt)
iip  = 1
nnp  = nord
CALL pack (v2(1),ifv,mcbvec(1))
CALL CLOSE (ifv,norew)
northo = northo + 1
ifn  = northo - nzero
IF (l16 /= 0) WRITE (io,610) ifn,mord,aii,db,d
IF (ifn >= mord) GO TO 630

!     IF NULL VECTOR GENERATED, RETURN TO OBTAIN A NEW SEED VECTOR

IF (db < depx*DABS(aii)) GO TO 630

!     A GOOD VECTOR IN V2. MOVE IT INTO 'PREVIOUS' VECTOR SPACE V1,
!     NORMALIZE V3 AND V2. LOOP BACK FOR MORE VECTORS.

dbi = 1.0D+0/db
DO  i = 1,nord
  v1(i) = v2(i)
  v3(i) = v3(i)*dbi
  v2(i) = v3(i)
END DO
GO TO 70

500 mord = ifn
WRITE (io,600) uwm,mord
GO TO 630

600 FORMAT (a25,' 2387, PROBLEM SIZE REDUCED TO',i5,' DUE TO -', /5X,  &
    'ORTHOGONALITY DRIFT OR NULL TRIAL VECTOR', /5X,  &
    'ALL EXISTING MODES MAY HAVE BEEN OBTAINED.  USE DIAG 16',  &
    ' TO DETERMINE ERROR BOUNDS',/)
610 FORMAT (5X,'TRIDIAGONAL ELEMENTS ROW (IFN)',i5, /5X,'MORD =',i5,  &
    ', AII,DB,D = ',1P,3D16.8)
620 FORMAT (11X,'ORTH ITER (IX)',i5,',  MAX PROJ (SDMAX)',1P,d16.8,  &
    ',  NORMAL FACT (DSQ)',1P,d16.8)

630 NAME(3) = NAME(5)
CALL conmsg (NAME,3,0)
RETURN
END SUBROUTINE fnxtvc
