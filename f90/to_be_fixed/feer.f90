SUBROUTINE feer
     
!     DRIVER FOR THE FEER (FAST EIGENVALUE EXTRACTION ROUTINE) METHOD.
!     THIS ROUTINE WAS CALLED FCNTL BEFORE
 
!     GIVEN A REAL SYMETRIC MATRIX, FEER WILL SOLVE FOR THE EIGENVALUES
!     AND EIGENVECTORS AROUND THE CENTER OF INTEREST
 
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
 
!     IFKAA(7) = 101, MATRIX GINO BLOCK FOR THE INPUT STIFFNESS MATRIX K
!     IFMAA(7) = 102, MATRIX GINO BLOCK FOR THE INPUT MASS MATRIX M
!     IFLELM(7)= 201, MATRIX GINO BLOCK FOR THE OUTPUT EIGENVALUES
!     IFLVEC(7)= 202, MATRIX GINO BLOCK FOR THE OUTPUT EIGENVECTORS
!            ? = 203
!     DMPFLE   = 204, EIGENVALUE SUMMARY FILE
!     SR1FLE-SR8FLE = 301-308, SCRATCH FILES REQUIRED INTERNALLY
!     XLMBDA   =  INPUT, CENTER OF RANGE OF INTEREST.
!                 (USER SPECIFIED SHIFT)
!     NEIG     =  NUMBER OF DESIRED EIGENVALUES AROUND THE CENTER
!                 OF INTEREST. (EIGENVALUES SPECIFIED BY USER)
!     NORD     =  PROBLEM SIZE (SET INTERNALLY USING THE DIMENSION OF
!                 THE STIFFNESS MATRIX)
!     MORD     =  ORDER OF THE REDUCED PROBLEM (SET INTERNALLY)
!     NORTHO   =  NO. OF ORTHOGONAL VECTORS IN PRESENT SET (INCLUDE
!                 PREVISOUSLY COMPUTED VECTORS)
!     EPXM     =  ZERO MASS CRITERIA TO DETERMINE RANK
!     EPX      =  ORTHOGONALITY CONVERGENCE CRITERIA
!     IBK      =  BUCKLING OPTION INDICATOR (SET INTERNALLY)
!     CRITF    =  THE USER SPECIFIED (OR DEFAULT) DESIRED THEORETICAL
!                 ACCURACY OF THE EIGENVALUES EXPRESSED AS A PERCENTAGE
!     LAMBDA   =  VALUE OF THE SHIFT ACTUALLY USED (D.P.)
!     CNDFLG   =  TERMINATION INDICATOR
!     ITER     =  NO. OF STARTING POINTS USED
!     IOPTF    =  SPECIFIED SHIFT OPTION INDICATOR, SET INTERNALLY
!     NOCHNG   =  THEORETICAL ERROR PARAMETER
!     IFSET    =  INTERNALLY COMPUTED SHIFT INDICTOR
!     NONUL    =  NO. OF VETOR ITERATIONS
!     MRANK    =  MATRIX RANK OF THE PROBLEM
!     IND,LMBDA,IDAIG = NOT ACTIVEATED
 
!     EIGENVALUES AND EIGENVECTORS WILL BE STORED ON THE ACTUAL SR1FLE
!     AND SR2FLE. THE SELECTION OF ACCURATE EIGENVALUES AND VECTORS WILL
!     PUT THEM ON IFLELM AND IFLVEC IN THE CORRECT SEQUENCE AT THE END
!     OF PROCESSING
 
!     IFLELM        CONTAINS (K+LAMBDA*M) OR KAA
!     IFLVEC        CONTAINS THE LOWER TRIANGLE L OR C
!     SR4FLE        IS USED AS SCRATCH IN SDCOMP
!     SR5FLE        IS USED AS SCRATCH IN SDCOMP
!     SR6FLE        IS USED AS SCRATCH IN SDCOMP
!     SR7FLE        CONTAINS THE VECTORS WHICH ARE USED TO ORTHOGONALIZE
!     SR8FLE        CONTAINS THE CONTITIONED MAA MATRIX
!     IFLRVA = 301
!     IFLRVC = 302
!     MCBLT         LOWER TRAINGULAR MATRIX L CONTROL BLOCK
!     MCBSMA        CONTITIONED MASTRIX M CONTROL BLOCK
!     MCBVEC        ORTHOGONAL VECTOR FILE CONTROL BLOCK
!     MCBRM         TRIAL VECTOR V OR C(INVERSE-TRANSPOSE)*V CONTROL
!                   BLOCK
 
 INTEGER :: sysbuf    ,cndflg   ,sr8fle   ,NAME(3)  ,  &
     dmpfle    ,iz(12)   ,timed    ,sturm    ,  &
     t1        ,t2       ,t3       ,timet    , mcb(7)    ,icr(2)   ,jcr(2)
 DOUBLE PRECISION :: lambda    ,lmbda    ,dz(1)    ,drsn     ,  &
     drsm      ,epxm     ,scale    ,dsm
 DIMENSION        tmt(4)    ,tml(4)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg / ufm       ,uwm      ,uim
 COMMON  /BLANK / iprob(2)  ,nummod   ,icase
 COMMON  /feercx/ ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7),  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle   ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle   ,  &
     dmpfle    ,nord     ,xlmbda   ,neig     ,  &
     mord      ,ibk      ,critf    ,northo   , iflrva    ,iflrvc
 COMMON  /feerxx/ lambda    ,cndflg   ,iter     ,timed    ,  &
     l16       ,ioptf    ,epx      ,nochng   ,  &
     ind       ,lmbda    ,ifset    ,nzero    ,  &
     nonul     ,idiag    ,mrank    ,istart   , nz3
 COMMON  /reigkr/ option(2)
 COMMON  /zzzzzz/ z(1)
 COMMON  /ntime / lntime    ,tcons(15)
 COMMON  /opinv / mcblt(7)  ,mcbsma(7),mcbvec(7),mcbrm(7)
 COMMON  /system/ ksystm(65)
 COMMON  /packx / itp1      ,itp2     ,iip      ,nnp      , incrp
 COMMON  /unpakx/ iprc      ,ii       ,nn       ,incr
 COMMON  /names / rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw
 COMMON  /sturmx/ sturm     ,shftpt   ,KEEP(2)
 EQUIVALENCE      (iz(1),z(1),dz(1))  ,(ksystm( 1),sysbuf),  &
     (ksystm(2),    io)  ,(ksystm(55),iprec ),  &
     (tcons(8) ,tmt(1))  ,(tcons(12) ,tml(4)), (ksystm(40), nbpw)
 DATA     NAME  / 4HFEER,2*2H   /     ,ibegn/ 4HBEGN   /
DATA     iend  / 4HEND         /     ,mode / 4HMODE   /
DATA     i1,i2 , i3,i4,i0      /  1H1,1H2,1H3,1H4,1H  /
DATA     icr   / 4HPASS,4HFAIL /, jcr/4HFREQ,4HBUCK   /

!     SET PRECISION DIGITS TO 12, ALL MACHINES (NEW 1/92)

it  = 12
epx = 10.**(2-it)
dsm = 10.0D0**(-2*it/3)
NAME(3)  = ibegn
CALL conmsg (NAME,3,0)
CALL feerdd

!     INITIALIZE FEERCX
!     DEFINITION OF INTERNAL PARAMETERS

ibk   = 0
IF (iprob(1) /= mode) ibk = 1
ioptf = ibk
timed = 0
timet = 0
CALL sswtch (16,l16)
IF (l16 == 1) WRITE (io,10)
10 FORMAT (//,' *** DIAG16 - ALL TERMS USED ARE DESCRIBED IN ',  &
    'PROGRAMMER MANUAL  P. 4.48-19I THRU K',/)
lambda = -xlmbda
IF (ibk    ==   0) GO TO 40
IF (xlmbda == 0.0) GO TO 30
CALL page2 (3)
WRITE  (io,20) uwm
20 FORMAT (a25,' 2388', /5X,'USER SPECIFIED RANGE NOT USED FOR FEER',  &
    ' BUCKLING. THE ROOTS OF LOWEST MAGNITUDE ARE OBTAINED')
30 lambda = 0.0D+0
40 ifset  = 0
IF (xlmbda == 0. .AND. ibk == 0) ifset = 1
IF (ifset == 1) ioptf = 1
cndflg = 0
nodcmp = 0
CALL rdtrl (ifkaa(1))
CALL rdtrl (ifmaa(1))
ifk   = ifkaa(1)
ifm   = ifmaa(1)
iprc  = iprec
nord  = ifkaa(2)
incr  = 1
incrp = incr
itp1  = iprc
itp2  = iprc
nz    = korsz(z)
ibuf1 = nz    - sysbuf
ibuf2 = ibuf1 - sysbuf
ntot  = iprc*(5*nord+1) + 4*sysbuf - nz
IF (ntot > 0) CALL mesage (-8,ntot,NAME)
CALL klock (istart)
mrank = 0
CALL gopen  (ifm,z(ibuf1),rdrew)
CALL makmcb (mcb,sr8fle,nord,6,iprc)
CALL gopen  (sr8fle,z(ibuf2),wrtrew)
mcb(2) = 0
mcb(6) = 0
IF (iprc == 2) GO TO 90
DO  j = 1,nord
  ii = 0
  CALL unpack (*60,ifm,z(1))
  nt = nn - ii + 1
  epxm = 0.0D+0
  IF (ii <= j .AND. nn >= j) epxm = z(j-ii+1)*dsm
  ntz = 0
  DO  jj = 1,nt
    IF (ABS(z(jj)) > epxm) CYCLE
    z(jj) = 0.
    ntz   = ntz + 1
  END DO
  IF (ntz < nt) mrank = mrank + 1
  GO TO 70
  60 ii  = 1
  nn  = 1
  nt  = 1
  z(1)= 0.
  70 iip = ii
  nnp = nn
  CALL pack (z(1),sr8fle,mcb(1))
END DO
GO TO 140
90 DO  j = 1,nord
  ii = 0
  CALL unpack (*110,ifm,dz(1))
  nt = nn - ii + 1
  epxm = 0.0D+0
  IF (ii <= j .AND. nn >= j) epxm = dz(j-ii+1)*dsm
  ntz = 0
  DO  jj = 1,nt
    IF (DABS(dz(jj)) > epxm) CYCLE
    dz(jj) = 0.0D+0
    ntz = ntz + 1
  END DO
  IF (ntz < nt) mrank = mrank + 1
  GO TO 120
  110 ii = 1
  nn = 1
  nt = 1
  dz(1) = 0.0D+0
  120 iip = ii
  nnp = nn
  CALL pack (dz(1),sr8fle,mcb(1))
END DO
140 CALL wrttrl (mcb)
mord = 2*(neig-northo) + 10
mrk  = mrank - northo
nzero= northo
IF (mord >   mrk) mord = mrk
IF (neig <= mrank) GO TO 160
CALL page2 (3)
WRITE  (io,150) uwm
150 FORMAT (a25,' 2385', /5X,'DESIRED NUMBER OF EIGENVALUES EXCEED ',  &
    'THE EXISTING NUMBER, ALL EIGENSOLUTIONS WILL BE SOUGHT.')
160 CALL CLOSE (sr8fle,norew)
CALL CLOSE (ifm,rew)
DO  i = 1,7
  mcbsma(i) = mcb(i)
  ifmaa(i)  = mcbsma(i)
END DO
ifm = ifmaa(1)
IF (ibk == 0) GO TO 180

!     SET UP TO DECOMPOSE KAA

iflelm(1) = ifkaa(1)
GO TO 210
180 IF (ifset == 0) GO TO 200

!     CALCULATE INITIAL SHIFT

CALL gopen (ifk,z(ibuf1),rdrew)
CALL gopen (ifm,z(ibuf2),rdrew)
CALL frmax (ifk,ifm,nord,iprc,drsn,drsm)
CALL CLOSE (ifk,rew)
CALL CLOSE (ifm,rew)
scale  = DBLE(FLOAT(nord))*10.0D0**(-it)*drsm
lambda = 10.0D0**(-it/3)*drsn
IF (lambda < scale) lambda = scale

!     CALL IN ADD LINK TO FORM  (K+LAMBDA*M)

200 NAME(2) = i1
CALL conmsg (NAME,3,0)
CALL feer1
NAME(3) = iend
CALL conmsg (NAME,3,0)

!     CALL IN SDCOMP TO DECOMPOSE THIS MATRIX

210 nodcmp  = nodcmp + 1
shftpt  = DABS(lambda)
NAME(2) = i2
NAME(3) = ibegn
CALL conmsg (NAME,3,0)
CALL feer2 (ising)
NAME(3) = iend
CALL conmsg (NAME,3,0)
ik = ibk   + 1
ij = ising + 1
IF (ising /= 1 .AND. l16 == 0) GO TO 230
CALL page2 (4)
WRITE  (io,220) jcr(ik),nord,mrank,mord,northo,neig,nzero,xlmbda,  &
    lambda,icr(ij)
220 FORMAT ('0*** DIAG 16 OUTPUT FOR FEER ANALYSIS, OPTION =',a4, /5X,  &
    'ORDER =',i5,',  MAX RANK =',i5,',  REDUCED ORDER =',i5,  &
    ',  ORTH VCT =',i5,',  NEIG =',i4,',  NZERO =',i4, /5X,  &
    'USER SHIFT =',1P,e16.8,',  INTERNAL SHIFT =',d16.8,  &
    ',  SINGULARITY CHECK ',a4)
230 IF (ising == 0) GO TO 300

!     SINGULAR MATRIX. ADJUST LAMBDA

IF (ibk == 1) GO TO 500
cndflg = cndflg + 1
IF (nodcmp == 3) GO TO 520
lambda = 100.0D0*lambda
GO TO 200

!     DETERMINE THE TIME REQUIRED TO COMPLETE FEER PROCESS

300 CALL tmtogo (t1)
xm  = mord
xmp = northo
xn  = nord
xi  = ifset
ifl = mcblt(1)
CALL gopen (ifl,z(ibuf1),rdrew)
ntms = 0
DO  i = 1,nord
  ii = 0
  CALL unpack (*310,ifl,z(1))
  ntms = ntms + nn - ii + 1
END DO
CALL CLOSE (ifl,rew)
xt = ntms
sp = (xt*(1.-xi)*(xm+xmp)+2.*xm) + xn*(2.+xi)*.5*(3.*xm**2+2.*xmp)  &
    + (16.+11.*xi*.5)*xn*xm + 14.*xm**2

!     OBTAIN TRIDIAGONAL REDUCTION

NAME(2) = i3
NAME(3) = ibegn
CALL conmsg (NAME,3,0)
CALL feer3
NAME(3) = iend
CALL conmsg (NAME,3,0)
IF (cndflg /= 3) GO TO 330
CALL page2 (3)
WRITE  (io,320) uwm
320 FORMAT (a25,' 2389', /5X,'PROBLEM SIZE REDUCED - NO MORE TRIAL ',  &
    'VECTORS CAN BE OBTAINED.')
330 IF (mord == 0) GO TO 350
CALL tmtogo (t2)
timet = t3 - t1

!     OBTAIN EIGENVALUES AND EIGENVECTORS

NAME(2) = i4
NAME(3) = ibegn
CALL conmsg (NAME,3,0)
CALL feer4 (it)
NAME(3) = iend
CALL conmsg (NAME,3,0)
CALL tmtogo (t3)
IF (l16 /= 0) WRITE (io,340) t1,t2,t3,sp
340 FORMAT (' FEER COMPLETE,  T1,T2,T3 =',3I9,',  SP = ',1P,e16.8)
IF (cndflg /= 4) GO TO 370
350 WRITE  (io,360) ufm
360 FORMAT (a23,' 2391, PROGRAM LOGIC ERROR IN FEER')
GO TO 540
370 IF (mord+nzero >= neig) GO TO 390
npr = neig - mord - nzero
CALL page2 (3)
WRITE  (io,380) uwm,npr,neig
380 FORMAT (a25,' 2390', /4X,i5,' FEWER ACCURATE EIGENSOLUTIONS THAN',  &
    ' THE',i5,' REQUESTED HAVE BEEN FOUND.')
cndflg = 1
GO TO 420
390 IF (mord+nzero == neig) GO TO 420
npr = mord + nzero - neig
CALL page2 (3)
WRITE  (io,400) uim,npr,neig
400 FORMAT (a29,' 2392', /4X,i5,' MORE ACCURATE EIGENSOLUTIONS THAN ',  &
    'THE',i5,' REQUESTED HAVE BEEN FOUND.')
IF (l16 == 0) WRITE (io,410)
410 FORMAT (5X,'USE DIAG 16 TO DETERMINE ERROR BOUNDS')
420 CALL gopen (dmpfle,z(ibuf1),wrtrew)

!    SET IZ(1) TO 2 (FOR INVPWR) THEN IZ(7) TO 1 (POINTS TO FEER METHOD)

iz(1) = 2
iz(2) = mord + nzero
iz(3) = iter
iz(4) = 0
iz(5) = nodcmp
iz(6) = nonul
iz(7) = 1
iz(8) = cndflg
iz(9) = 0
iz(10)= 0
iz(11)= 0
iz(12)= 0
CALL WRITE (dmpfle,iz,12,1)
CALL CLOSE (dmpfle,rew)
critf = xn*10.0**(-it)
NAME(2) = i0
CALL conmsg (NAME,3,0)
RETURN

500 WRITE  (io,510) ufm
510 FORMAT (a23,' 2436, SINGULAR MATRIX IN FEER BUCKLING SOLUTION.')
GO TO 540
520 WRITE  (io,530) ufm
530 FORMAT (a23,' 2386', /5X,'STIFFNESS MATRIX SINGULARITY CANNOT BE',  &
    ' REMOVED BY SHIFTING.')
540 CALL mesage (-37,0,NAME)
RETURN
END SUBROUTINE feer
