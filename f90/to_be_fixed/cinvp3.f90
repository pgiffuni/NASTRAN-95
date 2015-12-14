SUBROUTINE cinvp3
     
!     SUBROUTINE CINVP3, THE MAIN LINK OF CINVPR, SOLVES FOR THE
!     EIGENVALUES AND EIGENVECTORS OF (LAMBDA**2*M + LAMBDA*B*K)
 
!     TYPE DECLARATIONS
 
INTEGER :: comflg    ,REAL     ,END      ,filek
INTEGER :: fileb     ,filelm   ,scrfil
INTEGER :: switch    ,cdp      ,filel    ,sr1fil   ,  &
    fileu     ,sr3fil   ,sr4fil   ,sr8fil   ,  &
    sr9fil    ,sysbuf   ,NAME(2)  ,filevc   ,  &
    timed     ,timeit   ,sqr      ,FILE(7)  , sr2fil    ,s11fil   ,t1       ,t2
REAL :: maxmod
DOUBLE PRECISION :: dz(1)     ,lmbda    ,plus1(2)
DOUBLE PRECISION :: lambda    ,aln(2)   ,alnm1(2) ,eta(2)   ,  &
    etanm1(2) ,h2n(2)   ,h2nm1(2) ,lam1     ,  &
    lm1nm1(2) ,lam2(2)  ,lm2nm1(2),con1(2)  ,  &
    con2(2)   ,cn(2)    ,delta(2) ,xyz(2)
COMMON   /system/  sysbuf    ,nout
COMMON   /cinvpx/  filek(7)  ,filem(7) ,fileb(7) ,filelm(7),  &
    filevc(7) ,dudxx    ,scrfil(11)         , noreg     ,eps
COMMON   /cinfbx/  filel(7)  ,fileu(7)
COMMON   /cinvxx/  lambda(2) ,switch   ,comflg   ,lmbda(2) ,  &
    iterto    ,timed    ,nochng   ,rzero    ,  &
    ind       ,ivect    ,ireg     ,REAL     ,  &
    left      ,northo   ,noroot   ,nzero    ,  &
    lam1(2)   ,maxmod   ,nodes    ,noest    , istart    ,ind1     ,iter     ,isym
COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
    rew       ,norew    ,eofnrw   ,rsp      , rdp       ,csp      ,cdp      ,sqr
COMMON   /cdcmpx/  xxyy(20)  ,iofff
COMMON   /zzzzzz/  z(1)
EQUIVALENCE        (filek(2),ncol)     ,(scrfil(1),sr1fil) ,  &
    (scrfil(2),sr2fil)  ,(scrfil(3),sr3fil) , (scrfil(4),sr4fil)  ,  &
    (scrfil(8),sr8fil)  ,(scrfil(9),sr9fil) ,  &
    (dz(1),z(1))        ,(scrfil(11),s11fil)
DATA      NAME  /  4HCINV,4HP3   /
DATA      plus1 /  +1.d0 ,0.d0   /

!     DEFINITION OF LOCAL PARAMETERS

!     ITER     =
!     IRAPID   =
!     IEP2     =
!     NCOUNT   =
!     IEPCNT   =
!     SWITCH   =
!     A        =
!     EP1      =
!     EP2      =
!     EP3      =
!     GAMMA    =
!     II1      =    POINTER TO U(N)
!     II2      =    POINTER TO U(N-1) OR DELTA U(N)
!     JJ1      =    POINTER TO F(N)
!     JJ2      =    POINTER TO DELTA F(N-1)
!     JJ3      =    POINTER TO F(N-1) OR DELTA F(N)
!     JJ4      =
!     JJ5      =
!     KK1      =    POINTER TO V(N)
!     KK2      =    POINTER TO V(N-1)

10 CONTINUE
timeit= 0
nz    = korsz(z)
ncol2 = ncol  + ncol
ncol4 = ncol2 + ncol2

!     INITIALIZE

cn(1)  = 0.0D0
cn(2)  = 0.0D0
xyz(1) = 0.0D0
xyz(2) = 0.0D0
h2n(1) = 0.0D0
h2n(2) = 0.0D0
lam2(1)= 0.0D0
lam2(2)= 0.0D0
lam1(1)= 0.0D0
lam1(2)= 0.0D0
iter   = 0
CALL klock(t1)
irapid = 0
iep2   = 0
20 ncount = 0
iepcnt = 0
IF (switch == 1) GO TO 30
filel(1) = sr3fil
fileu(1) = sr4fil
GO TO 40
30 filel(1) = sr8fil
fileu(1) = sr9fil
40 filel(5) = cdp
filel(3) = filek(3)
fileu(7) = iofff
FILE(4)  = sqr
FILE(5)  = cdp

!     SET CONVERGENCE CRITERIA

a   = .1
CALL sswtch (12,idiag)
ep1 = .001
ep2 = .02
ep3 = .05
gamma = .01

!     INITILIZE POINTERS TO VECTORS

ii1 = 1
ii2 = ii1 + ncol2
jj1 = ii2 + ncol2
jj2 = jj1 + ncol2
jj3 = jj2 + ncol2
jj5 = jj3 + ncol2
kk1 = jj5 + ncol2
kk2 = kk1 + ncol2
kk3 = kk2 + ncol2
kk4 = kk3 + ncol2
kk5 = kk4 + ncol2
kk6 = kk5 + ncol2
ll1 = kk6 + ncol2
ll2 = ll1 + ncol2
END = (ll2 + ncol2)*2
iobuf = nz - sysbuf + 1
ibuf1 = iobuf - sysbuf
ibuf2 = ibuf1 - sysbuf
ibuf3 = ibuf2 - sysbuf
!     IBUF4 = IBUF3 - SYSBUF
!     IBUF5 = IBUF4 - SYSBUF
!     IBUF6 = IBUF5 - SYSBUF
!     IF (END .GE. IBUF6) GO TO 240
!     NZZ = IBUF6 - END
IF (END >= ibuf3) GO TO 600
nzz = ibuf3 - END
!     IFILE = FILEL(1)
!     CALL OPEN (*610,FILEL,Z(IBUF4),0)
!     IFILE = FILEU(1)
!     CALL OPEN (*610,FILEU,Z(IBUF5),0)
!     IFILE = FILEM(1)
!     CALL OPEN (*610,FILEM,Z(IBUF6),0)

!     GENERATE A STARTING VECTOR

!     FORM U0

IF (left == 1) GO TO 500
lam1(1) = 0.0D0
lam1(2) = 0.0D0
GO TO 70
50 lam1(1) = 0.0D0
lam1(2) = 0.0D0
IF (northo == 0) GO TO 70

!     TEST FOR INSUFFICIENT TIME

CALL klock  (icurnt)
CALL tmtogo (iijjkk)
navg = (icurnt-istart)/northo
IF (iijjkk >= 2*navg) GO TO 70
60 comflg = 8
GO TO 490
70 CONTINUE
IF (ivect == 1) GO TO 90
k = IABS (ind)
DO  i = 1,ncol2,2
  dz(i) = (MOD(k,13)+1)*(1+5*i/ncol)
  k = k + 1
  dz(i+1) = 0.d0
  dz(i  ) = 1.d0/dz(i)
END DO
CALL cnorm1 (dz(ii1),ncol)

!     FORM V0 = LAMBDA*U0

GO TO 110

!     USE PREVIOUSLY STORED VECTOR FOR STARTING VECTOR

90 ifile = filevc(1)
CALL gopen  (filevc,z(iobuf),rd)
CALL bckrec (filevc(1))
in1 = 1
IF (comflg /= 1) GO TO 100
in1 = jj5
CALL bckrec (filevc(1))
100 CALL fread  (filevc,dz(in1),ncol4,1)
IF (comflg == 1) GO TO 140
CALL bckrec (filevc(1))
CALL CLOSE  (filevc(1),norew)
ivect = 0
110 CONTINUE
DO  iu = 1,ncol2,2
  j = kk1 + iu - 1
  dz(j  ) = dz(iu)*lambda(1) - dz(iu+1)*lambda(2)
  dz(j+1) = dz(iu)*lambda(2) + dz(iu+1)*lambda(1)
END DO
IF (northo == 0) GO TO 150
CALL ortho (dz(ii1),dz(kk1),dz(kk2),dz(kk3),dz(kk4),dz(kk5),  &
    dz(kk6),nzz,z(iobuf),z(ibuf1),z(ibuf2),z(ibuf3))
IF (fileb(1) /= 0) GO TO 150
DO  iu = 1,ncol2,2
  j = kk1 + iu - 1
  dz(j  ) = dz(iu)*lambda(1) - dz(iu+1)*lambda(2)
  dz(j+1) = dz(iu)*lambda(2) + dz(iu+1)*lambda(1)
END DO
GO TO 150

!     PICK UP LAST ITERATED VECTOR FOR A STARTING VECTOR

140 CALL fread  (filevc,dz,ncol4,1)
CALL skprec (filevc,-2)
CALL CLOSE  (filevc(1),norew)
GO TO 110
150 CONTINUE
CALL cm tim u (dz(ii1),dz(jj1),0,z(iobuf))
IF (fileb(1) == 0) GO TO 160
FILE(1) = fileb(1)
CALL cm tim u (dz(ii1),dz(kk2),FILE,z(ibuf1))
con1(1) = 2.0D0*lambda(1)
con1(2) = 2.0D0*lambda(2)
CALL cdivid (plus1,con2,con1,2)
con2(1) = -con2(1)
con2(2) = -con2(2)
CALL csub (dz(jj1),dz(kk2),dz(jj1),plus1,con2)
160 CONTINUE
CALL cx trn y (dz(ii1),dz(jj1),aln(1))
CALL csqrtx (aln(1),aln(1))

!     COMPUTE THE R.H.S. OF THE SYSTEM OF EQUATIONS

170 FILE(1) = sr2fil
IF (switch == 1) FILE(1) = s11fil
CALL cm tim u (dz(ii1),dz(ll1),FILE(1),z(iobuf))
CALL cm tim u (dz(kk1),dz(ll2),0      ,z(iobuf))
CALL csub (dz(ll1),dz(ll2),dz(ll2),plus1(1),plus1(1))

!     SHIFT POINTERS

ii  = ii1
ii1 = ii2
ii2 = ii
ii  = jj1
jj1 = jj2
jj2 = jj3
jj3 = ii

!     SAVE THE N-1 VECTOR

IF (switch /= 0) GO TO 190
ixx = jj5 + ncol2 - 1
ixz = ii2
DO  i = jj5,ixx
  dz(i) = dz(ixz)
  ixz = ixz + 1
END DO
190 CONTINUE
CALL tmtogo (ixx)
IF (ixx <= 0) GO TO 60

!     SHIFT PARAMETERS

alnm1(1)  = aln(1)
alnm1(2)  = aln(2)
etanm1(1) = xyz(1)
etanm1(2) = xyz(2)
h2nm1(1)  = h2n(1)
h2nm1(2)  = h2n(2)
lm1nm1(1) = lam1(1)
lm1nm1(2) = lam1(2)
lm2nm1(1) = lam2(1)
lm2nm1(2) = lam2(2)

!     CALL CINFBS TO MAKE ONE ITERATION

CALL cinfbs (dz(ll2),dz(ii1),z(iobuf))
iterto = iterto + 1
iter   = iter + 1
iepcnt = iepcnt + 1
CALL cnorm (dz(ii1),cn(1),dz(ii2))

IF (idiag == 0) GO TO 210
kkkk = ii1 + ncol2 - 1
WRITE  (nout,200) iterto,iter,cn,timed,timeit,(dz(kx),kx=ii1,kkkk)
200 FORMAT (15H iterto =        ,i5,10H iter =     ,i5,' CN =    ',  &
    2D15.5,10H timed =    ,i5,10H timeit=       ,i5,  &
    //,20H iterater vector      , //,(10D12.4))
210 CONTINUE

!     COMPUTE V(N)BAR

con1(1) =-cn(1)/(cn(1)**2 + cn(2)**2)
con1(2) = cn(2)/(cn(1)**2 + cn(2)**2)
CALL csub (dz(ii1),dz(ii2),dz(kk1),lambda,con1)
IF (northo == 0) GO TO 220

!     ORTHOGONALIZE CURRENT ITERANT WITH RESPECT TO VECTORS FOUND IN
!     THE CURRENT AND PREVIOUS REGIONS

CALL ortho (dz(ii1),dz(kk1),dz(kk2),dz(kk3),dz(kk4),dz(kk5),  &
    dz(kk6),nzz,z(iobuf),z(ibuf1),z(ibuf2),z(ibuf3))
220 CONTINUE
IF (fileb(1) /= 0) GO TO 230

!     COMPUTE V(N)

CALL csub (dz(ii1),dz(ii2),dz(kk1),lambda,con1(1))
230 CONTINUE

!     BEGIN TESTING CONVERGENCE CRITERIA

!     COMPUTE F(N)

CALL cm tim u (dz(ii1),dz(jj1),0,z(iobuf))
IF (fileb(1) == 0) GO TO 240
FILE(1) = fileb(1)
CALL cm tim u (dz(ii1),dz(kk2),FILE,z(ibuf1))
con1(1) = 2.0D0*lambda(1)
con1(2) = 2.0D0*lambda(2)
CALL cdivid (plus1,con2,con1,2)
con2(1) = -con2(1)
con2(2) = -con2(2)
CALL csub (dz(jj1),dz(kk2),dz(jj1),plus1,con2)
240 CONTINUE

!     COMPUTE ALPHA(N)

CALL cx trn y (dz(ii1),dz(jj1),aln(1))
CALL csqrtx (aln(1),aln(1))

!     COMPUTE DELTA U(N)

con1(1) = aln(1)/(aln(1)**2 + aln(2)**2)
con1(2) =-aln(2)/(aln(1)**2 + aln(2)**2)
con2(1) = alnm1(1)/(alnm1(1)**2 + alnm1(2)**2)
con2(2) =-alnm1(2)/(alnm1(1)**2 + alnm1(2)**2)
CALL csub (dz(ii1),dz(ii2),dz(ii2),con1(1),con2(1))

!     COMPUTE DELTA F(N)

CALL csub (dz(jj1),dz(jj3),dz(jj3),con1(1),con2(1))
con1(1) = cn(1)*aln(1) - cn(2)*aln(2)
con1(2) = cn(2)*aln(1) + cn(1)*aln(2)
lam1(1) = (alnm1(1)*con1(1) + alnm1(2)*con1(2))/(con1(1)**2 + con1(2)**2)
lam1(2) = (alnm1(2)*con1(1) - alnm1(1)*con1(2))/(con1(1)**2 + con1(2)**2)
IF (irapid == 1) GO TO 410
CALL cx trn y (dz(ii2),dz(jj3),eta(1))
CALL csqrtx (eta(1),xyz(1))

IF (idiag == 0) GO TO 260
WRITE  (nout,250) lam1,xyz,aln
250 FORMAT (12H lambda =    ,2D15.5,12H  eta =        ,2D15.5,  &
    12H alpha =     ,2D15.5)
260 CONTINUE
IF (iter == 1) GO TO 170

!     RAPID CONVERGENCE TEST

!     IF (ETA.GE.A*EPS*GAMMA*(1.+LAMBDA/LAM1)

con1(1) = (lambda(1)*lam1(1) + lambda(2)*lam1(2))/(lam1(1)**2 + lam1(2)**2)
con1(2) = (lambda(2)*lam1(1) - lambda(1)*lam1(2))/(lam1(1)**2 + lam1(2)**2)
IF (DSQRT(xyz(1)**2+xyz(2)**2) >= a*eps*gamma*DSQRT(1.+con1(1)**2  &
    + con1(1)**2+con1(2)**2)) GO TO 280
270 irapid = 1
GO TO 170
280 CONTINUE
IF (DSQRT(etanm1(1)**2+etanm1(2)**2) >= 1.e-06) GO TO 290
IF (DSQRT(xyz(1)**2+xyz(2)**2)-1.01*DSQRT(etanm1(1)**2+  &
    etanm1(2)**2)) 290,290,270

!     EPSILON 2 TEST

290 IF (iep2 == 1) GO TO 310
CALL cx trn y (dz(ii2),dz(jj2),con1(1))
con2(1) = con1(1)*lam1(1) - con1(2)*lam1(2)
con1(2) = con1(1)*lam1(2) + con1(2)*lam1(1)
con1(1) = con2(1)
lam2(1) = (con1(1)*eta(1) + con1(2)*eta(2))/(eta(1)**2 +eta(2)**2)
lam2(2) = (con1(2)*eta(1) - con1(1)*eta(2))/(eta(1)**2 +eta(2)**2)
con1(1) = lam2(1) - lm2nm1(1)
con1(2) = lam2(2) - lm2nm1(2)
h2n(1)  = (con1(1)*lambda(1) + con1(2)*lambda(2))/(lambda(1)**2 +  &
    lambda(2)**2)
h2n(2)  = (con1(2)*lambda(1) - con1(1)*lambda(2))/(lambda(1)**2 +  &
    lambda(2)**2)
IF (iter < 4) GO TO 310
IF (ep2 > DSQRT(h2n(1)**2 + h2n(2)**2).AND.  &
    DSQRT(h2n(1)**2+h2n(2)**2) > DSQRT(h2nm1(1)**2+h2nm1(2)**2)) GO TO 300
GO TO 310
300 iep2 = 1
lam2(1) = lm2nm1(1)
lam2(2) = lm2nm1(2)
310 con1(1) = 1. - (lam2(1)*lam1(1) + lam2(2)*lam1(2))/(lam1(1)**2 +  &
    lam1(2)**2)
con1(2) = (lam2(2)*lam1(1) - lam2(1)*lam1(2))/(lam1(1)**2 + lam1(2)**2)
con2(1) = con1(1)*con1(1) - con1(2)*con1(2)
con1(2) = 2.*con1(2)*con1(1)
con1(1) = con2(1)
con1(1) = DMIN1(DSQRT(con1(1)**2+con1(2)**2),10.0D0)
delta(1)= eta(1)/con1(1)
delta(2)= eta(2)/con1(1)

IF (idiag == 0) GO TO 330
WRITE  (nout,320)lam2,h2n,delta
320 FORMAT (12H  lambda =      ,2D15.5,     12H  h2n =     ,2D15.5,  &
    12H delta =      ,2D15.5)
330 CONTINUE

!     VECTOR CONVERGENCE TEST

IF (DSQRT(delta(1)**2+delta(2)**2) <= (a*eps)**2) GO TO 410
IF (iter <= 3) GO TO 170

!     EPSILON 1 TEST

IF (iepcnt >= 100) GO TO 520
IF (iepcnt >=  10) GO TO 340
IF (DSQRT((lam1(1)-lm1nm1(1))**2+(lam1(2)-lm1nm1(2))**2)  &
    /DSQRT((lambda(1) + DABS(lam1(1)))**2+(lambda(2) +DABS(lam1(2))  &
    )**2) >= ep1) GO TO 170
iepcnt = 0
340 CONTINUE

!     SHIFT DECISION

CALL klock (t2)
timeit = t2-t1
IF (idiag == 0) GO TO 360
WRITE  (nout,350) t2,t1,timeit
350 FORMAT (3I15)
360 CONTINUE
k = DLOG(DSQRT(delta(1)**2+delta(2)**2)/(a*eps)**2)/DABS(DLOG(  &
    DSQRT(lam1(1)**2+lam1(2)**2)/DSQRT(lam2(1)**2+lam2(2)**2)))+1.
k = k/2
IF (idiag == 0) GO TO 380
WRITE  (nout,370) k
370 FORMAT (i5)
380 CONTINUE
ir1 = FLOAT(k-3)*FLOAT(timeit)/FLOAT(iter)
IF (timed >= ir1) GO TO 170
lambda(1) = lambda(1) + lam1(1)
lambda(2) = lambda(2) + lam1(2)

!     STORE THE LAST VECTOR BEFORE A SHIFT FOR USE AS A STARTING VECTOR

IF (switch == 1) GO TO 390
in1 = ii1
GO TO 400
390 in1 = jj5
400 ifile = filevc(1)
CALL gopen (ifile,z(iobuf),wrt)
CALL WRITE (ifile,dz(in1),ncol4,1)
ivect  = 1
comflg = 1

!     STORE  THE CURRENT VECTOR ON THE EIGENVECTOR FILE SO IT CAN BE
!     USED AS THE STARTING VECTOR

CALL WRITE (ifile,dz(ii1),ncol4,1)
CALL CLOSE (ifile,eofnrw)
GO TO 490

!     M  RAPID CONVERGENCE MAKE SURE LAMD1 PASSES EP1 TEST

410 CONTINUE
IF (DSQRT((lam1(1)-lm1nm1(1))**2+(lam1(2)-lm1nm1(2))**2)  &
    /DSQRT((lambda(1) + DABS(lam1(1)))**2+(lambda(2) +DABS(lam1(2))  &
    )**2) >= ep1) GO TO 170

!     CONVERGENCE ACHIEVED, NORMALIZE THE VECTOR

!     STORE THE EIGENVECTOR AND EIGENVALUE ON THE OUTPUT FILES

CALL cnorm1 (dz(ii1),ncol)
lam1(1) = lam1(1) + lambda(1)
lam1(2) = lam1(2) + lambda(2)
inu = ii1 + ncol2 - 1
IF (idiag == 0) GO TO 430
WRITE  (nout,420) lam1,(dz(i),i=ii1,inu)
420 FORMAT (1H1, 20H convergence           ,//,' LAMBDA = ',2D15.5,  &
    //,(10D12.4))
430 CONTINUE
ifile = filevc(1)
CALL gopen (ifile,z(iobuf),wrt)
CALL WRITE (ifile,dz(ii1),ncol4,1)
CALL CLOSE (ifile,eofnrw)
ifile = filelm(1)
CALL gopen (ifile,z(iobuf),wrt)
CALL WRITE (ifile,lam1(1),4,1)
CALL CLOSE (ifile,eofnrw)
northo = northo + 1
noroot = noroot + 1
iep2   = 0
irapid = 0
nochng = 0
comflg = 0
IF (switch == 0) GO TO 440
switch = 0
lambda(1) = lmbda(1)
lambda(2) = lmbda(2)
GO TO 450
440 CONTINUE
ivect = 0
IF (iter <= 5) GO TO 460
450 in1   = jj5
ifile = filevc(1)
CALL gopen (ifile,z(iobuf),wrt)
CALL WRITE (ifile,dz(in1),ncol4,1)
CALL CLOSE (ifile,eofnrw)
ivect = 1
460 iter  = 0

!     COMPUTE PSEUDO LEFT VECTOR

CALL cm tim u (dz(ii1),dz(jj3),0,z(iobuf))
IF (fileb(1) == 0) GO TO 470
CALL cmtimu (dz(ii1),dz(jj2),fileb,z(ibuf1))
con1(1) = 2.0D0*lam1(1)
con1(2) = 2.0D0*lam1(2)
con2(1) =-1.0D0
con2(2) = 0.0D0
CALL csub (dz(jj3),dz(jj2),dz(jj3),con1,con2)
470 IF (isym == 1) GO TO 480

!     LEFT = RIGHT FINISH JOB

CALL cx trn y (dz(ii1),dz(jj3),con1)
CALL cdivid (dz(ii1),dz(jj3),con1,ncol2)

!     PUT SCALED VECTOR ON LEFT VECTOR FILE

480 ifile = scrfil(10)
CALL gopen (ifile,z(ibuf1),wrt)
CALL WRITE (ifile,dz(jj3),ncol4,1)
CALL CLOSE (ifile,eofnrw)
left = 1
IF (isym == 0) GO TO 10

! 490 CALL CLOSE (FILEL,1)
!     CALL CLOSE (FILEU,1)
!     CALL CLOSE (FILEM,1)
490 RETURN

!     RETURN TO MAIN DRIVER TO COMPUTE THE LEFT EIGENVECTOR


!     ENTRY POINT UPON RETURNING FROM OBTAINING THE LEFT VECTOR

500 left = 0
IF (nodes  <= noroot) GO TO 550
IF (noroot >= 3*noest) GO TO 540
aaa = DSQRT((lambda(1) - lam1(1))**2 + (lambda(2)-lam1(2))**2)
IF (aaa <= rzero) GO TO 570
IF (ireg ==    0) GO TO 510
IF (ind > 0) THEN
  GO TO   530
END IF
510 IF (nodes <= noroot) GO TO 550
IF (lam1(1)**2+lam1(2)**2 >= maxmod) GO TO 560

!     GET NEW STARTING  POINT

520 comflg = 0
ind    =-ind
GO TO 490

!     GENERATE NEW ARBITRARY STARTING VECTOR

530 ind   =-(ind+1)
ivect = 0
IF (ind == -13) ind = -1
GO TO 50

!     3*NOEST FOUND

540 comflg = 4
GO TO 490

!     ALL ROOTS IN PROBLEM FOUND

!     COMFLG = 5
!     GO TO 176

!     NO. DES. ROOTS FOUND IN REGION OF CONVERGENCE OUTSIDE REGION

550 comflg = 6
GO TO 490

!     ONE OR MORE ROOTS OUTSIDE REGION

560 comflg = 7
GO TO 490

!     FOUND ROOT OUTSIDE REGION OF CURRENT START POINT

570 ind  = IABS(ind)
ireg = 1
IF (eps*rzero/DSQRT((lam1(1)-lambda(1))**2+(lam1(2)-lambda(2))**2)  &
     < ep3) GO TO 20

!     CURRENT SHIFT POINT IS TOO CLOSE TO AN EIGENVALUE

IF (comflg /= 2) GO TO 580
comflg = 9
got o 490
580 lambda(1) = lambda(1) + .02*rzero
lambda(2) = lambda(2) + .02*rzero
comflg = 2
GO TO 490

!     ERROR EXITS

600 j = -8
GO TO 620
! 610 J = -1
620 CALL mesage (j,ifile,NAME)
RETURN
END SUBROUTINE cinvp3
