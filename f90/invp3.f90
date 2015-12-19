SUBROUTINE invp3 (norm1,sub,mtimsu,xtrnsy)
     
!     SUBROUTINE INVP3, THE MAIN LINK OF INVPWR, SOLVES FOR THE
!     EIGENVALUES AND EIGENVECTORS OF (K-LAMBDA*M)
!     THIS ROUTINE HANDLES BOTH SINGLE AND DOUBLE PRECISION VERSIONS
 
 
 EXTERNAL          norm1     ,sub      ,mtimsu   ,xtrnsy
INTEGER :: filek     ,END      ,sysbuf   ,sr2fil   ,  &
    sr3fil    ,filel    ,filelt   ,NAME(2)  ,  &
    sr7fil    ,comflg   ,timeit   ,timed    ,  &
    sr8fil    ,switch   ,t1       ,t2       ,  &
    option    ,opt2     ,filem    ,filevc   , filelm    ,mcbvc(7) ,dmpfil
INTEGER :: rew       ,eofnrw
REAL :: lammin    ,lammax
DOUBLE PRECISION :: dz(1)     ,aln      ,alnm1    ,cn       ,  &
    dtemp     ,lam1     ,lm1nm1   ,eta      ,  &
    etanm1    ,lam2     ,lm2nm1   ,h2n      ,  &
    h2nm1     ,delta    ,lambda   ,lmbda    , lam1d     ,freq
COMMON   /zzzzzz/ z(1)
COMMON   /unpakx/ itu       ,iiu      ,jju      ,incru
COMMON   /packx / itp1      ,itp2     ,iip      ,jjp      , incrp
COMMON   /invpwx/ filek(7)  ,filem(7) ,sr1fil(7),sr2fil(7),  &
    filelm    ,filevc   ,sr3fil   ,sr4fil   ,  &
    sr5fil    ,sr6fil   ,sr7fil   ,sr8fil   ,  &
    dmpfil    ,lammin   ,lammax   ,noest    ,  &
    ndplus    ,ndmnus   ,eps      ,northo
COMMON   /system/ ksystm(65)
COMMON   /infbsx/ filel(7)  ,filelt(7)
COMMON   /fbsx  / lfile(7)
COMMON   /names / rd        ,rdrew    ,wrt      ,wrtrew   ,  &
    rew       ,norew    ,eofnrw
COMMON   /invpxx/ lambda    ,comflg   ,iterto   ,timed    ,  &
    nopos     ,rzero    ,neg      ,nochng   ,  &
    ind       ,lmbda    ,switch   ,nzero    ,  &
    noneg     ,ivect    ,ireg     ,istart
COMMON   /reigkr/ option
COMMON   /dcompx/ dumx(20)  ,iofff
COMMON   /trdxx / idummy(27),iopen
COMMON   /regean/ idum40(40),ibuck
EQUIVALENCE       (dz(1)     ,z(1)  ) ,(ksystm( 1),sysbuf ),  &
    (ksystm( 2),ioutpt) ,(ksystm( 9),nlpp   ),  &
    (ksystm(12),nlns  ) ,(ksystm(55),iprec  )
DATA      NAME  / 4HINVP, 4H3       / ,opt2   / 4HUINV    /

!     DEFINITION OF LOCAL PARAMETERS

!     ITER     =  NUMBER OF ITERATIONS FROM THE CURRENT SHIFT POINT
!     IRAPID   =  1 = RAPID CONVERGENCE DO ONE MORE ITERATION
!     IEP2     =  1 = EPSILON 2 TEST FAILED
!     A        =  CONVERGENCE SAFETY FACTOR
!     EP1      =  EPSILON FOR DETERMINING IF IT IS POSSIBLE TO SHIFT
!     EP2      =  EPSILON TO DETERMINE IF LAMBDA 2 IS VALID
!     EP3      =  EPSILON TO DETERMINE IF EIGENVALUE IS TOO CLOSE TO SHI
!     GAMMA    =  CLOSE ROOT CRITERION
!     II1      =  POINTER TO U(N)
!     II2      =  POINTER TO U(N-1) OR DELTA U(N)
!     JJ1      =  POINTER TO F(N)
!     JJ2      =  POINTER TO DELTA F(N-1)
!     JJ3      =  POINTER TO F(N-1) OR DELTA F(N)
!     ALN      =  ALPHA(N)
!     ALNM1    =  ALPHA(N-1)
!     CN       =  NORMALIZATION FACTOR FOR LAST EIGENVECTOR

iopen = -10
CALL sswtch (16,l16)
khr   = 0
nz    = korsz(z)
ncol  = filek(2)
ncol2 = ncol*iprec
CALL makmcb (mcbvc,filevc,ncol,2,iprec)
itu   = iprec
iiu   = 1
jju   = ncol
incru = 1
itp1  = iprec
itp2  = iprec
iip   = 1
jjp   = ncol
incrp = 1

!     INITIALIZE

iter  = 0
irapid= 0
iep2  = 0
kep2  = 0
kold  =-1
kount = 0
GO TO 30

10 IF (northo == 0) GO TO 30
CALL klock  (icurnt)
CALL tmtogo (iijjkk)
navg = (icurnt-istart)/northo
IF (iijjkk >= 2*navg) GO TO 30
20 comflg = 8
GO TO 1140

30 iepcnt = 0
IF (switch == 1) GO TO 40
filel(1)  = sr2fil(1)
filelt(1) = sr3fil
GO TO 50
40 filel(1)  = sr7fil
filelt(1) = sr8fil

50 DO  i = 2,7
  lfile(i) = filek(i)
  filel(i) = filek(i)
END DO
lfile(5) = iprec
filel(5) = iprec
lfile(1) = filel(1)
filelt(7)= iofff

!     SET CONVERGENCE CRITERIA

a   = .1
ep1 = .003
ep2 = .00001
ep2 = .02
ep3 = .05
gamma = .01
IF (l16 == 0 .OR. khr /= 0) GO TO 100
CALL page1
nlns = nlns + 10
WRITE  (ioutpt,70)
70 FORMAT (85H0D i a g   1 6   o u t p u t   f r o m    r o u t i n e  &
    i n v p 3   f o l l o w s . ,//)
WRITE  (ioutpt,80) rzero,eps,gamma,a,ep1,ep2,ep3
80 FORMAT (8H0RZERO =,1P,e13.5,4X,5HEPS =,1P,e13.5,4X,7HGAMMA =,  &
    1P,e13.5,4X,3HA =,1P,e13.5, /,8H ep1   =,1P,e13.5,4X,  &
    5HEP2 =,1P,e13.5,4X,7HEP3   =,1P,e13.5)
WRITE  (ioutpt,90)
90 FORMAT (5H0ITER,5H cflg,11X,3HSTP,11X,3HSHP,10X,4HLAM1,10X,4HLAM2,  &
    11X,3HETA,9X,5HDELTA,4X,1HK,11X,3HH2N,9X,5HLAM1D,/1X,126(1H=))

!     INITIALIZE POINTERS TO VECTORS

100 ii1   = 1
ii2   = ii1 + ncol2
jj1   = ii2 + ncol2
jj2   = jj1 + ncol2
jj3   = jj2 + ncol2
jj4   = jj3 + ncol2
jj5   = jj4 + ncol2
END   = jj5 + ncol2
iend  = END
END   = iend + ncol
ibuf1 = nz   - sysbuf
ibuf2 = ibuf1- sysbuf
ibuf3 = ibuf2- sysbuf
iobuf = ibuf3- sysbuf
IF (END >= iobuf) GO TO 1300

!     GET ORTHOGONALITY FLAGS FOR PREVIOUS EIGENVECTORS

IF (iterto /= 0) GO TO 160
IF (northo == 0) GO TO 170
CALL gopen (filevc,z(iobuf),rdrew)
CALL gopen (filem ,z(ibuf1),rdrew)
DO  i = 1,northo
  ix = iend + i - 1
  z(ix) = 1.0
  CALL unpack (*110,filevc,z(ii1))
  GO TO 140
  110 j = ncol2
  IF (iprec == 2) GO TO 130
  120 z(j) = 0.0
  j = j - 1
  IF (j > 0) GO TO 120
  GO TO 140
  130 dz(j) = 0.0D0
  j = j - 1
  IF (j > 0) GO TO 130
  140 CALL mtimsu (z(ii1),z(jj1),z(ibuf1))
  CALL xtrnsy (z(ii1),z(jj1),dtemp)
  IF (dtemp < 0.0D0) z(ix) = -1.0
END DO
CALL CLOSE (filem ,rew)
CALL CLOSE (filevc,rew)
GO TO 170
160 IF (northo == 0) GO TO 170
ifile = dmpfil
CALL gopen (dmpfil,z(iobuf),rdrew)
CALL READ  (*1310,*1320,dmpfil,z(iend),northo,1,idum)
CALL CLOSE (dmpfil,1)
170 ifile = filem(1)
CALL gopen (ifile,z(ibuf3),rdrew)
ifile = filel(1)
CALL gopen (ifile,z(ibuf1),rdrew)
!WKBNB 1/95    FILELT NOT NEEDED FOR SMCOMP OR SDCOMP - ONLY DECOMP
IF (option /= opt2) GO TO 171
ifile = filelt(1)
CALL gopen (ifile,z(ibuf2),rdrew)
171 CONTINUE
!WKBNE 1/95

!     GENERATE A STARTING VECTOR

IF (ivect == 1) GO TO 240
ksave = k
k = IABS(ind)
IF (iprec == 2) GO TO 210
DO  i = 1,ncol
  z(i) = 1.0/FLOAT((MOD(k,13)+1)*(1+5*i/ncol))
  k = k + 1
END DO
GO TO 230
210 DO  i = 1,ncol
  dz(i) = 1.0D0/FLOAT((MOD(k,13)+1)*(1+5*i/ncol))
  k = k + 1
END DO
230 k = ksave
GO TO 310

!      USE PREVIOUSLY STORED VECTOR AS A STARTING VECTOR

240 ifile = filevc
CALL gopen (filevc,z(iobuf),rd)
CALL bckrec (filevc)
in1 = 1
IF (comflg-1 == 0.0) THEN
  GO TO   250
ELSE
  GO TO   260
END IF
250 in1 = jj5
CALL bckrec (filevc)
260 CALL unpack (*270,filevc,z(in1))
GO TO 300
270 j = in1 + ncol2
IF (iprec == 2) GO TO 290
280 j = j - 1
z(j) = 0.0
IF (j > in1) GO TO 280
GO TO 300
290 j = j - 1
dz(j) = 0.0D0
IF (j   >  in1) GO TO 290
300 IF (comflg == 1) GO TO 320
CALL bckrec (filevc)
CALL CLOSE  (filevc,norew)
ivect = 0
310 CONTINUE
intsub = 1
GO TO 490

!     PICK UP THE LAST ITERATED VECTOR FOR A STARTING VECTOR

320 CONTINUE
CALL unpack (*330,filevc,z)
GO TO 360
330 j = ncol2
IF (iprec == 2) GO TO 350
340 z(j) = 0.0
j = j - 1
IF (j > 0) GO TO 340
GO TO 360
350 dz(j) = 0.0D0
j = j - 1
IF (j > 0) GO TO 350
360 CALL bckrec (filevc)
CALL bckrec (filevc)
CALL CLOSE  (filevc,norew)
GO TO 310

!     SHIFT POINTERS TO VECTORS

400 ii  = ii1
ii1 = ii2
ii2 = ii
ii  = jj1
jj1 = jj2
jj2 = jj3
jj3 = ii
IF (l16 == 0 .OR. khr == 0) GO TO 420
IF (nlns >= nlpp) CALL page1
nlns = nlns + 1
WRITE (ioutpt,410) iterto,comflg, lmbda,lambda,lam1,lam2,eta,delta,k,h2n,lam1d
410 FORMAT (2I5,6(1P,d14.5),i5,2(1P,d14.5))
420 khr = 1

!     SAVE N-1 VECTOR

IF (switch /= 0) GO TO 460
ixx = jj5 + ncol2 - 1
ixz = ii2
IF (iprec /= 2) GO TO 440
DO  i = jj5,ixx,2
  z(i  ) = z(ixz  )
  z(i+1) = z(ixz+1)
  ixz = ixz + 2
END DO
GO TO 460
440 DO  i = jj5,ixx
  z(i) = z(ixz)
  ixz  = ixz + 1
END DO
460 CONTINUE

!     SHIFT PARAMETERS

alnm1  = aln
etanm1 = eta
h2nm1  = h2n
lm1nm1 = lam1
lm2nm1 = lam2

!     CALL INVFBS TO MAKE ONE ITERATION

CALL klock (t1)
IF (option /= opt2) GO TO 470
IF (filel(5) == 2) CALL invfbs (z(jj3),z(ii1),z(iobuf))
IF (filel(5) == 1) CALL intfbs (z(jj3),z(ii1),z(iobuf))
GO TO 480
470 CALL fbsinv (z(jj3),z(ii1),z(iobuf))
480 iterto = iterto + 1
iter   = iter   + 1
iepcnt = iepcnt + 1
CALL tmtogo (ijkk)
IF (ijkk <= 0) GO TO 20
intsub = 2
490 CONTINUE
IF (northo == 0) GO TO 550

!     NORMALIZE CURRENT ITERANT WITH RESPECT TO VECTORS FOUND IN THE
!     CURRENT AND PREVIOUS SEARCH REGIONS

CALL mtimsu (z(ii1),z(jj1),z(iobuf))
ifile = filevc
CALL gopen (filevc,z(iobuf),rdrew)
DO  i = 1,northo
  CALL unpack (*500,filevc,z(jj4))
  GO TO 530
  500 j = jj4 + ncol2
  IF (iprec == 2) GO TO 520
  510 j = j - 1
  z(j) = 0.0
  IF (j > jj4) GO TO 510
  GO TO 530
  520 j = j - 1
  dz(j) = 0.0D0
  IF (j > jj4) GO TO 520
  530 CALL xtrnsy (z(jj4),z(jj1),dtemp)
  ix = iend + i - 1
  dtemp = -dtemp*z(ix)
  CALL sub (z(jj4),z(ii1),dtemp,-1.0D0)
END DO
CALL CLOSE (filevc,norew)
550 CALL norm1 (z(ii1),cn)

!     BEGIN TESTING CONVERGENCE CRITERIA

!     COMPUTE F(N)

CALL mtimsu (z(ii1),z(jj1),z(iobuf))

!     COMPUTE ALPHA(N)

CALL xtrnsy (z(ii1),z(jj1),aln)
aln = DSQRT(DABS(aln))

!     COMPUTE DELTA U(N)

SELECT CASE ( intsub )
  CASE (    1)
    GO TO 400
  CASE (    2)
    GO TO 600
END SELECT
600 CALL sub (z(ii1),z(ii2),1.0D0/aln,1.0D0/alnm1)

!     COMPUTE DELTA F(N)

CALL sub (z(jj1),z(jj3),1.0D0/aln,1.0D0/alnm1)
lam1 = alnm1/(cn*aln)
IF (irapid == 1) GO TO 900
CALL xtrnsy (z(ii2),z(jj3),eta)
eta = DSQRT(DABS(eta))

!     RAPID CONVERGENCE TEST

IF (eta >= a*eps*gamma*DABS(1.0D0+lambda/lam1)) GO TO 620
610 irapid = 1
GO TO 400
620 IF (iter   ==     1) GO TO 400
IF (etanm1 >= 1.e-6) GO TO 700
IF (eta - 1.01*etanm1 > 0.0) THEN
  GO TO   610
END IF

!     EPSILON 2 TEST

700 IF (iep2  ==  1) GO TO 720
IF (eta == 0.d0) GO TO 910
CALL xtrnsy (z(ii2),z(jj2),dtemp)
lam2 = lam1*dtemp/eta**2
h2n  = (lam2-lm2nm1)/lambda
!WKBI 3/94 THE FOLLOWING LINE ADDED TO GET AROUND AN APPARENT COMPILER BUG ON
!          ULTRIX
IF ( eta == 0.d0)PRINT *,' invp3,lam1,dtemp,eta=',lam1,dtemp,eta
IF (iter < 4) GO TO 720
IF (ep2 > DABS(h2n) .AND. DABS(h2n) > DABS(h2nm1)) GO TO 710
GO TO 720
710 CONTINUE
iep2 = 1
lam2 = lm2nm1
720 deltm1 = delta
delta  = eta**2/DMIN1((1.0D0-lam2/lam1)**2,10.0D0)

!     VECTOR CONVERGENCE TEST

IF (DSQRT(delta) <= a*eps) GO TO 910
IF (iter <= 3) GO TO 400

!     EPSILON 1 TEST

IF (iepcnt >= 100) GO TO 1270
IF (iepcnt >=  10) GO TO 800
lam1d = DABS(lam1-lm1nm1)/rzero
IF (lam1d >= DBLE(ep1)) GO TO 400
800 CONTINUE

!     SHIFT DECISION

IF (iepcnt > 5 .AND. delta > deltm1) GO TO 850
IF (DABS(lam2/lam1) > 1.) GO TO 820
IF (kep2 < 0) THEN
  GO TO   850
END IF
810 kep2 = -1
GO TO 400
820 kep2 = 0
CALL klock (t2)
timeit = t2 - t1
k= DLOG(DSQRT(DABS(delta))/(a*eps))/DABS(DLOG(DABS(lam2/lam1)))+1.
k= MIN0(k,9999)
IF (k /= kold) GO TO 830
kount = kount + 1
IF (kount >= 6) GO TO 850
GO TO 840
830 kold  = k
kount = 0
840 CONTINUE
850 lambda= lambda + lam1
k     = 0
kold  =-1
kount = 0
iepcnt= 0
IF (l16  ==    0) GO TO 870
IF (nlns >= nlpp) CALL page1
nlns = nlns + 3
WRITE  (ioutpt,860) lambda
860 FORMAT (18H0NEW shift point =,1P,d14.5,/)

!     STORE THE LAST VECTOR BEFORE A SHIFT FOR USE AS A STARTING VECTOR

870 IF (switch == 1) GO TO 880
in1 = ii1
GO TO 890
880 in1 = jj5
890 ifile = filevc
CALL gopen (filevc,z(iobuf),wrt)
CALL pack  (z(in1),filevc,mcbvc)
ivect  = 1
comflg = 1

!     STORE THE CURRENT VECTOR ON THE EIGENVECTOR FILE SO IT CAN BE
!     USED AS A STARTING VECTOR

CALL pack  (z(ii1),filevc,mcbvc)
CALL CLOSE (filevc,eofnrw)
GO TO 1140

!     MAKE EPSILON 1 TEST

900 IF (DABS (lam1-lm1nm1)/rzero >= DBLE(ep1)) GO TO 400

!     CONVERGENCE ACHIEVED, NORMALIZE THE EIGENVECTOR

910 CALL mtimsu (z(ii1),z(jj1),z(iobuf))
CALL xtrnsy (z(ii1),z(jj1),dtemp)
ix = iend + northo
z(ix) = 1.0
IF (dtemp < 0.0D0) z(ix) = -1.0
dtemp = 1.0D0/DSQRT(DABS(dtemp))
j = ii1
klocal = ii1 + ncol2 - 1
IF (iprec /= 2) GO TO 930
j = (j+1)/2
klocal = klocal/2
DO  i = j,klocal
  dz(i) = dz(i)*dtemp
END DO
GO TO 950
930 DO   i = j,klocal
  z(i) = z(i)*dtemp
END DO
950 CONTINUE

!     STORE THE EIGENVECTOR AND EIGENVALUE ON THE OUTPUT FILES

lam1 = lam1 + lambda
IF (l16  == 0) GO TO 1010
IF (nlns >= nlpp) CALL page1
nlns = nlns + 3
freq = (1.0D0/(8.0D0*DATAN(1.0D0)))*DSQRT(DABS(lam1))
WRITE  (ioutpt,1000) lam1,freq
1000 FORMAT (32H0CONVERGENCE achieved AND lam1 =,1P,d14.5,  &
    7X,'FREQ =',1P,d14.5,'HZ',/)
1010 ifile = filevc
CALL gopen (filevc,z(iobuf),wrt)
CALL pack  (z(ii1),filevc,mcbvc)
CALL CLOSE (filevc,eofnrw)
CALL gopen (filelm,z(iobuf),wrt)
CALL WRITE (filelm,lam1,2,1)
CALL CLOSE (filelm,eofnrw)
CALL CLOSE (sr7fil,eofnrw)
CALL CLOSE (filel,rew)
CALL CLOSE (filelt,rew)
CALL CLOSE (filem,rew)
northo = northo + 1
iep2   = 0
irapid = 0
nochng = 0
IF (lam1 < 0) THEN
  GO TO  1020
ELSE
  GO TO  1030
END IF
1020 IF (ibuck /= 3) GO TO 1030
IF (lam1 >= lammin) noneg = noneg + 1
GO TO 1040
1030 IF (lam1 <= lammax) nopos = nopos + 1
1040 IF (nopos >= ndplus .AND. noneg >= ndmnus) GO TO 1230
IF (northo >= ncol-nzero) GO TO 1220
IF (northo >=    3*noest) GO TO 1210
comflg = 0
IF (switch == 0) GO TO 1050
switch = 0
lambda = lmbda
GO TO 1060
1050 CONTINUE
ivect = 0
IF (iter <= 5) GO TO 1070
1060 in1 = jj5
CALL gopen (filevc,z(iobuf),wrt)
CALL pack  (z(in1),filevc,mcbvc)
CALL CLOSE (filevc,eofnrw)
ivect = 1
1070 iter  = 0

!     TEST IF REGION IS EXHAUSTED

IF (neg < 0) THEN
  GO TO  1120
ELSE IF (neg == 0) THEN
  GO TO  1100
ELSE
  GO TO  1110
END IF

!     NO NEGATIVE REGION

1100 IF (lam1 > lammax) GO TO 1240
GO TO 1130

!     ON POSITIVE SIDE

1110 IF (nopos < ndplus .AND. lam1 <= lammax) GO TO 1130

!     SWITCH TO NEGATIVE SIDE

comflg = 3
GO TO 1140

!     ON NEGATIVE SIDE

1120 IF (noneg >= ndmnus .OR. lam1 < lammin) GO TO 1240

!     CONTINUE ON SAME SIDE

1130 IF (lam1 <= lambda+rzero .AND. lam1 >= lambda-rzero) GO TO 1250
IF (ireg /= 0 .AND. ind > 0) GO TO 1200
comflg = 0
ind = -ind
1140 CALL CLOSE (filel,rew)
CALL CLOSE (filelt,rew)
CALL CLOSE (filem,rew)
CALL wrttrl (mcbvc)
IF (l16  ==    0) GO TO 1150
IF (nlns >= nlpp) CALL page1
nlns = nlns + 1
WRITE (ioutpt,410) iterto,comflg,lmbda,lambda, lam1,lam2,eta,delta,k,h2n,lam1d
1150 IF (northo == 0) RETURN

CALL gopen (dmpfil,z(iobuf),wrtrew)
CALL WRITE (dmpfil,z(iend),northo,1)
CALL CLOSE (dmpfil,1)
RETURN

1200 ind = -(ind+1)
ivect = 0
IF (ind == -13) ind = -1
GO TO 1260
1210 comflg = 4
GO TO 1140
1220 comflg = 5
GO TO 1140
1230 comflg = 6
GO TO 1140
1240 comflg = 7
GO TO 1140
1250 ind = IABS(ind)
ireg = 1
xxx = lam1 - lambda
IF (eps*ABS(rzero) >= ep3*ABS(xxx)) GO TO 1270
1260 IF (northo == 0) GO TO 10
CALL gopen (dmpfil,z(iobuf),wrtrew)
CALL WRITE (dmpfil,z(iend),northo,1)
CALL CLOSE (dmpfil,1)
GO TO 10

!     CURRENT SHIFT POINT TOO CLOSE TO THE EIGENVALUE

1270 IF (comflg /= 2) GO TO 1280
comflg = 9
GO TO 1140
1280 CONTINUE
xxx = lam1 - lambda
lambda = lambda + SIGN(.02,xxx)*rzero
comflg = 2
GO TO 1140

!     ERROR EXITS

1300 no = -8
ifile = END - iobuf
GO TO 1330
1310 no = -2
GO TO 1330
1320 no = -3
1330 CALL mesage (no,ifile,NAME(1))

RETURN
END SUBROUTINE invp3
