SUBROUTINE givens
     
!     DRIVER FOR GIVENS-HOUSEHOLDER METHOD
 
 INTEGER :: sysbuf   ,eigr(4)  ,icore(1) ,option   ,FILE    ,  &
NAME(4)  ,END      ,phia     ,t        ,ix(7)   ,  &
    scr1     ,scr2     ,scr3     ,scr4     ,scr5    , scr6     ,scr7
REAL :: lfreq    ,mb(1)
DOUBLE PRECISION :: dcore(1) ,dlmdas   ,dalpha(2),dbeta(2)
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /xmssg /  ufm      ,uwm      ,uim      ,sfm
COMMON /saddx /  nomat    ,llcore   ,mcba(7)  ,itypa    ,alpha(4),  &
    mcbb(7)  ,itypb    ,beta(4)  ,mcbcde(36),mcbx(7)
COMMON /mgivxx/  dlmdas
COMMON /BLANK /  iprob(2) ,nummod   ,icase    ,xlmdas  &
    /givn  /  dum(100) ,n        ,lfreq    ,order    ,d1      ,  &
    hfreq    ,d2       ,nv       ,d3       ,d4      , nfr  &
    /ntime /  lntime   ,tcons(15)  &
    /regean/  im(7)    ,ik(7)    ,iev(7)   ,scr1     ,scr2    ,  &
    scr3     ,scr4     ,scr5     ,lcore    ,rmax    ,  &
    rmin     ,mz       ,nev      ,epsi     ,rminr   ,  &
    NE       ,nit      ,nevm     ,scr6     ,scr7    ,  &
    nfound   ,lamda    ,ibuck    ,nsym /reigkr/  option  &
    /system/  sysbuf   ,nout     ,nogo     ,ksys(51) ,jprec /zzzzzz/  core(1)
COMMON /condas/  pi       ,twopi    ,radeg    ,degra    ,fps
COMMON /packx /  itp1     ,itp2     ,iip      ,jjp      ,incrp
COMMON /unpakx/  itu      ,iiu      ,jju      ,incru
EQUIVALENCE      (beta(1),dbeta(1)) ,(alpha(1),dalpha(1)),  &
    (slmdas  ,dlmdas ) ,(tcons(8),mb(1)   ),  &
    (core(1),dcore(1)) ,(icore(1),core(1) ),  &
    (tcons(4),apc    ) ,(tcons(5),apu     )
DATA     NAME /  4HGIVE   ,4HNS     ,4HBEGI   ,4HNS    /,  &
mgiv /  4HMGIV   /         ,END   /   4HENDS  /,  &
    icr1 ,  icr2     /301, 302 /


CALL conmsg (NAME,4,0)
i   = 0
kaa = icore(  1)
maa = icore(i+2)
phia= icore(i+3)
DO  i = 1,4
  eigr(i) = icore(i+3)
END DO
nnv = nv
nz  = korsz(core(1))
ibuf1 = nz - 3 - sysbuf
ibuf2 = ibuf1  - sysbuf
ix(1) = kaa
CALL rdtrl (ix)
IF (ix(1) > 0) GO TO 70
WRITE  (nout,60) sfm,ix,kaa,maa,phia
60 FORMAT (a25,' FROM GIVENS.  FILE ERROR,  TRAIL =',5I5,2I8, /5X,  &
    'KAA,MAA,PHIA = ',3I5)
CALL errtrc ('GIVENS  ',60)
70 an = ix(2)

!     CHECK THE CORE SIZE REQUIREMENT FOR WILVEC/VALVEC BEFORE GOING
!     BLINDLY INTO EIGENVALUE COMPUTATION AND EVENTUALLY STOP DUE TO
!     INSUFFICIENT CORE IN THOSE ROUTINES.
!     PRESENTLY CDC IS USING D.P. IN GIVENS COMPUTATION.  IF CDC VERSION
!     IS MODIFIED TO USE S.P., 19 IN THE FOLLOWING FORMULA SHOULD CHANGE
!     TO 10. (COMMENT FROM G.CHAN/UNISYS)

n   = (9*jprec+1)*ix(2) + 2*sysbuf - nz
IF (n > 0) GO TO 120
az  = nz - (3*jprec+1)*ix(2) - 2*sysbuf
az  = az/jprec
am  = SQRT(2.0*az)
ak  = an - am
an2 = an**2
amb = mb(jprec)
av  = nv
anv = an*av
av2 = av**2
t1  = amb*an*(3.0*(an2+anv) + av2)
t23 = apu*(10.0*an2 + 5.0*anv)
t2  = apc*( 5.0*an2 + 3.0*anv + av2) + t23
t3  = 0
IF (am < an) t3 = t23+.5*(apc+apu)*ak*(an2-ak*(an+.5+ak/3.)+an)
t   = (t1+t2+t3)*1.0E-6
n   = an
m   = am
WRITE  (nout,100) uim,t,n,m
100 FORMAT (a29,' 2016, GIVENS TIME ESTIMATE IS ',i8,' SECONDS.',  &
    /36X,'PROBLEM SIZE IS',i8,', SPILL WILL OCCUR FOR THIS ',  &
    'CORE AT A PROBLEM SIZE OF',i8,2H .)
IF (t > 2000 .OR. n > 1000) WRITE (nout,110) uim
110 FORMAT (a29,', FEER METHOD WOULD BE MORE EFFICIENT FOR PROBLEM ',  &
    'OF THIS SIZE',/)
CALL tmtogo (i)
IF (i >= t) GO TO 200
ip1  =-50
FILE = t
GO TO 180
120 WRITE  (nout,150) uim,ix(2),ix(2),n
150 FORMAT (a29,' 3008, INSUFFICIENT CORE FOR GIVENS METHOD.', /5X,  &
    'MATRIX SIZE IS',i5,3H by,i5,'.  ADDITIONAL CORE OF',i7,  &
    ' WORDS IS NEEDED.', /5X,'OR SWITCH TO INVPWR OR FEER ', 'METHOD.')
CALL mesage (-37,0,NAME)
180 CALL mesage (ip1,FILE,NAME)

!     CHOLESKI DECOMPOSE  MAA

200 IF (option /= mgiv) GO TO 250
nomat   = 2
mcba(1) = kaa
mcbb(1) = maa
CALL rdtrl (mcba)
CALL rdtrl (mcbb)
mcbx(1) = icr1
mcbx(2) = mcba(2)
mcbx(3) = mcba(3)
mcbx(4) = mcba(4)
mcbx(5) = jprec
mcbx(6) = 0
mcbx(7) = 0
dalpha(1) = 0.0D0
dalpha(2) = 0.0D0
dbeta(1)  = 0.0D0
dbeta(2)  = 0.0D0
IF (jprec == 2) GO TO 210
slmdas = xlmdas
alpha(1) = 1.0
beta(1)  = slmdas
itypa = 1
itypb = 1
GO TO 220
210 dlmdas = xlmdas
dalpha(1) = 1.0D0
dbeta(1)  = dlmdas
itypa = 2
itypb = 2
220 llcore = nz
CALL sadd (core,core)
CALL wrttrl (mcbx)
ifile1 = icr1
ifile2 = maa
GO TO 260
250 ifile1 = maa
ifile2 = kaa
260 CALL factor (ifile1,scr3,-scr4,scr5,scr6,scr7)

!     C  IS ON SCR3

!     CHANGE SIGNS OF THE OFF-DIAGONAL TERMS OF C AS SDCOMP HAS THEM
!     REVERSED.

ip1   = -5
FILE  = scr3
ix(1) = scr3
CALL rdtrl (ix)
ix(5) = jprec
itp1  = ix(5)
itp2  = itp1
itu   = itp1
incrp = 1
incru = 1
ncol  = ix(2)
ix(1) = scr7
ix(2) = 0
ix(6) = 0
ix(7) = 0
CALL gopen (scr3,core(ibuf1+1),0)
CALL gopen (scr7,core(ibuf2+1),1)
DO  l = 1,ncol
  iiu = 1
  jju = ncol
  CALL unpack (*180,scr3,core)
  IF (itu == 2) GO TO 320
  DO  k = 1,ncol
    core(k) = -core(k)
  END DO
  core(l) = -core(l)
  GO TO 350
  320 DO  k = 1,ncol
    dcore(k) = -dcore(k)
  END DO
  dcore(l) = -dcore(l)
  350 iip = iiu
  jjp = jju
  CALL pack (core,scr7,ix)
END DO
CALL CLOSE (scr3,1)
CALL CLOSE (scr7,1)
CALL wrttrl (ix)

!     C IS NOW ON SCR7

!     INVERT  C

CALL invert (scr7,scr5,scr6)

!     C INVERSE IS ON SCR5


!     GET C INVERSE TRANSPOSE ON SCR6

!     CALL TRANP1 (SCR5,SCR6,4,SCR4,SCR3,SCR7,ICR1,0,0,0,0)
!     GINO UNITS    308, 305,   304, 303, 204, 301
!                   ARE THESE UNITS AVAILABEL?    , 306, 307, 309
!                                                  SCR1,SCR2, EMPTY

!     TRANP1 SHOULD BE 60 PERCENT FASTER BY ADDING 3 MORE SCRATCH FILES

CALL tranp1 (scr5,scr6,7,scr4,scr3,scr7,icr1,scr1,scr2, 309,0)

!     COMPUTE  J

CALL ssg2b (ifile2,scr6,0,scr5,0,jprec,1,scr4)
CALL ssg2b (scr6  ,scr5,0,scr4,1,jprec,1,scr3)

!     J IS ON SCR4

!     EXTRACT EIGENVALUES

CALL valvec

!     TRANSFORM

CALL ssg2b (scr6,scr5,0,scr4,0,jprec,1,scr7)

!     MERGE MODES AND FREE BODY MODES

CALL read6 (icr2,scr4,nfr,phia)
icore(1)= nnv
NAME(3) = END
CALL conmsg (NAME,3,0)
RETURN
END SUBROUTINE givens
