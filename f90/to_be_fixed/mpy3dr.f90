SUBROUTINE mpy3dr (z)
     
!     SECONDARY DRIVER IF MPY3DR IS CALLED BY MPY3
!     PRIMARY   DRIVER IF CALLED BY OTHERS (COMB2 AND MRED2 GROUP)
 
!     SETS UP OPEN CORE AND DETERMINES SOLUTION METHOD.
 
 
 INTEGER, INTENT(IN OUT)                  :: z(1)
 IMPLICIT INTEGER (a-z)
 EXTERNAL         andf,orf,complf,lshift
 LOGICAL :: e
 INTEGER :: mpy(3),mcb(7,3),NAME(2)
 REAL :: rhoa,rhob,rhoe,tcol,timcon,timem,timem1,timem2, timem3
 DOUBLE PRECISION :: dd,nn,mm,pp,xx
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBI 4/94
 COMMON /logout/  lout
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /mpy3tl/  filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr3,lkore,code,prec,lcore,scr(7),buf1,buf2, buf3,buf4,e
 COMMON /mpy3cp/  itrl,icore,n,ncb,m,nk,d,maxa,dumcp(34)
 COMMON /ntime /  timcon(16)
 COMMON /system/  sysbuf,nout,dum1(22),diag,dum2(32),meth
 COMMON /mpyadx/  mfilea(7),mfileb(7),mfilee(7),mfilec(7),mcore,  &
     mt,signab,signc,mprec,mscr,timem
 EQUIVALENCE      (ac,filea(2)), (ar,filea(3)),  &
     (bc,fileb(2)), (br,fileb(3)), (bf,fileb(4)), (ec,filee(2)),  &
     (er,filee(3)), (ef,filee(4))
 EQUIVALENCE      (mcb(1,1),filea(1))
 DATA    NAME  /  4HMPY3,4HDR          /
 DATA    mpy   /  4HMPY3,4H    ,4H     /
DATA    jbegn ,  jend  /4HBEGN,4HEND  /

!     RETURN IF EITHER A OR B IS PURGED

IF (filea(1) < 0) RETURN
IF (fileb(1) < 0) RETURN

!     TEST FOR MATRIX COMPATABILITY.

mpy(3) = jbegn
CALL conmsg (mpy,3,0)

scr(1) = scr3
IF (code /= 0) GO TO 5
IF (bf == 2  .OR.  bf ==  7) GO TO 901
5 IF (ar /= br .AND. code == 1) GO TO 902
IF (ar /= bc .AND. code /= 1) GO TO 903
IF (filee(1) <= 0) GO TO 15
e = .true.
IF (code /= 0) GO TO 10
IF (ef == 2 .OR. ef == 7) GO TO 905
10 IF (ec /= bc .AND. code == 1) GO TO 909
IF (ec /= ac .AND. code /= 1) GO TO 904
IF (er /= ac .AND. code == 1) GO TO 910
IF (er /= br .AND. code == 2) GO TO 906
GO TO 30

15 e = .false.
DO  i = 1,7
  filee(i) = 0
END DO

!     CORE ALLOCATION.

30 buf1 = lkore - sysbuf
buf2 = buf1  - sysbuf
buf3 = buf2  - sysbuf
buf4 = buf3  - sysbuf
lcore= buf4  - 1
IF (lcore < 0) GO TO 2008

!     IF REQUESTED CALCULATE THE OUTPUT PRECISION

IF (prec >= 1 .AND. prec <= 4) GO TO 46
iprc = 1
ityp = 0
DO  i = 1,3
  IF (mcb(5,i) == 2 .OR. mcb(5,i) == 4) iprc = 2
  IF (mcb(5,i) >= 3) ityp = 2
END DO
prec = ityp + iprc
IF (prec <= 2) filec(5) = prec
46 CONTINUE

!     DETERMINE NK, THE NUMBER OF COLUMNS OF B MATRIX ABLE TO BE HELD
!     IN CORE.

n   = fileb(3)
ncb = fileb(2)
m   = filea(2)
d   = filea(7) + 1
maxa= filea(6)/filea(5)

!     (NCB SHOULD BE USED IN THE ABOVE EQUATION INSTEAD OF N. SEE
!     MPY3IC)

dd = d
nn = ncb
mm = m
pp = 1 + prec
xx = dd*pp*nn*mm/10000.d0
ixx= xx + 0.5D0
nk = (lcore - 2*ncb - ixx - prec*m - prec - (2+prec)*maxa)/ (2+prec*n)

!     SET UP CONSTANTS IN MPYADX COMMON

mscr  = scr2
mcore = lkore
mprec = 0
signab= 1
signc = 1

!     CALCULATE PROPERTIES OF THE MATRICES

rhoa  = (filea(7)+1)/10000.
rhob  = (fileb(7)+1)/10000.
rhoe  = (filee(7)+1)/10000.
aelms = ar*ac*rhoa
belms = br*bc*rhob
eelms = er*ec*rhoe

!     CALCULATE MPY3 TIME ESTIMATE - REMEMBER NO COMPLEX FOR MPY3

CALL sswtch (19,l19)
timem3 = 1.0E+10
IF (prec >= 3) GO TO 100
IF (code == 1) GO TO 100
timem3 = (rhoa + 2./FLOAT(m))*FLOAT(m)*FLOAT(n)*  &
    (FLOAT(m) + FLOAT(n))*timcon(8+prec) +  &
    (FLOAT(n)**2 + FLOAT(m)**2 + rhoa*FLOAT(m)*  &
    FLOAT(n)*(2. + FLOAT(m)))*timcon(5)
timem3 = timem3/1.0E6
!WKBR 4/94 IF (L19 .NE. 0) WRITE (NOUT,50) FILEA(1),AR,AC,AELMS,RHOA,
IF (l19 /= 0) WRITE (lout,50) filea(1),ar,ac,aelms,rhoa,  &
    fileb(1),br,bc,belms,rhob, filee(1),er,ec,eelms,rhoe,  &
    code,lcore,nk,timem3
50 FORMAT (50H0(a mat  rows  cols   terms    dens) (b mat  rows ,  &
    50H cols   terms    dens) (e mat  rows  cols   terms ,  &
    32H   dens) c  core    nk      time / 3(i6,i7,i6,i9,f7.4,1X),i2,i6,i6,f10.1 )

IF (nk >= 3 .OR. code == 2) GO TO 70
DO  i = 1,7
  mfilea(i) = filea(i)
  mfilee(i) = filee(i)
END DO
CALL makmcb (mfileb,scr1,br,2,prec)
mfileb(2) = ac
tcol = FLOAT(belms)*FLOAT(aelms)/FLOAT(ar)/FLOAT(ac)
mfileb(6) = tcol + 1.0
mfileb(7) = tcol/br*1.0E+4
mfilec(1) = -1
mfilec(5) = prec
mt = 1
CALL mpyad (z(1),z(1),z(1))
timem3 = timem3 + timem

!WKBR 4/94 70 WRITE  (NOUT,80) UIM,TIMEM3
70 WRITE  (lout,80) uim,timem3
80 FORMAT (a29,' 6525, TRIPLE MULTIPLY TIME ESTIMATE FOR MPY3 = ',  &
    f10.1,' SECONDS.')

!     CALCULATE MPYAD TIME ESTIMATE FOR (AT*B)*A + E

100 timem1 = 1.0E+10
IF (code == 2) GO TO 200
DO  i = 1,7
  mfilea(i) = filea(i)
  mfileb(i) = fileb(i)
  IF (code == 1) mfilee(i) = filee(i)
  IF (code /= 1) mfilee(i) = 0
END DO
CALL makmcb (mfilec,-1,ac,2,prec)
mt = 1
CALL mpyad (z(1),z(1),z(1))
timem1 = timem
IF (code == 1) GO TO 130

DO  i = 1,7
  mfileb(i) = mfilea(i)
  mfilea(i) = mfilec(i)
  mfilee(i) = filee(i)
END DO
mfilea(1) = scr1
mfilea(2) = bc
tcol = FLOAT(belms)*FLOAT(aelms)/FLOAT(ar)/FLOAT(bc)
mfilea(6) = tcol + 1.0
mfilea(7) = tcol/ac*1.0E+4
mt = 0
CALL mpyad (z(1),z(1),z(1))
timem1 = timem1 + timem

!WKBR 4/94  130 WRITE  (NOUT,140) UIM,TIMEM1
130 WRITE  (lout,140) uim,timem1
140 FORMAT (a29,' 6525, TRIPLE MULTIPLY TIME ESTIMATE FOR MPYAD - ',  &
    '(AT*B)*A + E = ',f10.1,' SECONDS.')

!     CALCULATE MPYAD TIME ESTIMATE FOR AT*(B*A) + E

200 timem2 = 1.0E+10
IF (code == 1) GO TO 290
DO  i = 1,7
  mfilea(i) = fileb(i)
  mfileb(i) = filea(i)
  IF (code == 2) mfilee(i) = filee(i)
  IF (code /= 2) mfilee(i) = 0
END DO
CALL makmcb (mfilec,-1,br,2,prec)
mt = 0
CALL mpyad (z(1),z(1),z(1))
timem2 = timem
IF (code == 2) GO TO 230

DO  i = 1,7
  mfilea(i) = mfileb(i)
  mfileb(i) = mfilec(i)
  mfilee(i) = filee(i)
END DO
mfileb(1) = scr1
mfileb(2) = ac
tcol = FLOAT(belms)*FLOAT(aelms)/FLOAT(ar)/FLOAT(ac)
mfileb(6) = tcol + 1.0
mfileb(7) = tcol/br*1.0E+4
mt = 1
CALL mpyad (z(1),z(1),z(1))
timem2 = timem2 + timem

!WKBR 4/94 230 WRITE  (NOUT,240) UIM,TIMEM2
230 WRITE  (lout,240) uim,timem2
240 FORMAT (a29,' 6525, TRIPLE MULTIPLY TIME ESTIMATE FOR MPYAD - ',  &
    'AT*(B*A) + E = ',f10.1,' SECONDS.')

!     CHOOSE METHOD BASED ON THE BEST TIME ESTIMATE OR USER REQUEST

290 CALL tmtogo (ttg)
IF (FLOAT(ttg) <= 1.2*AMIN1(timem3,timem1,timem2)) GO TO 908
diag  = andf(diag,complf(lshift(1,18)))
kmeth = meth
jmeth = meth
meth  = 0
IF (jmeth < 1 .OR. jmeth > 3) jmeth = 0
IF (jmeth == 1 .AND. code == 2) jmeth = 0
IF (jmeth == 2 .AND. code == 1) jmeth = 0
IF (jmeth == 3 .AND. code == 1) jmeth = 0
IF (jmeth /= 0) THEN
   SELECT CASE ( jmeth )
    CASE (    1)
      GO TO 400
    CASE (    2)
      GO TO 500
    CASE (    3)
      GO TO 300
  END SELECT
END IF
filec(4) = fileb(4)

IF (timem3 < timem1 .AND. timem3 < timem2) GO TO 300
IF (timem1 < timem2) GO TO 400
GO TO 500

!     PERFORM MULTIPLY WITH MPY3

300 IF (nk < 3) GO TO 310
icore = 0
CALL mpy3ic (z(1),z(1),z(1))
GO TO 9999

!     OUT OF CORE PROCESSING FOR MPY3

310 icore = 1
!WKBR 4/94      WRITE  (NOUT,320) UIM
WRITE  (lout,320) uim
320 FORMAT (a29,' 6526,  THE CENTER MATRIX IS TOO LARGE FOR', /5X,  &
    'IN-CORE PROCESSING.  OUT-OF-CORE PROCESSING WILL BE ', 'PERFORMED.')

nk = (lcore - 4*ncb - prec*m - (2+prec)*maxa)/(2+prec*n)
CALL mpy3oc (z(1),z(1),z(1))
filec(4) = fileb(4)
GO TO 9999

!     PERFORM MULTIPLY WITH MPYAD DOING (AT * B) FIRST

400 DO  i = 1,7
  mfilea(i) = filea(i)
  mfileb(i) = fileb(i)
  IF (code == 1) mfilee(i) = filee(i)
  IF (code /= 1) mfilee(i) = 0
END DO
CALL makmcb (mfilec,scr1,ac,2,prec)
IF (code == 1) mfilec(1) = filec(1)
mt = 1
CALL mpyad (z(1),z(1),z(1))
IF (code == 1) GO TO 425
CALL wrttrl (mfilec)

DO  i = 1,7
  mfileb(i) = mfilea(i)
  mfilea(i) = mfilec(i)
  mfilee(i) = filee(i)
END DO
CALL makmcb (mfilec,filec(1),ac,fileb(4),prec)
mt = 0
CALL mpyad (z(1),z(1),z(1))
425 DO  i = 1,7
  filec(i) = mfilec(i)
END DO
GO TO 9999

!     PERFORM MULTIPLY WITH MPYAD DOING (B*A) FIRST

500 DO  i = 1,7
  mfilea(i) = fileb(i)
  mfileb(i) = filea(i)
  IF (code == 2) mfilee(i) = filee(i)
  IF (code /= 2) mfilee(i) = 0
END DO
CALL makmcb (mfilec,scr1,br,2,prec)
IF (code == 2) mfilec(1) = filec(1)
mt = 0
CALL mpyad (z(1),z(1),z(1))
IF (code == 2) GO TO 525
CALL wrttrl (mfilec)

DO  i = 1,7
  mfilea(i) = mfileb(i)
  mfileb(i) = mfilec(i)
  mfilee(i) = filee(i)
END DO
CALL makmcb (mfilec,filec(1),ac,fileb(4),prec)
mt = 1
CALL mpyad (z(1),z(1),z(1))
525 DO  i = 1,7
  filec(i) = mfilec(i)
END DO
GO TO 9999

!    ERROR MESSAGES.

901 WRITE (nout,9001) ufm
GO TO 1001
902 WRITE (nout,9002) ufm
GO TO 1001
903 WRITE (nout,9003) ufm
GO TO 1001
904 WRITE (nout,9004) ufm
GO TO 1001
905 WRITE (nout,9005) ufm
GO TO 1001
906 WRITE (nout,9006) ufm
GO TO 1001
908 WRITE (nout,9008) ufm
GO TO 1001
909 WRITE (nout,9009) ufm
GO TO 1001
910 WRITE (nout,9010) ufm
1001 CALL mesage (-37,0,NAME)
2008 CALL mesage ( -8,0,NAME)
9001 FORMAT (a23,'6551, MATRIX B IN MPY3 IS NOT SQUARE FOR A(T)BA + E',  &
    ' PROBLEM.')
9002 FORMAT (a23,' 6552, NO. OF ROWS OF MATRIX A IN MPY3 IS UNEQUAL TO'  &
    ,      /5X,'NO. OF ROWS OF MATRIX B FOR A(T)B + E PROBLEM.')
9003 FORMAT (a23,' 6553, NO. OF ROWS OF MATRIX A IN MPY3 IS UNEQUAL TO'  &
    /5X,'NO. OF COLUMNS OF MATRIX B FOR A(T)BA + E PROBLEM.')
9004 FORMAT (a23,' 6554, NO. OF COLUMNS OF MATRIX E IN MPY3 IS UNEQUAL'  &
    ,      /5X,'TO NO. OF COLUMNS OF MATRIX A FOR A(T)BA +E PROBLEM.')
9005 FORMAT (a23,' 6555, MATRIX E IN MPY3 IS NOT SQUARE FOR A(T)BA + ',  &
    'E PROBLEM.')
9006 FORMAT (a23,' 6556, NO. OF ROWS OF MATRIX E IN MPY3 IS UNEQUAL TO'  &
    ,      /5X,'NO. OF ROWS OF MATRIX B FOR BA + E PROBLEM.')
9008 FORMAT (a23,' 6558, INSUFFICIENT TIME REMAINING FOR MPY3 ', 'EXECUTION.')
9009 FORMAT (a23,' 6524, NO. OF COLUMNS OF MATRIX E IN MPY3 IS UNEQUAL'  &
    ,      ' TO',/5X,'NO. OF COLUMNS OF MATRIX B FOR A(T)B + E ', 'PROBLEM.')
9010 FORMAT (a23,' 6559, NO. OF ROWS OF MATRIX E IN MPY3 IS UNEQUAL TO'  &
    ,      /5X,'NO. OF COLUMNS OF MATRIX A FOR A(T)B + E PROBLEM.')

!     RETURN

9999 diag = orf(diag,lshift(l19,18))
meth = kmeth
mpy(3) = jend
CALL conmsg (mpy,3,0)
RETURN
END SUBROUTINE mpy3dr
