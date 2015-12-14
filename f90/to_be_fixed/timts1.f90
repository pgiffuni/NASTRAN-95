SUBROUTINE timts1
     
!     TIMTS1 TIME TESTS GINO AND THE PACK ROUTINES
 
 EXTERNAL        andf
 INTEGER :: sysbuf, output, files(2), f1, f2, buf1, buf2,  &
     END, rd(4), wrt(4), bck(4), mcb(7), eol, eor,  &
     bld(16), INT(16), pak(16), unp(16), TYPE, p,  &
     typin1, typou1, typou2, isubr(2), andf, opt1,  &
     opt2, NAME(4), mask( 9), ablk(15), bblk(15), get(16), put(16)
 REAL :: x(1), z(1)
 DOUBLE PRECISION :: zd, xd
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm, uim, sfm
 COMMON /BLANK / n, m, TYPE, opt1, opt2
 COMMON /system/ sysbuf, output
 COMMON /zblpkx/ zd(2), iz
 COMMON /zntpkx/ xd(2), ix, eol, eor
 COMMON /packx / typin1, typou1, i1, j1, incr1
 COMMON /unpakx/ typou2, i2, j2, incr2
 COMMON /zzzzzz/ a(1)
 EQUIVALENCE     (zd(1),z(1)),(xd(1),x(1))
 DATA    files / 301, 302 / , rd   / 1H , 1H , 1H , 4HREAD /
 DATA    i1000 / 1000     / , i1001/ 1001 /  ,  &
     wrt   / 1H , 1H , 4H   w, 4HRITE /  ,  &
     bck   / 4H   b, 4HACKW,   4HARD , 4HREAD /
 DATA    bld   / 1H , 4HBLDP, 4HK( r, 4HSP ) ,  &
     1H , 4HBLDP, 4HK( r, 4HDP ) , 1H , 4HBLDP, 4HK( c, 4HSP ) ,  &
     1H , 4HBLDP, 4HK( c, 4HDP ) /
 DATA    INT   / 1H , 4HINTP, 4HK( r, 4HSP ) ,  &
     1H , 4HINTP, 4HK( r, 4HDP ) , 1H , 4HINTP, 4HK( c, 4HSP ) ,  &
     1H , 4HINTP, 4HK( c, 4HDP ) /
 DATA    pak   / 1H , 4H pac, 4HK( r, 4HSP ) ,  &
     1H , 4H pac, 4HK( r, 4HDP ) , 1H , 4H pac, 4HK( c, 4HSP ) ,  &
     1H , 4H pac, 4HK( c, 4HDP ) /
 DATA    unp   / 1H , 4HUNPA, 4HK( r, 4HSP ) ,  &
     1H , 4HUNPA, 4HK( r, 4HDP ) , 1H , 4HUNPA, 4HK( c, 4HSP ) ,  &
     1H , 4HUNPA, 4HK( c, 4HDP ) /
 DATA    put   / 4H   p, 4HUTST, 4HR( r, 4HSP ) ,  &
     4H   p, 4HUTST, 4HR( r, 4HDP ) , 4H   p, 4HUTST, 4HR( c, 4HSP ) ,  &
     4H   p, 4HUTST, 4HR( c, 4HDP ) /
 DATA    get   / 4H   g, 4HETST, 4HR( r, 4HSP ) ,  &
     4H   g, 4HETST, 4HR( r, 4HDP ) , 4H   g, 4HETST, 4HR( c, 4HSP ) ,  &
     4H   g, 4HETST, 4HR( c, 4HDP ) /
 DATA    nmask / 9 /
 DATA    isubr / 4HTIMT, 4HS1  /
 
!     INITIALIZE
 
 CALL page1
 f1   = files(1)
 f2   = files(2)
 buf1 = korsz(a) - sysbuf
 buf2 = buf1 - sysbuf
END  = n*m
IF (END >= buf1-1) CALL mesage (-8,0,isubr)
DO  i = 1,END
a(i) = i
END DO
n10  = n*10
m10  = m/10
IF (m10 <= 0) m10 = 1
fn = n
fm = m
p  = 4*(TYPE-1) + 1
mask(1) = 1
DO  i = 2,nmask
  mask(i) = 2*mask(i-1)
END DO
WRITE  (output,11) n, m, TYPE, opt1, opt2
11 FORMAT (1H  , 20X, 25HNASTRAN time test c   n =, i4, 5H, m =, i4 ,  &
    8H, TYPE =,i4, 8H, opt1 =,i4, 8H, opt2 =,i4)

!     WRITE TEST

IF (andf(opt2,mask(1)) == 0) GO TO 50
CALL OPEN (*901,f1,a(buf1),1)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL WRITE (f1,a,m,1)
END DO
CALL cputim (t2,t1,1)
CALL CLOSE  (f1,1)
CALL OPEN (*901,f2,a(buf2),1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL WRITE (f2,a,m10,1)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 30 TO iret
NAME(1) = wrt(1)
NAME(2) = wrt(2)
NAME(3) = wrt(3)
NAME(4) = wrt(4)
GO TO 100

!     READ TEST

30 CONTINUE
IF (andf(opt2,mask(2)) == 0) GO TO 40
CALL OPEN (*901,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL READ (*902,*903,f1,a(i1000),m,1,flag)
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,2)
CALL OPEN (*901,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL READ (*902,*903,f2,a(i1000),m10,1,flag)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,2)
ASSIGN 40 TO iret
NAME(1) = rd(1)
NAME(2) = rd(2)
NAME(3) = rd(3)
NAME(4) = rd(4)
GO TO 100

!     BACKWARD READ TEST

40 CONTINUE
IF (andf(opt2,mask(3)) == 0) GO TO 50
CALL OPEN (*901,f1,a(buf1),2)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL bckrec (f1)
  CALL READ (*902,*903,f1,a(i1000),m,1,flag)
  CALL bckrec (f1)
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
CALL OPEN (*901,f2,a(buf2),2)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL bckrec (f2)
  CALL READ (*902,*903,f2,a(i1000),m10,1,flag)
  CALL bckrec (f2)
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 50 TO iret
NAME(1) = bck(1)
NAME(2) = bck(2)
NAME(3) = bck(3)
NAME(4) = bck(4)
GO TO 100

!     BLDPK TEST

50 CONTINUE
IF (andf(opt2,mask(4)) == 0) GO TO 70
CALL OPEN (*901,f1,a(buf1),1)
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
CALL OPEN (*901,f2,a(buf2),1)
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
ASSIGN 60 TO iret
NAME(1) = bld(p)
NAME(2) = bld(p+1)
NAME(3) = bld(p+2)
NAME(4) = bld(p+3)
GO TO 100

!     INTPK TEST

60 CONTINUE
IF (andf(opt2,mask(5)) == 0) GO TO 70
CALL OPEN (*901,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL intpk (*902,f1,0,TYPE,0)
  DO  j = 1,m
    CALL zntpki
    IF (ix  /= j) GO TO 110
    IF (eol == 0) CYCLE
    IF (ix  /= m) GO TO 110
  END DO
  IF (eol == 0) GO TO 110
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
CALL OPEN (*901,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL intpk (*902,f2,0,TYPE,0)
  DO  j = 1,m10
    CALL zntpki
    IF (ix  /=   j) GO TO 110
    IF (eol ==   0) CYCLE
    IF (ix  /= m10) GO TO 110
  END DO
  IF (eol == 0) GO TO 110
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 70 TO iret
NAME(1) = INT(p)
NAME(2) = INT(p+1)
NAME(3) = INT(p+2)
NAME(4) = INT(p+3)
GO TO 100

!     PACK TEST

70 CONTINUE
IF (andf(opt2,mask(6)) == 0) GO TO 90
CALL makmcb (mcb,f1,m,2,TYPE)
typin1 = TYPE
typou1 = TYPE
i1 = 1
j1 = m
incr1 = 1
mx = m*TYPE
DO  i = 1,mx
  a(i+1000) = i
END DO
CALL OPEN (*901,f1,a(buf1),1)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL pack (a(i1001),f1,mcb)
END DO
CALL cputim (t2,t2,1)
CALL wrttrl (mcb)
CALL CLOSE  (f1,1)
CALL makmcb (mcb,f2,m10,2,TYPE)
j1 = m10
CALL OPEN (*901,f2,a(buf2),1)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL pack (a(i1001),f2,mcb)
END DO
CALL cputim (t4,t4,1)
CALL wrttrl (mcb)
CALL CLOSE  (f2,1)
ASSIGN 80 TO iret
NAME(1) = pak(p)
NAME(2) = pak(p+1)
NAME(3) = pak(p+2)
NAME(4) = pak(p+3)
GO TO 100

!     UNPACK TEST

80 CONTINUE
IF (andf(opt2,mask(7)) == 0) GO TO 90
typou2 = TYPE
i2 = 1
j2 = m
incr2 = 1
CALL OPEN (*901,f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  CALL unpack (*902,f1,a(i1001))
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
j2 = m10
CALL OPEN (*901,f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  CALL unpack (*902,f2,a(i1001))
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,2)
ASSIGN 90 TO iret
NAME(1) = unp(p)
NAME(2) = unp(p+1)
NAME(3) = unp(p+2)
NAME(4) = unp(p+3)
GO TO 100
90 CONTINUE

!     PUTSTR TEST


IF (andf(opt2,mask(8)) == 0) GO TO 220
kerr = 1
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
    91 CALL putstr (ablk)
    IF( nbrstr == 0) GO TO 910
    ablk(7) = MIN0(ablk(6),nbrstr)
    ablk(4) = ablk(4) + ablk(7) + 4
    mm = ablk(7)*nwds
    DO  k = 1,mm
      x(1) = a(k)
    END DO
    IF (ablk(7) == nbrstr) GO TO 93
    CALL endput (ablk)
    nbrstr = nbrstr - ablk(7)
    GO TO 91
    93 IF (j == 10) ablk(8) = 1
    CALL endput (ablk)
  END DO
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
m100 = MAX0(m10/10,1)
CALL gopen (f2,a(buf2),1)
kerr = 2
bblk(1) = f2
bblk(2) = TYPE
bblk(3) = 1
CALL cputim (t3,t3,1)
DO  i = 1,n10
  bblk(4) = 0
  bblk(8) = -1
  DO  j = 1,10
    nbrstr = m100
    202 CALL putstr (bblk)
    IF (nbrstr == 0) GO TO 910
    bblk(7) = MIN0(bblk(6),nbrstr)
    bblk(4) = bblk(4) + bblk(7) + 4
    mm = bblk(7)*nwds
    DO  k = 1,mm
      x(1) = a(k)
    END DO
    IF (bblk(7) == nbrstr) GO TO 206
    nbrstr = nbrstr - bblk(7)
    GO TO 202
    206 IF (j == 10) bblk(8) = 1
    CALL endput (bblk)
  END DO
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 210 TO iret
NAME(1) = put(p)
NAME(2) = put(p+1)
NAME(3) = put(p+2)
NAME(4) = put(p+3)
GO TO 100

!     GETSTR TEST

210 IF (andf(opt2,mask(9)) == 0) GO TO 220
CALL gopen (f1,a(buf1),0)
CALL cputim (t1,t1,1)
DO  i = 1,n
  ablk(8) = -1
  211 CALL getstr (*214,ablk)
  mm = ablk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endget (ablk)
  GO TO 211
END DO
CALL cputim (t2,t2,1)
CALL CLOSE  (f1,1)
CALL gopen  (f2,a(buf2),0)
CALL cputim (t3,t3,1)
DO  i = 1,n10
  bblk(8) = -1
  215 CALL getstr (*219,bblk)
  mm = bblk(6)*nwds
  DO  k = 1,mm
    x(1) = a(k)
  END DO
  CALL endget (bblk)
  GO TO 215
END DO
CALL cputim (t4,t4,1)
CALL CLOSE  (f2,1)
ASSIGN 220 TO iret
NAME(1) = get(p)
NAME(2) = get(p+1)
NAME(3) = get(p+2)
NAME(4) = get(p+3)
GO TO 100

220 CONTINUE
RETURN


!     INTERNAL ROUTINE TO WRITE OUTPUT ONTO THE OUTPUT FILE

100 CONTINUE
time1 = t2 - t1
time2 = t4 - t3
tprrec = 1.0E6*(time2 - time1)/(9.0*fn)
tprwrd = (1.0E6*time1 - fn*tprrec)/(fn*fm)

WRITE  (output,111) NAME, time1, NAME, time2, tprwrd, tprrec
111 FORMAT (1H0, 4A4,  &
    '   N     M-WORD RECORDS -- TIME1 = ', e12.5, ' SECONDS'/ 1X , 4A4,  &
    ' 10*N M/10-WORD RECORDS -- TIME2 = ', e12.5, ' SECONDS'/  &
1H0,'IF THE MODEL IS TIME = (N*M)TPRWRD + N*TPRREC, THEN'/ 1H0, 16X,  &
      '     -- TIME PER WORD   (TPRWRD) = ', e12.5, ' MICRO',  &
      'SECONDS  --  DATA FOR USE IN COMMON /NTIME/'/ 1X , 16X,  &
      '     -- TIME PER RECORD (TPRREC) = ', e12.5, ' MICRO', 'SECONDS')
  
  GO TO iret, (30,40,50,60,70,80,90,210,220)
  
!     INTERNAL ROUTINE CALLED FOR AN ABORT IN THE INTPK TEST
  
  110 CONTINUE
  WRITE  (output,121) sfm,INT(p),INT(p+1),INT(p+2),INT(p+3)
  121 FORMAT (a25,' 2197, ABORT CALLED DURING TIME TEST OF ',4A4)
  
!     ABNORMAL RETURNS FROM GINO--ALL FATAL ERRORS
  
  901 CONTINUE
  902 CONTINUE
  903 CALL mesage (-61,0,0)
  910 WRITE  (output,911) kerr
  911 FORMAT (23H0*** timts1 fatal error,i4 )
  GO TO 903
END SUBROUTINE timts1
