SUBROUTINE ssg3a (a,lll,b,x,sr1,sr2,itr1,res)
     
!     SSG3A SOLVES AX = B USING A = L*LT
 
!     ON OPTION COMPUTES RESIDUAL VECTOR RES = A*X - B
!     AND EPSI= X(T)*RES/B(T)*X
 
 
 INTEGER, INTENT(IN OUT)                  :: a
 INTEGER, INTENT(IN)                      :: lll
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN)                      :: x
 INTEGER, INTENT(IN OUT)                  :: sr1
 INTEGER, INTENT(IN OUT)                  :: sr2
 INTEGER, INTENT(IN OUT)                  :: itr1
 INTEGER, INTENT(IN)                      :: res
 INTEGER :: fill,     fillt,    filb, filx,     prec, sysbuf,  &
     NAME(2)
 DOUBLE PRECISION :: dcore(1), dnum,     dnom,     dx
 COMMON /BLANK /  n,        ires,     nskip,    iepsi
 COMMON /fbsx  /  fill(7),  fillt(7), filb(7),  filx(7),  &
     nz,       prec,     ISIGN
 COMMON /zzzzzz/  core(1)
 COMMON /system/  ksystm(55)
 COMMON /unpakx/  itb,      i,        j,        incur
 COMMON /zntpkx/  dx(2),    ik,       ieol,     IEOR
 EQUIVALENCE      (core(1),dcore(1)), (ksystm(1),sysbuf), (ksystm(55),iprec)
 DATA    NAME  /  4HSSG3,   4HA       /
 
 fill(1) = lll
 CALL rdtrl (fill)
 IF (fill(1) <= 0) CALL mesage (-1,lll,NAME)
 filb(1) = b
 CALL rdtrl (filb)
 nload = filb(2)
 nlen  = filb(3)
 ISIGN = 1
 prec  = 2
 nz    = korsz(core)
 DO  i = 2,7
   filx(i) = filb(i)
 END DO
 filx(1) = x
 
!     SAVE DISPLACEMENT VECTOR IN DOUBLE PRECISION
 
 filx(5) = 1
 IF (filb(5) > 2) filx(5) = 3
 filx(5) = filx(5) + iprec - 1
 CALL fbs (core,core)
 CALL wrttrl (filx)
 IF (itr1 < 0) GO TO 130
 fill(1) = res
 CALL rdtrl (fill)
 IF (fill(1) <= 0) GO TO 130
 
!     COMPUTE RESIDUAL VECTOR
 
 CALL ssg2b (a,x,b,res,0,2,-2,sr1)
 
!     COMPUTE EPSI
 
 nz = nz - sysbuf
 CALL gopen (x,core(nz+1),0)
 nz = nz - sysbuf
 CALL gopen (res,core(nz+1),0)
 nz = nz - sysbuf
 CALL gopen (b,core(nz+1),0)
 IF (nz < 2*nlen) GO TO 180
 itb = 2
 incur = 1
 i = 1
 j = nlen
 DO  l = 1,nload
   CALL unpack (*80,x,core)
   dnum = 0.0D0
   dnom = 0.0D0
   CALL intpk (*90,res,0,2,0)
   20 IF (ieol == 0) THEN
     GO TO    30
   ELSE
     GO TO    40
   END IF
   30 CALL zntpki
   dnum = dnum + dx(1)*dcore(ik)
   GO TO 20
   40 CALL intpk (*100,b,0,2,0)
   50 IF (ieol == 0) THEN
     GO TO    60
   ELSE
     GO TO    70
   END IF
   60 CALL zntpki
   dnom = dnom + dx(1)*dcore(ik)
   GO TO 50
   70 epsi = dnum/dnom
   GO TO 110
   80 CALL fwdrec (*160,res)
   90 CALL fwdrec (*170,b)
   100 epsi = 0.0
   110 CALL mesage (35,nskip+l-1,epsi)
   IF (ABS(epsi) < 1.0E-3) CYCLE
   iepsi = -1
   CALL mesage (58,1.0E-3,nskip+l-1)
 END DO
 CALL CLOSE (x,1)
 CALL CLOSE (res,1)
 CALL CLOSE (b,1)
 130 RETURN
 
 150 CALL mesage (-1,ipm,NAME)
 160 ipm = res
 GO TO 150
 170 ipm = b
 GO TO 150
 180 CALL mesage (-8,0,NAME)
 RETURN
 
END SUBROUTINE ssg3a
