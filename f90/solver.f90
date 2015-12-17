SUBROUTINE solver (lower,x,b,in,out,eps,ifl,scr)
     
!    SOLVER PERFORMS THREE OPERATIONS--
!    1. SOLVES FOR B BY FORWARD-BACKWARD SUBSTITUTION
!    2. COMPUTES OUT = IN + B(T)*X
!    3. IF REQUESTED, COMPUTES EPSILON = NORM(OUT)/NORM(IN)
 
 
 INTEGER, INTENT(IN)                      :: lower
 INTEGER, INTENT(IN OUT)                  :: x
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN)                      :: in
 INTEGER, INTENT(IN OUT)                  :: out
 REAL, INTENT(OUT)                        :: eps
 INTEGER, INTENT(IN OUT)                  :: ifl
 INTEGER, INTENT(IN)                      :: scr
 INTEGER :: filel ,fileu ,fileb ,filex ,  &
     prec  ,SIGN  ,filee ,filef ,fileg ,fileh ,t      ,  &
     signc ,signab,precx ,eol   ,eor   ,sysbuf,scrtch ,  scr1  ,NAME(2)
 DOUBLE PRECISION :: ad    ,num   ,denom
 DIMENSION       filel(7)     ,fileu(7)     ,fileb(7)     ,  &
     filex(7)     ,filee(7)     ,filef(7)     , fileg(7)     ,fileh(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm
 COMMON /zzzzzz/ z(1)
 COMMON /fbsx  / filel ,fileu ,fileb ,filex ,nz    ,prec  ,SIGN   , scr1
 COMMON /mpyadx/ filee ,filef ,fileg ,fileh ,nzz   ,t     ,signab ,  &
     signc ,precx ,scrtch
 COMMON /zntpkx/ ad(2) ,i     ,eol   ,eor
 COMMON /system/ ksystm(65)
 EQUIVALENCE     (ksystm(1),sysbuf)  ,(ksystm(55),iprec)  , (ksystm(2),ioutpt)
 
!     INITIALIZE MATRIX CONTROL BLOCKS FOR FORWARD-BACKWARD SOLUTION
 
 nz = korsz(z)
 filel(1) = lower
 CALL rdtrl (filel)
 fileb(1) = b
 CALL rdtrl (fileb)
 CALL makmcb (filex,x,fileb(3),fileb(4),iprec)
 prec = iprec
 SIGN = -1
 
!     SOLVE A*X = -B FOR X WHERE A HAS BEEN FACTORED
 
 scr1 = scr
 CALL fbs (z,z)
 CALL wrttrl (filex)
 
!     INITIALIZE MATRIX CONTROL BLOCKS FOR MPYAD OPERATION
 
 DO  k = 1,7
   filee(k) = fileb(k)
   filef(k) = filex(k)
 END DO
 fileg(1) = in
 CALL rdtrl (fileg)
 CALL makmcb (fileh,out,fileg(3),fileg(4),iprec)
 nzz = nz
 t   = 1
 signab = 1
 signc  = 1
 precx  = iprec
 scrtch = scr
 
!     COMPUTE OUT = IN + B(T)*X
 
 CALL mpyad  (z,z,z)
 CALL wrttrl (fileh)
 
!     IF REQUESTED,COMPUTE EPS = NORM(OUT) / NORM(IN)
 
 IF (ifl == 0) RETURN
 n1 = nz - sysbuf
 n2 = n1 - sysbuf
 CALL gopen (out,z(n1+1),0)
 CALL gopen ( in,z(n2+1),0)
 num   = 0.0D0
 denom = 0.0D0
 ncol  = fileg(2)
 DO  k = 1,ncol
   CALL intpk (*110,out,0,2,0)
   100 CALL zntpki
   num = num + DABS(ad(1))*DABS(ad(1))
   IF (eol == 0) GO TO 100
   110 CALL intpk (*130,in,0,2,0)
   120 CALL zntpki
   denom = denom + DABS(ad(1))*DABS(ad(1))
   IF (eol == 0) GO TO 120
   130 CONTINUE
 END DO
 IF (denom == 0.0D0) GO TO 160
 eps = DSQRT(num/denom)
 GO TO 180
 160 CALL fname (in,NAME)
 WRITE  (ioutpt,170) uwm,NAME
 170 FORMAT (a25,' 2401, ',2A4,' MATRIX IS NULL.  AN ARBITRARY VALUE ',  &
     'OF 1.0 IS THEREFORE ASSIGNED TO', /5X,  &
     'THE RIGID BODY ERROR RATIO (EPSILON SUB E).')
 eps = 1.0
 180 CALL CLOSE (in, 1)
 CALL CLOSE (out,1)
 RETURN
END SUBROUTINE solver
