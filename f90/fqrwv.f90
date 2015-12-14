SUBROUTINE fqrwv (m,e,er,a,b,w,p,q,xm,INT,zb,srfle, mcbc )
!                                                  SR5FLE SR4FLE
 
 
 INTEGER, INTENT(IN)                      :: m
 DOUBLE PRECISION, INTENT(OUT)            :: e(1)
 DOUBLE PRECISION, INTENT(OUT)            :: er(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: a(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: b(2)
 DOUBLE PRECISION, INTENT(OUT)            :: w(1)
 DOUBLE PRECISION, INTENT(OUT)            :: p(1)
 DOUBLE PRECISION, INTENT(OUT)            :: q(1)
 DOUBLE PRECISION, INTENT(OUT)            :: xm(1)
 LOGICAL, INTENT(OUT)                     :: INT(1)
 REAL, INTENT(IN OUT)                     :: zb(1)
 INTEGER, INTENT(IN OUT)                  :: srfle
 INTEGER, INTENT(IN)                      :: mcbc(7)
 
 
 DOUBLE PRECISION :: pprc     ,zerr    ,ss     ,DSIGN
 DOUBLE PRECISION :: prc    ,hov     ,  &
     sqrt2    ,tol     ,bmax   ,tmax   ,scale   ,  &
     delta    ,eps     ,t      ,x      ,y       ,  &
     s        ,e1      ,e2     ,shift  ,c       ,  &
     gg       ,base    ,sum    ,erf    ,z       ,  &
     f        ,ev      ,x1     ,lambda ,dim     ,  &
     dimf     ,ratio   ,sumx   ,epx    ,epx2    , emax
 DIMENSION  mcb(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg / ufm      ,uwm
 COMMON  /machin/ machx
 COMMON  /lhpwx / lhpw(3)  ,iacc
 COMMON  /feerxx/ lambda   ,cndflg  ,iter   ,timed  ,l16
 COMMON  /names / rd       ,rdrew   ,wrt    ,wrtrew ,rew     , norew   ,eofnrw
 COMMON  /system/ ksystm(65)
 COMMON  /packx / itp1     ,itp2    ,iip    ,nnp    ,incrp
 COMMON  /unpakx/ iprc     ,ii      ,nn     ,incr
 EQUIVALENCE      (ksystm(2),io)    ,(ksystm(55),iprec)
 DATA     ilim  , iexp    ,base /    120, 60, 2.d0      /
 
!     IACC =  MACHINE ACCURACY CONTROL (EPSILON)
!     IACC IS USED TO CONTROL NUMBER UNDERFLOW
!     IEXP AND BASE ARE USED TO CONTROL NUMBER OVERFLOW
 
 IF (m == 1) RETURN
 iprc = 2
 CALL makmcb (mcb(1),srfle,m,2,iprc)
 icf  = mcbc(1)
 incr = 1
 incrp= 1
 itp1 = iprc
 itp2 = iprc
 it   = iacc*iprec
 prc  = 10.d0**(-it)
 pprc = 10.d-4
 jerr = 0
 epx  = 10.d0**(2-it)
 epx2 = epx**2
 hov  = base**iexp
 IF ((machx >= 5 .AND. machx <= 11) .OR. machx == 13 .OR.  &
     machx == 21) hov = base**(iexp-10)
 sqrt2= DSQRT(base)
 m1   = m - 1
 DO  i = 1,m
   e(i) = a(i)
 END DO
 tol  = prc/(10.d0*DBLE(FLOAT(m)))
 bmax = 0.d0
 tmax = 0.d0
 w(m+1) = 0.d0
 DO  i = 1,m
   IF (bmax < DABS(b(i))) bmax = DABS(b(i))
   IF (tmax < DABS(a(i))) tmax = DABS(a(i))
 END DO
 IF (tmax < bmax) tmax = bmax
 scale = 1.d0
 DO  i = 1,ilim
   IF (scale*tmax > hov) GO TO 50
   scale = scale*2.d0
 END DO
 50 IF (bmax == 0.d0) GO TO 170
 DO  i = 1,m
   e(i) =  a(i)*scale
   w(i) = (b(i)*scale)**2
 END DO
 delta= tmax*scale*tol
 eps  = delta*delta
 k    = m
 70 l    = k
 IF (l <= 0) GO TO 140
 l1 = l - 1
 DO  i = 1,l
   k1 = k
   k  = k - 1
   IF (w(k1) <= eps) EXIT
 END DO
 90 IF (k1 /= l) GO TO 100
 w(l) = 0.d0
 GO TO 70
 100 t  = e(l) - e(l1)
 x  = w(l)
 y  = .5D0*t
 s  = DSQRT(x)
 IF (DABS(t) > delta) s = (x/y)/(1.d0+DSQRT(1.d0+x/y**2))
 e1 = e(l ) + s
 e2 = e(l1) - s
 IF (k1 /= l1) GO TO 110
 e(l ) = e1
 e(l1) = e2
 w(l1) = 0.d0
 GO TO 70
 110 shift = e1
 IF (DABS(t) < delta .AND. DABS(e2) < DABS(e1)) shift = e2
 s  = 0.d0
 c  = 1.d0
 gg = e(k1) - shift
 GO TO 130
 120 c  = f/t
 s  = x/t
 x  = gg
 gg = c*(e(k1) - shift) - s*x
 e(k) = (x - gg) + e(k1)
 130 IF (DABS(gg) < delta) gg = gg + c*delta*DSIGN(1.d0,gg)
 f  = gg**2/c
 k  = k1
 k1 = k + 1
 x  = w(k1)
 t  = x + f
 w(k) = s*t
 IF (k < l) GO TO 120
 e(k) = gg + shift
 GO TO 70
 140 DO  i = 1,m
   e(i) = e(i)/scale
 END DO
 DO  l = 1,m1
   k = m - l
   DO  i = 1,k
     IF (e(i) > e(i+1)) CYCLE
     x = e(i)
     e(i  ) = e(i+1)
     e(i+1) = x
   END DO
 END DO
 DO  l = 1,m1
   k = m - l
   DO  i = 1,k
     IF (DABS(e(i)) > DABS(e(i+1))) CYCLE
     x = e(i)
     e(i  ) = e(i+1)
     e(i+1) = x
   END DO
 END DO
 170 IF (m == 0) RETURN
 
!     COMPUTE EIGENVECTORS BY INVERSE ITERATION
 
 erf  = b(m+1)
 mvec = m
 f    = scale/hov
 DO  i = 1,m
   a(i) = a(i)*f
   b(i) = b(i)*f
 END DO
 x1   = 0.d0
 dimf = 10.d0**(-it/3)
 DO  nv = 1,mvec
   ij   = nv
   sumx = 0.d0
   irp  = 0
   IF (nv == 1) GO TO 200
   ratio= DABS(e(nv)/e(nv-1) - 1.d0)
   dim  = .02D0*DABS(1.d0-lambda*e(nv))
   IF (ratio < dim .OR. ratio < dimf) GO TO 220
   nrp = 0
   GO TO 225
   200 nrp = 0
   w(i) = 1.d0
   iip = 1
   nnp = m
   GO TO 330
   
!     MULTIPLE EIGENVALUES
   
   220 nrp = nrp + 1
   225 IF (nv /= 2) GO TO 230
   CALL gopen (srfle,zb(1),wrtrew)
   mcb(2) = 0
   mcb(6) = 0
   GO TO 240
   230 CALL gopen (srfle,zb(1),wrt)
   240 iip = 1
   nnp = m
   CALL pack  (w(1),srfle,mcb(1))
   CALL CLOSE (srfle,norew)
   sum = 0.d0
   ss  = 1.0D0
   DO  i = 1,m
     ss  =-ss
     ij  = ij + 1
     p(i)= FLOAT(MOD(ij,3)+1)/(3.0*FLOAT((MOD(ij,13)+1)*(1+5*i/m)))
     p(i)= p(i)*ss
     sum = sum + p(i)**2
   END DO
   sum = 1.d0/DSQRT(sum)
   DO  i = 1,m
     p(i) = p(i)*sum
     q(i) = p(i)
   END DO
   CALL gopen (srfle,zb(1),rdrew)
   j   = 0
   260 sum = 0.d0
   j   = j + 1
   DO  i = 1,m
     sum = sum + w(i)*p(i)
   END DO
   DO  i = 1,m
     q(i) = q(i) - sum*w(i)
   END DO
   IF (j == (nv-1)) GO TO 290
   ii = 1
   nn = m
   CALL unpack (*290,srfle,w(1))
   GO TO 260
   290 CALL CLOSE (srfle,norew)
   sum = 0.d0
   DO  i = 1,m
     sum = sum + q(i)**2
   END DO
   sum = 1.d0/DSQRT(sum)
   DO  i = 1,m
     q(i) = q(i)*sum
     w(i) = q(i)
   END DO
   330 ev = e(nv)*f
   x  = a(1) - ev
   y  = b(2)
   DO  i = 1,m1
     c  = a(i+1) - ev
     s  = b(i+1)
     IF (DABS(x) >= DABS(s)) GO TO 340
     p(i) = s
     q(i) = c
     INT(i) = .true.
     z = -x/s
     x = y + z*c
     IF (i < m1) y = z*b(i+2)
     GO TO 350
     340 IF (DABS(x) < tol) x = tol
     p(i) = x
     q(i) = y
     INT(i) = .false.
     z = -s/x
     x = c + z*y
     y = b(i+2)
     350 xm(i) = z
   END DO
   IF (DABS(x) < tol) x = tol
   niter = 0
   360 niter = niter + 1
   w(m)  = w(m)/x
   emax  = DABS(w(m))
   DO  l = 1,m1
     i = m-l
     y = w(i) - q(i)*w(i+1)
     IF (INT(i)) y = y - b(i+2)*w(i+2)
     w(i) = y/p(i)
     IF (DABS(w(i)) > emax) emax = DABS(w(i))
   END DO
   sum = 0.d0
   DO  i = 1,m
!WKBR W(I) = (W(I)/EMAX)/EPX
     IF ( emax /= 0.0 ) w(i) = (w(i)/emax)/epx
     IF (DABS(w(i)) < epx2) w(i) = epx2
     sum = sum + w(i)**2
   END DO
   s   = DSQRT(sum)
   DO  i = 1,m
     w(i) = w(i)/s
   END DO
   IF (niter >= 4) GO TO 402
   DO  i = 1,m1
     IF (INT(i)) GO TO 390
     w(i+1) = w(i+1) + xm(i)*w(i)
     CYCLE
     390 y = w(i)
     w(i  ) = w(i+1)
     w(i+1) = y + xm(i)*w(i)
   END DO
   GO TO 360
   402 IF (nv == 1) GO TO 410
   
!     MULTIPLE EIGENVALUES AND ORTHOGONALIZATION
   
   irp = irp + 1
   CALL gopen (srfle,zb(1),rdrew)
   DO  i = 1,m
     q(i) = w(i)
   END DO
   sumx = 0.d0
   jrp  = nv - 1
   DO  i = 1,jrp
     ii = 1
     nn = m
     CALL unpack (*408,srfle,p(1))
     sum = 0.d0
     DO  j = 1,m
       sum = sum + p(j)*q(j)
     END DO
     IF (DABS(sum) > sumx) sumx = DABS(sum)
     DO  j = 1,m
       w(j) = w(j) - sum*p(j)
     END DO
   END DO
   408 CALL CLOSE (srfle,norew)
   410 CONTINUE
   
!     LOGIC SETTING SUM (BY G.CHAN/UNISYS  7/92)
   
!     SUM = PRC*PREC COULD PRODUCE UNDERFLOW (IT=16, PRC=10.**-32)
!     SUM = ZERO, COULD CAUSE DIVIDED BY ZERO AFTER 420 FOR NULL VECTOR
!     SO, WE CHOOSE SUM A LITTLE SMALLER THAN PRC
   
!     SUM = PRC*PRC
!     SUM = 0.0D+0
   sum = prc*1.0D-2
   
   DO  i = 1,m
     IF (DABS(w(i)) >= prc) sum = sum + w(i)*w(i)
   END DO
   sum = 1.d0/DSQRT(sum)
   DO  i = 1,m
     w(i) = w(i)*sum
   END DO
   IF (sumx > 0.9D0 .AND. irp < 3) GO TO 330
   IF (l16 /= 0) WRITE (io,435) nv,niter,irp,sumx
   435 FORMAT (10X,18H feer qrw element ,i5,6H iter ,2I3,6H proj ,d16.8)
   IF (jerr > 0) GO TO 450
   zerr = DABS(w(1))
   DO  i = 2,m
     IF (DABS(w(i)) > zerr) zerr = DABS(w(i))
   END DO
   zerr = (DABS(w(m)))/zerr
   IF (zerr > pprc) jerr = nv - 1
   IF (jerr /=    0) WRITE (io,445) uwm,jerr
   445 FORMAT (a25,' 2399', /5X,'ONLY THE FIRST',i5,' EIGENSOLUTIONS ',  &
       'CLOSEST TO THE SHIFT POINT (F1 OR ZERO) PASS THE FEER ',  &
       'ACCURACY TEST FOR EIGENVECTORS.')
   450 CONTINUE
   CALL pack (w(1),icf,mcbc(1))
   er(nv) = DABS(w(m)*erf/e(nv))
 END DO
 RETURN
END SUBROUTINE fqrwv
