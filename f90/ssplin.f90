SUBROUTINE ssplin(ni,xyi,nd,xyd,kx,ky,kd,kt,dz,g,ncore,isng)
     
 INTEGER, INTENT(IN)                      :: ni
 REAL, INTENT(IN)                         :: xyi(1)
 INTEGER, INTENT(IN)                      :: nd
 REAL, INTENT(IN)                         :: xyd(1)
 INTEGER, INTENT(IN OUT)                  :: kx
 INTEGER, INTENT(IN OUT)                  :: ky
 INTEGER, INTENT(IN)                      :: kd
 INTEGER, INTENT(IN OUT)                  :: kt
 REAL, INTENT(IN)                         :: dz
 REAL, INTENT(OUT)                        :: g(1)
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(OUT)                     :: isng
 LOGICAL :: lx,ly,lone,ikd,ikt
 DIMENSION  NAME(2)
 REAL :: det
 DATA NAME/4HSSPL,4HIN  /
 
 lone = .true.
 lx = .true.
 ly = .true.
 ikt = .false.
 ikd = .false.
 IF(ky < 0.OR.kx < 0) lone = .false.
 IF(ky < 0.OR.kx > 0) lx   = .false.
 IF(ky > 0.OR.kx < 0) ly   = .false.
 n = ni
 IF(lone) n=n+1
 IF(lx  ) n=n+1
 IF(ly  ) n=n+1
 ex = FLOAT(kx)
 ey = FLOAT(ky)
 IF(kt == 1) ikt = .true.
 IF(kd == 1)  ikd = .true.
 nb = nd*(1+kd)
 
!     CORE NEEDED
 
!                A          G         INVERS
 needed =           nb*ni      + 3*n
!                               B         C
 IF(ikt) needed = needed + nb*n  + ni*n
!                                      C          A OR B
 IF(.NOT.ikt) needed = needed + ni*n + MAX0(n*n,nb*n)
 IF(needed > ncore) CALL mesage(-8,0,NAME)
 is = ncore - 3*n -1
 ig = 1
 
!     IF  KT = 1 COMPUTE B THEN A THEN C IN A SPACE
 
!     IF KT = 0 COMPUTE C THEN A THEN B IN A SPACE
 
 nt = 2*ni
 IF(.NOT.ikt) GO TO 65
 GO TO 95
 
!     COMPUTE TO A MATRIX
 
 1 k = ia
 
!     ZERO A
 
 ii = k+1
 ik = ii + n*n
 DO  i = ii,ik
   g(i) = 0.0
 END DO
 ii = 1
 ik = 0
 DO  i = 1,nt,2
   k = k+ik
   jj = i/2
   DO  j = i,nt,2
     k = k+1
     jj = jj +1
     sum = 0.0
     xm = (xyi(i) - xyi(j)) **2
     xp = (xyi(i) + xyi(j)) **2
     ym = (xyi(i+1) - xyi(j+1)) **2
     yp = (xyi(i+1) + xyi(j+1)) **2
     t1 = xm+ym
     t2 = xp+ym
     t3 = xm+yp
     t4 = xp+yp
     IF(t1 /= 0.0) sum = t1 * ALOG(t1)
     IF(t2 /= 0.0.AND.kx /= 0)sum = sum + (t2*ALOG(t2)*ex)
     IF(t3 /= 0.0.AND.ky /= 0) sum = sum + (t3 * ALOG(t3) * ey)
     IF(t4 /= 0.0.AND.ky /= 0.AND.kx /= 0)sum=sum+(t4*ALOG(t4)*ex*ey)
     IF(j == i) GO TO 10
     g(k) = sum
     
!     SYMETRY TERM
     
     kk = k + (n-1)*(jj-ii)
     g(kk) = sum
     CYCLE
     10 g(k) = sum + dz
     kk = k
   END DO
   inr = 0
   IF(.NOT.lone) GO TO 30
   inr = inr +1
   g(k+inr) = 1.0
   g(kk+inr*n) = 1.0
   30 IF(.NOT.lx) GO TO 40
   inr = inr +1
   g(k+inr) = xyi(i)
   g(kk+inr*n) = xyi(i)
   40 IF(.NOT.ly) GO TO 50
   inr = inr +1
   g(k+inr) = xyi(i+1)
   g(kk+inr*n) = xyi(i+1)
   50 ik = ii + inr
   ii = ii +1
 END DO
 
!     CALL INVERS FOR A-1 C  OR A-1 B
 
!     REPLACE CALLS TO INVAER WITH CALLS TO INVERS
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 isng = -1
 CALL invers(n,g(ia+1),n,g(mp),nc,det,isng,g(is))
 IF(isng == 2) GO TO 1000
 IF(.NOT.ikt) GO TO 100
 ic = ia
 k = ic+1
 GO TO 70
 
!     C MATRIX COLUMN STORED
 
 65 ic = nb*ni
 mp = ic+1
 70 DO  i = 1,ni
   DO  j = 1,n
     ic = ic+1
     g(ic) = 0.0
     IF(i == j) g(ic) = 1.0
   END DO
 END DO
 IF(ikt) GO TO 170
 nc = ni
 ia = ic
 GO TO 1
 
!     B MATRIX COLUMN STORED
 
 95 ib = nb*ni
 mp = ib +1
 GO TO 110
 100 ib = ia
 110 nr = 2*nd
 k = ib +1
 DO  j = 1,nr,2
   DO  i = 1,nt,2
     ib = ib+1
     alt1 = 0.0
     alt2 = 0.0
     alt3 = 0.0
     alt4 = 0.0
     xm = xyd(j)- xyi(i)
     xp = xyi(i) + xyd(j)
     ym = xyd(j+1) - xyi(i+1)
     yp = xyi(i+1) + xyd(j+1)
     t1 = xm*xm + ym*ym
     t2 = xp*xp + ym*ym
     t3 = xm*xm + yp*yp
     t4 = xp*xp +yp*yp
     IF(t1 /= 0.0) alt1 = ALOG(t1)
     IF(t2 /= 0.0.AND.kx /= 0)alt2 = ALOG(t2)
     IF(t3 /= 0.0.AND.ky /= 0)alt3 = ALOG(t3)
     IF(t4 /= 0.0.AND.kx /= 0.AND.ky /= 0) alt4 = ALOG(t4)
     g(ib) = t1*alt1 + t2*alt2*ex + t3*alt3*ey + t4*alt4*ex*ey
     IF(.NOT.ikd) CYCLE
     ik = ib + n
     g(ik) = 2.0*( xm*(1.0+alt1) + xp*(1.0+alt2)*ex +  &
         xm*(1.0+alt3)*ey + xp*(1.0+alt4)*ex*ey)
   END DO
   inr = 0
   IF(.NOT.lone) GO TO 130
   inr = inr +1
   g(ib+inr) = 1.0
   IF(ikd) g(ib+inr+n) = 0.0
   130 IF(.NOT.lx) GO TO 140
   inr = inr +1
   g(ib+inr) = xyd(j)
   IF(ikd) g(ib+inr+n) = 1.0
   140 IF(.NOT.ly) GO TO 150
   inr = inr +1
   g(ib+inr) = xyd(j+1)
   IF(ikd) g(ib+inr+n) = 0.0
   150 ib = ib+inr + n*kd
 END DO
 IF(.NOT.ikt) GO TO 180
 ia = ib
 nc = nb
 GO TO 1
 170 CONTINUE
 
!     GMMATS WANTS ROW STORED SO INVERT ROWS AND COLUMNS AND INVERT
!     MULTIPLICATION ORDER
 
 CALL gmmats(g(mp),nb,n,0,g(k),ni,n,1,g(ig))
 GO TO 1000
 180 CONTINUE
 CALL gmmats(g(mp),ni,n,0,g(k),nb,n,1,g(ig))
 1000 RETURN
END SUBROUTINE ssplin
