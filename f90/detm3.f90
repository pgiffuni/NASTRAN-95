SUBROUTINE detm3 (*,*,*)
     
!     RMAX   = APPROXIMATE MAGNITUDE OF LARGEST EIGENVALUE OF INTEREST
!     RMIN   = LOWEST  NON-ZERO  EIGENVALUE
!     MZ     = NUMBER OF ZERO EIGENVALUES
!     NEV    = NUMBER OF NON-ZERO EIGENVALUES IN RANGE OF INTEREST
!     EPSI   = CONVERGENCE CRITERION
 
!     NEVM   = MAXIMUM NUMBER OF EIGENVALUES DESIRED
!     IS     = STARTING SET COUNTER
!     IC     = COUNTER FOR CHANGE OF CONVERGENCE CRITERIA
!     NFOUND = THE NUMBER OF EIGENVALUES FOUND TO DATA
!     IM     = MASS MATRIX CONTROL BLOCK
!     IK     = K MATRIX CONTROL BLOCK
!     IEV    = EIGENVECTOR CONTROL BLOCK
 
!     A      = M + P*K
 
 INTEGER :: prec,u2,scr1,scr2,scr3,scr4,scr5,scr6,scr7
 DOUBLE PRECISION :: p,detx,ps1,det1,psave(1),det(1),ps(1),aa,hk1,hk,  &
     lamdak,deltak,gk,root,root1,lamdk1,a,xlamsv,ptry,  &
     detry,srp,h1,h2,h3,gk1,hkp1,t1,t2,dist,dsave, temp2
 DIMENSION        ipdet(1)
 COMMON /detmx /  p(4),detx(4),ps1(4),det1(4),n2ev,ipsav,ips,idet,  &
     ipdeta,prec,nstart,u2,ic,l1,l2,is,nd,iadd,sml1,  &
     ipdetx(4),ipdet1(4),ifail,k,fact1,iffnd,nfail, npole
 COMMON /regean/  im(7),ik(7),iev(7),scr1,scr2,scr3,scr4,scr5,  &
     lcore,rmax,rmin,mz,nev,epsi,rminr,NE,nit,nevm, scr6,scr7,nfound,lama
 COMMON /zzzzzz/  psave
 EQUIVALENCE      (psave(1),ps(1),det(1),ipdet(1))
 
 CALL arrm (ps1,det1,ipdet1)
 aa = ps1(3) - ps1(2)
 dsave = 1.0E38
 
!     COPY INTO INTERATION BLOCK
 
 DO  n = 1,3
   detx(n) = det1(n)
   p(n) = ps1(n)
   ipdetx(n) = ipdet1(n)
 END DO
 
!     START INTERATION LOOP
 
 k     = 1
 igoto = 1
 40 hk1   = p(2) - p(1)
 hk    = p(3) - p(2)
 lamdak= hk/hk1
 IF (DABS(hk) <= DABS(epsi*100.0*p(3))) GO TO 240
 
!     CHECK FOR EARLY CONVERGENCE
 
 deltak = 1.0D0 + lamdak
 
!     COMPUTE  GK
 
 CALL summ (t1,it1,detx(1)*lamdak*lamdak,ipdetx,detx(2)*deltak*  &
     deltak,ipdetx(2),-1)
 CALL summ (gk,igk,t1,it1,detx(3)*(lamdak+deltak),ipdetx(3),1)
 
!     COMPUTE ROOT1
 
 CALL summ (t1,it1,detx(1)*lamdak,ipdetx(1),detx(2)*deltak, ipdetx(2),-1)
 CALL summ (t2,it2,t1,it1,detx(3),ipdetx(3),1)
 CALL summ (root1,iroot1,gk*gk,2*igk,-4.0*deltak*lamdak*detx(3)*t2,  &
     ipdetx(3)+it2,1)
 
!     COMPUTE ROOT = DSQRT (ROOT1)
 
 CALL sqrtm (root,iroot,root1,iroot1)
 a   = -2.0*detx(3)*deltak
 gk1 = gk
 DO  n = 1,2
   IF (root1 < 0.0) GO TO 50
   temp2 = root
   IF (gk1 /= 0.0D0) temp2 = DSIGN(root,gk1)
   
   CALL summ (t1,it1,gk,igk,temp2,iroot,1)
   
   lamdk1 = a/t1
   ilmk   = ipdetx(3) - it1
   lamdk1 = lamdk1*10.0**ilmk
   GO TO 60
   
!     T1= GK*GK + DABS(ROOT1)
   
   50 CALL summ (t1,it1,gk*gk,igk+igk,DABS(root1),iroot1,1)
   lamdk1 = a*gk/t1
   ilmk   = ipdetx(3) + igk - it1
   lamdk1 = lamdk1*10.0**ilmk
   GO TO 100
   60 IF (k /= 1) GO TO 100
   
!     IF (K .EQ. 1) RECALC LK1 TO MINIMIZE DIST
   
   dist = 0.0D0
   DO  i = 1,3
     dist = DABS(ps1(i)-ps1(3)-lamdk1*aa) + dist
   END DO
   IF (dist >= dsave) GO TO 80
   dsave  = dist
   xlamsv = lamdk1
   80 gk1 = -gk1
 END DO
 lamdk1 = xlamsv
 100 hkp1 = lamdk1*hk
 ptry = p(3) + hkp1
 
!     RANGE CHECKS
 
 IF (ptry > rmax) GO TO 120
 IF (is == n2ev-1) GO TO 110
 nnp = is + ips
 IF (ptry > 0.45*ps(nnp+2)+0.55*ps(nnp+3)) GO TO 120
 110 IF (ptry < rminr) GO TO 111
 GO TO 140
 
!     INCREASE POLE  AT LOWEST  E. V. GEOMETRICALLY
 
 111 npole1 = npole + 1
 npole  = 2*npole + 1
 
!     SWEEP PREVIOUSLY EVALUATED STARTING POINTS BY POLES
 
 n2ev2 = nd + iadd
 DO  n = 1,n2ev2
   nnd   = n + idet
   nnp   = n + ips
   nni   = n + ipdeta
   ptry  = 1.0D0
   iptry = 0
   DO  i = 1,npole1
     ptry = ptry*(ps(nnp)-rminr)
     CALL detm6 (ptry,iptry)
   END DO
   det(nnd)   = det(nnd)/ptry
   ipdet(nni) = ipdet(nni) - iptry
   CALL detm6 (det(nnd),ipdet(nni))
 END DO
 GO TO 120
 
!     NEW STARTING SET
 
 120 ifail = 0
 119 is    = is + 1
 IF (is  >= n2ev) GO TO 130
 IF (nstart == 0) iadd = iadd + 1
 RETURN 1
 
!      LOOK AT OLD STARTING SETS AGAIN
 
 130 IF (iffnd /= 1) RETURN 2
 iffnd  = 0
 is     = 1
 nstart = nstart + 1
 RETURN 1
 
!     TRY FOR CONVERGENCE
 
 140 CALL tmtogo (iptry)
 IF (iptry <= 0) RETURN 3
 CALL eadd (-ptry,prec)
 CALL detdet (detry,iptry,ptry,sml1,detx(3),ipdetx(3))
 IF (detry /= 0.0D0) GO TO 145
 igoto = 2
 GO TO 180
 
!     BEGIN CONVERGENCE TESTS
 
 145 IF (k <= 2) GO TO 170
 srp = DSQRT(DABS(p(3)))
 h1  = DABS(hk1)/srp
 h2  = DABS(hk)/srp
 h3  = DABS(hkp1)/srp
 150 fact1 = epsi*SQRT(rmax)
 IF (h1 > 2.e7*fact1) GO TO 200
 IF (h2 > 2.e4*fact1) GO TO 200
 IF (h3 > h2) GO TO 160
 IF (h3 > 2.*fact1) GO TO 200
 igoto = 2
 GO TO 180
 160 IF (h2 > 20.*fact1) GO TO 200
 igoto = 2
 GO TO 180
 
!     INTERATE AGAIN
 
 170 k = k + 1
 180 DO  i = 1,2
   p(i) = p(i+1)
   ipdetx(i) = ipdetx(i+1)
   detx(i)   = detx(i+1)
 END DO
 ipdetx(3) = iptry
 detx(3)   = detry
 p(3) = ptry
 SELECT CASE ( igoto )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 240
 END SELECT
 
!     FAIL TEST
 
 200 k = k + 1
 IF (k-nit < 0) THEN
   GO TO   180
 ELSE IF (k-nit == 0) THEN
   GO TO   210
 ELSE
   GO TO   220
 END IF
 210 IF (ifail == 1 .AND. ic < NE) GO TO 230
 220 ifail = 1
 nfail = nfail + 1
 GO TO 119
 230 epsi = 10.0*epsi
 ic = ic + 1
 GO TO 150
 
!     ACCEPT PK
 
 240 iffnd = 1
 
 RETURN
END SUBROUTINE detm3
