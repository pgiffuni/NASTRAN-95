SUBROUTINE detm1(*)
     
!     RMAX = APPROXIMATE MAGNITUDE OF LARGEST EIGENVALUE OF INTEREST
 
!     RMIN = LOWEST  NON-ZERO  EIGENVALUE
 
!     MZ = NUMBER OF ZERO EIGENVALUES
 
!     NEV = NUMBER OF NON-ZERO EIGENVALUES IN RANGE OF INTEREST
 
!     EPSI = CONVERGENCE CRITERION
 
!     RMINR = LOWEST EIGENVALUE OF INTEREST
 
!     NE   =  NUMBER OF PERMISSIBLE CHANGES OF EPSI
 
!     NIT = INTERATIONS TO AN EIGENVALUE
 
!     NEVM = MAXIMUM NUMBER OF EIGENVALUES DESIRED
 
!     IS  = STARTING SET COUNTER
 
!     IC  = COUNTER FOR CHANGE OF CONVERGENCE CRITERIA
 
!     NFOUND  = THE NUMBER OF EIGENVALUES FOUND TO DATA
!     NSTART = NUMBER OF TIMES THROUGH THE STARTING VALUES
 
 
!      IM = MASS MATRIX CONTROL BLOCK
 
!      IK = K MATRIX CONTROL BLOCK
 
!        A = M +P*K
 
!     IEV = EIGENVECTOR CONTROL BLOCK
 
 
 , INTENT(IN)                             :: *
 DOUBLE PRECISION :: p,detx,ps1,det1,psave(1),det(1),ps(1),fact
 DOUBLE PRECISION :: f1, f2, x, y, ratio
 
 INTEGER :: prec,u2,scr1,scr2,scr3,scr4,scr5,scr6,scr7,NAME(2)
 
 DIMENSION ipdet(1)
 
 COMMON /detmx/p(4),detx(4),ps1(4),det1(4),n2ev,ipsav,ips,idet,  &
     ipdeta,prec,nstart,u2,ic,l1,l2,is,nd,iadd,sml1,ipdetx(4),ipdet1(4)  &
     ,ifail,k,fact1,iffnd,nfail,npole,isng
 COMMON  /regean/ im(7),ik(7),iev(7),scr1,scr2,scr3,scr4,scr5,lcore  &
     , rmax,rmin,mz,nev,epsi,rminr,NE,nit,nevm,scr6,scr7 ,nfound,lama
 COMMON /zzzzzz/psave
 
 EQUIVALENCE  (psave(1),ps(1),det(1),ipdet(1))
 
 DATA NAME/4HDETM,4H1   /
 
! ----------------------------------------------------------------------
 
 ic = 0
 
!     CALCULATE THE NUMBER OF STARTING POINTS TO BE USED
 
 n2ev = 2*nev
 nn = n2ev
 srrmin = SQRT (rmin)
 srrmax = SQRT (rmax)
 fact = (srrmax - srrmin)/n2ev
 f1 = srrmin
 i = 0
 120 i = i + 1
 f2 = f1 + fact
 x = DLOG10 (f2/f1)
 IF (x < 1.0D0) GO TO 140
 ix = x
 y = ix
 IF (x /= y) ix = ix + 1
 n2ev = n2ev + ix - 1
 f1 = f2
 140 IF (i < nn) GO TO 120
 
!     CHECK AVAILABILITY OF CORE
 
 lc = 2*(korsz(psave)/2)
 ipsav = lc/2-nevm
 ips = ipsav -n2ev-1
 idet = ips-n2ev-1
 ipdeta = 2*idet -n2ev-2
 IF(ipdeta <= 0) GO TO 80
 lcore = lc-ipdeta+1
 
!     COMPUTE THE STARTING POINTS
 
 nn = ips + 1
 ps(nn) = rmin
 f1 = srrmin
 i = 0
 160 f2 = f1 + fact
 ratio = f2/f1
 x = DLOG10 (ratio)
 IF (x < 1.0D0) GO TO 200
 ix = x
 y = ix
 IF (x /= y) ix = ix + 1
 ratio = ratio**(1.0D0/ix)
 n = 0
 180 n = n + 1
 i = i + 1
 nn = nn + 1
 ps(nn) = ps(nn-1)*ratio*ratio
 IF (n < ix) GO TO 180
 GO TO 220
 200 i = i + 1
 nn = nn + 1
 ps(nn) = f2**2
 220 f1 = f2
 IF (i < n2ev) GO TO 160
 is=1
 nd=3
 iadd=0
 isng = 0
 rmax = 1.05*rmax
 fact1 = epsi*SQRT(rmax)
 
!     CALCULATE DETERMINANTE OF FIRST 3 STARTING VALUES
 
 ENTRY detm2
 IF(nstart /= 0) GO TO 40
 DO  n = 1, nd
   nn = n+iadd
   nnp = nn+ips
   nnd = nn+idet
   nni = nn+ipdeta
   CALL eadd(-ps(nnp),prec)
   CALL detdet(det(nnd),ipdet(nni),ps(nnp),sml1,0.0D0,1)
 END DO
 IF(nd == 3.AND.isng == 3)RETURN 1
 IF(is == 1) iadd=2
 nd = 1
 
!     CALCULATE THE INITAL GUESS
 
 
!     PERMUT VALUES TO ORDER BY DETERMINANT
 
 40 DO  n=1,3
   ns = n-1+is
   nnd = ns+idet
   nni = ns+ipdeta
   nnp = ns+ips
   det1(n) = det(nnd)
   ipdet1(n) = ipdet(nni)
   ps1(n) = ps(nnp)
 END DO
 RETURN
 80 CALL mesage (-8, 0, NAME)
 RETURN
END SUBROUTINE detm1
