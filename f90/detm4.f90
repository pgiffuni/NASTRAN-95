SUBROUTINE detm4
     DOUBLE PRECISION :: p,detx,ps1,det1,psave(1),det(1),ps(1)
 INTEGER :: prec,u1,u2,scr1,scr2,scr3,scr4,scr5,scr6,scr7
 DIMENSION ipdet(1)
 COMMON /detmx/p(4),detx(4),ps1(4),det1(4),n2ev,ipsav,ips,idet,  &
     ipdeta,prec, u1,u2,ic,nsmove,l2,is,nd, iadd,sml1,ipdetx(4),  &
     ipdet1(4), ifail,k,fact1
 COMMON  /regean/ im(7),ik(7),iev(7),scr1,scr2,scr3,scr4,scr5,lcore  &
     , rmax,rmin,mz,nev,epsi,rminr,NE,nit,nevm,scr6,scr7 ,nfound,lama
 COMMON /zzzzzz/psave
 EQUIVALENCE  (psave(1),ps(1),det(1),ipdet(1))
 
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
 
!     NSMOVE = THE NUMBER OF TIMES THE STATTING POINTS HAVE BEEN MOVED
 
!      IM = MASS MATRIX CONTROL BLOCK
 
!      IK = K MATRIX CONTROL BLOCK
 
!        A = M +P*K
 
!     IEV = EIGENVECTOR CONTROL BLOCK
 
 nn = ipsav+nfound
 psave(nn) = p(3)
 eps1 = fact1*DSQRT(DABS(p(3)))
 DO  n=1,3
   nn = n  + iadd -2
   nnp = nn+ips
   IF(DABS(ps(nnp)-p(3)) >= 400.*eps1) CYCLE
   10 ps(nnp) = ps(nnp) +2.e3*eps1
   nsmove = nsmove+1
   IF(nfound == 1) GO TO 30
   nfnd = nfound-1
   DO  i=1,nfnd
     nnz = ipsav+i
     IF(DABS(ps(nnp)-psave(nnz)) > 400.*eps1) CYCLE
     GO TO 10
   END DO
   30 nnd = nn+idet
   nni = nn+ipdeta
   CALL eadd(-ps(nnp),prec)
   CALL detdet(det(nnd),ipdet(nni),ps(nnp),sml1,0.0D0,1)
 END DO
 n2ev2 = iadd + nd
 DO  i=1,n2ev2
   nnd = i+idet
   nnp = i+ips
   nni = i+ipdeta
   det(nnd) = det(nnd)/(ps(nnp)-p(3))
   CALL detm6(det(nnd),ipdet(nni))
 END DO
 DO  i=1,3
   det1(i) = det1(i)/(ps1(i)-p(3))
 END DO
 RETURN
END SUBROUTINE detm4
