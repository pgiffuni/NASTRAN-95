SUBROUTINE detm
     DOUBLE PRECISION :: p,detx,ps1,det1
 INTEGER :: prec,scr1,scr2,scr3,scr4,scr5,scr6,scr7
 INTEGER :: NAME(2)
 COMMON /detmx/p(4),detx(4),ps1(4),det1(4),n2ev,ipsav,ips,idet,  &
     ipdeta,prec,nstart,ndcmp,ic, nsmove,iterm,is,nd,iadd,sml1,  &
     ipdetx(4),ipdet1(4), ifail,k,fact1,iffnd,nfail , npole,isng
 COMMON  /regean/ im(7),ik(7),iev(7),scr1,scr2,scr3,scr4,scr5,lcore  &
     , rmax,rmin,mz,nev,epsi,rminr,NE,nit,nevm,scr6,scr7 ,nfound,lama
 DATA NAME/4HDETE,4HRM  /
 
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
 
!      IM = MASS MATRIX CONTROL BLOCK
 
!      IK = K MATRIX CONTROL BLOCK
 
!        A = M +P*K
 
!     IEV = EIGENVECTOR CONTROL BLOCK
 
 nstart =0
 lcore=0
 ndcmp = 0
 nsmove =0
 npole =0
 iterm = 1
 iffnd = 0
 nfail =0
!*****
 prec = ik(5)
!*****
 iscr7 = scr7
 IF(mz > nevm) GO TO 40
 IF (im(1) > 0) GO TO 5
 
!     MASS MATRIX PURGED -- ASSUME IDENTITY
 
 im(1) =  ik(1)
 CALL rdtrl(im(1))
 im(4) =8
 5 CONTINUE
 CALL detm1(*60)
 10 CALL klock(itime1)
 CALL detm3(*30,*40,*11)
 nfound = nfound+1
 CALL fdvect(sml1,p(3))
 idone=nfound+1
 IF(mz > 0)  idone=idone+mz
 CALL detm4
 IF(idone > nevm) GO TO 50
 CALL klock(itime2)
 CALL tmtogo(itleft)
 IF(2*(itime2-itime1) <= itleft) GO TO 10
 
!     INSUFFICIENT TIME TO FIND ANOTHER E. V.
 
 11 CONTINUE
 CALL mesage(45,nevm-idone,NAME)
 iterm = 3
 GO TO 50
 20 RETURN
 30 CALL detm2
 GO TO 10
 40 iterm = 2
 50 scr7 = iscr7
 CALL detm5
 GO TO 20
 
!     SINGULAR MATRIX EVERYWHERE
 
 60 iterm = 4
 GO TO 50
END SUBROUTINE detm
