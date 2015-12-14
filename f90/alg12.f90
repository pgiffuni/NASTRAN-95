SUBROUTINE alg12
     
 REAL :: loss,lami,lamip1,lamim1
 
 DIMENSION pstat(32),xx(32)
 
 COMMON /ud300c/ nstns,nstrms,nmax,nforce,nbl,ncase,nsplit,nread,  &
     npunch,npage,nset1,nset2,istag,icase,ifailo,ipass,i,ivfail,iffail,  &
     nmix,ntrans,nplot,iloss,lnct,itub,imid,ifail,iter,log1,log2,log3,  &
     log4,log5,log6,iprint,nmany,nstplt,neqn,nspec(30),nwork(30),  &
     nloss(30),ndata(30),nterp(30),nmach(30),nl1(30),nl2(30),ndimen(30)  &
     ,is1(30),is2(30),is3(30),neval(30),ndiff(4),ndel(30),nliter(30),  &
     nm(2),nrad(2),ncurve(30),nwhich(30),nout1(30),nout2(30),nout3(30),  &
     nblade(30),dm(11,5,2),wfrac(11,5,2),r(21,30),xl(21,30),x(21,30),  &
     h(21,30),s(21,30),vm(21,30),vw(21,30),tbeta(21,30),diff(15,4),  &
     fdhub(15,4),fdmid(15,4),fdtip(15,4),terad(5,2),datac(100),  &
     data1(100),data2(100),data3(100),data4(100),data5(100),data6(100),  &
     data7(100),data8(100),data9(100),flow(10),speed(30),spdfac(10),  &
     bblock(30),bdist(30),wblock(30),wwbl(30),xstn(150),rstn(150),  &
     delf(30),delc(100),delta(100),title(18),drdm2(30),rim1(30),  &
     xim1(30),work(21),loss(21),taneps(21),xi(21),vv(21),delw(21),  &
     lami(21),lamim1(21),lamip1(21),phi(21),cr(21),gama(21),sppg(21),  &
     cppg(21),hkeep(21),skeep(21),vwkeep(21),delh(30),delt(30),visk,  &
     shape,sclfac,ej,g,tolnce,xscale,pscale,plow,rlow,xmmax,rconst,  &
     fm2,hmin,c1,pi,contr,conmx
 
 xmax=x(1,nstns)
 xmin=x(1,1)
 DO  j=2,nstrms
   IF(x(j,1) < xmin)xmin=x(j,1)
   IF(x(j,nstns) > xmax)xmax=x(j,nstns)
 END DO
 IF(xmin < 0.0)xmin=xmin-1.0
 l1=xmin-1.0
 xmin=FLOAT(l1)
 l1=xmax+1.0
 xmax=FLOAT(l1)
 delx=(xmax-xmin)/xscale+0.01
 xx(nstns+1)=xmin
 xx(nstns+2)=xscale
 IF(nplot == 2)GO TO 134
 pstat(nstns+1)=plow
 pstat(nstns+2)=pscale
 j=1
 k=1
 110   DO  i=1,nstns
   hs=h(j,i)-(vw(j,i)**2+vm(j,i)**2)/(2.0*g*ej)
   IF(hs < hmin)hs=hmin
   pstat(i)=alg4(hs,s(j,i))/sclfac**2
   xx(i)=x(j,i)
 END DO
 IF(j == nstrms)GO TO 130
 k=k+1
 IF(j == imid)j=nstrms
 IF(j == 1)j=imid
 GO TO 110
 130   CONTINUE
 IF(nplot == 1)GO TO 180
 134   CONTINUE
 pstat(nstns+1)=rlow
 pstat(nstns+2)=xscale
 DO  j=1,nstrms
   DO  i=1,nstns
     xx(i)=x(j,i)
     pstat(i)=r(j,i)
   END DO
 END DO
 pstat(nstrms+1)=rlow
 pstat(nstrms+2)=xscale
 xx(nstrms+1)=xmin
 xx(nstrms+2)=xscale
 DO  i=1,nstns
   DO  j=1,nstrms
     pstat(j)=r(j,i)
     xx(j)=x(j,i)
   END DO
 END DO
 180   RETURN
END SUBROUTINE alg12
