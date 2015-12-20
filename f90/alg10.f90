SUBROUTINE alg10
     
 REAL :: loss,lami,lamip1,lamim1
 
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
 
 IF(i > 1)GO TO 130
 v5=visk**0.2
 vinh=0.0
 vint=0.0
 IF(wwbl(1) > 0.0)GO TO 100
 c1h=0.0
 c1t=0.0
 delh(1)=0.0
 delt(1)=0.0
 GO TO 150
 100   IF(istag > 0)GO TO 110
 x1=wwbl(1)*xl(nstrms,1)*(cppg(1)+cppg(nstrms))/4.0
 delh(1)=x1
 delt(1)=x1
 x1=x1/(sclfac*shape)
 c1h=(x1*vm(1,1)**3.4/v5)**1.25
 c1t=(x1*vm(nstrms,1)**3.4/v5)**1.25
 GO TO 150
 110   delh(1)=0.0
 c1h=0.0
 IF(ABS(phi(nstrms)) > pi/2.0-0.00015 .AND. ABS(phi(nstrms)) < &
      pi/2.0+0.00015)GO TO 120
 x1=(r(nstrms,1)-SQRT(r(nstrms,1)**2-COS(phi(nstrms))*cppg(nstrms)*  &
     wwbl(1)*(r(nstrms,1)+r(1,1))*xl(nstrms,i)))/COS(phi(nstrms))
 delt(1)=x1
 c1t=(x1/(shape*sclfac*v5)*vm(nstrms,1)**3.4)**1.25
 GO TO 150
 120   delt(1)=wwbl(1)*xl(nstrms,1)/cppg(nstrms)
 c1t=(delt(1)*vm(nstrms,1)**3.4/(v5*sclfac*shape))**1.25
 GO TO 150
 130   vint=vint+SQRT((x(nstrms,i)-x(nstrms,i-1))**2+(r(nstrms,i)&
           -r(nstrms,i-1))**2)*((vm(nstrms,i)+vm(nstrms,i-1))&
           /2.0)**4/sclfac
 delt(i)=v5*(c1t+0.016*vint)**0.8/vm(nstrms,i)**3.4*sclfac*shape
 delh(i)=0.0
 IF(i <= istag)GO TO 140
 vinh=vinh+SQRT((x(1,i)-x(1,i-1))**2+(r(1,i)-r(1,i-1))**2)*((vm(1,i  &
     )+vm(1,i-1))/2.0)**4/sclfac
 delh(i)=v5*(c1h+0.016*vinh)**0.8/vm(1,i)**3.4*sclfac*shape
 140   wwbl(i)=0.5*wwbl(i)+0.5*(((2.0*r(nstrms,i)-delt(i)*COS(phi(nstrms)))&
              *delt(i)/cppg(nstrms)+(2.0*r(1,i)+delh(i)*COS(phi(1)))*delh(i)&
              /cppg(1))/((r(nstrms,i)+r(1,i))*xl(nstrms,i)))
 IF(wwbl(i) > 0.3)wwbl(i)=0.3
 IF(wwbl(i) < 0.0)wwbl(i)=0.3
 150   CONTINUE
 
 RETURN
END SUBROUTINE alg10
