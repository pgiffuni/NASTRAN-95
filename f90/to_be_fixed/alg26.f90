SUBROUTINE alg26
     
 REAL :: loss,lami,lamip1,lamim1
 
 DIMENSION xx1(21),dsdm(21),drvwdm(21),dl(21),dsdl(21),fx1(21),fx2(  &
     21),vvold(21),afun(20),bfun(20),hs(20),xm2(20),dvmdvm(20),dvmdm(21  &
     ),tbip1(21),teip1(21)
 
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
 
 itmax=20
 lpmax=10
 k=1
 IF(i == istag)k=2
 xn=speed(i)*spdfac(icase)*pi/(30.0*sclfac)
 IF(i == 1)GO TO 234
 DO  j=1,nstrms
   lami(j)=lamip1(j)
   lamip1(j)=1.0
 END DO
 IF(i == nstns)GO TO 234
 IF(ndata(i+1) == 0)GO TO 210
 l1=ndimen(i+1)+1
 SELECT CASE ( l1 )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 150
   CASE (    4)
     GO TO 170
 END SELECT
 110   DO  j=1,nstrms
   xx1(j)=r(j,i+1)
 END DO
 GO TO 190
 130   DO  j=1,nstrms
   xx1(j)=r(j,i+1)/r(nstrms,i+1)
 END DO
 GO TO 190
 150   DO  j=1,nstrms
   xx1(j)=xl(j,i+1)
 END DO
 GO TO 190
 170   DO  j=1,nstrms
   xx1(j)=xl(j,i+1)/xl(nstrms,i+1)
 END DO
 190   l1=is2(i+1)
 CALL alg01(datac(l1),data4(l1),ndata(i+1),xx1,xx1,x1,nstrms,nterp (i+1),0)
 DO  j=1,nstrms
   lamip1(j)=1.0-xx1(j)
 END DO
 210   DO  j=1,nstrms
   x1=SQRT((r(j,i+1)-r(j,i))**2+(x(j,i+1)-x(j,i))**2)
   x2=SQRT((r(j,i)-rim1(j))**2+(x(j,i)-xim1(j))**2)
   x3=ATAN2(r(j,i+1)-r(j,i),x(j,i+1)-x(j,i))
   x4=ATAN2(r(j,i)-rim1(j),x(j,i)-xim1(j))
   phi(j)=(x3+x4)/2.0
   cr(j)=(x3-x4)/(x1+x2)*2.0
   dsdm(j)=0.0
   drvwdm(j)=0.0
   dvmdm(j)=0.0
   IF(ipass == 1)CYCLE
   dsdm(j)=((s(j,i+1)-s(j,i))/x1+(s(j,i)-s(j,i-1))/x2)/2.0*g*ej
   drvwdm(j)=((r(j,i+1)*vw(j,i+1)-r(j,i)*vw(j,i))/x1+(r(j,i)*vw(j,i)-  &
       rim1(j)*vw(j,i-1))/x2)/(2.0*r(j,i))
   dvmdm(j)=((vm(j,i+1)-vm(j,i))/x1+(vm(j,i)-vm(j,i-1))/x2)*0.5
 END DO
 IF(ipass == 1.OR.ndata(i) == 0.OR.neqn == 3.OR.nwork(i) /= 0.OR.nw  &
     ork(i+1) == 0)GO TO 390
 l1=ndimen(i)+1
 SELECT CASE ( l1 )
   CASE (    1)
     GO TO 221
   CASE (    2)
     GO TO 223
   CASE (    3)
     GO TO 225
   CASE (    4)
     GO TO 227
 END SELECT
 221   DO  j=1,nstrms
   teip1(j)=r(j,i)
 END DO
 GO TO 229
 223   DO  j=1,nstrms
   teip1(j)=r(j,i)/r(nstrms,i)
 END DO
 GO TO 229
 225   DO  j=1,nstrms
   teip1(j)=xl(j,i)
 END DO
 GO TO 229
 227   DO  j=1,nstrms
   teip1(j)=xl(j,i)/xl(nstrms,i)
 END DO
 229   l1=is2(i)
 CALL alg01(datac(l1),data3(l1),ndata(i),teip1,teip1,x1,nstrms,nte rp(i),0)
 x1=speed(i+1)*spdfac(icase)*pi/(30.0*sclfac)
 DO  j=1,nstrms
   teip1(j)=TAN(teip1(j)/c1)
   tbip1(j)=(vw(j,i)-x1*r(j,i))/vm(j,i)
 END DO
 GO TO 390
 234   DO  j=1,nstrms
   dvmdm(j)=0.0
   dsdm(j)=0.0
   drvwdm(j)=0.0
   cr(j)=0.0
 END DO
 IF(i == 1)GO TO 244
 DO  j=1,nstrms
   phi(j)=ATAN2(r(j,i)-rim1(j),x(j,i)-xim1(j))
 END DO
 GO TO 390
 244   DO  j=1,nstrms
   phi(j)=ATAN2(r(j,2)-r(j,1),x(j,2)-x(j,1))
 END DO
 DO  j=1,nstrms
   xi(j)=h(j,1)
   lami(j)=1.0
   lamip1(j)=1.0
 END DO
 IF(ndata(2) == 0)GO TO 390
 l2=ndimen(2)+1
 SELECT CASE ( l2 )
   CASE (    1)
     GO TO 290
   CASE (    2)
     GO TO 310
   CASE (    3)
     GO TO 330
   CASE (    4)
     GO TO 350
 END SELECT
 290   DO  j=1,nstrms
   xx1(j)=r(j,2)
 END DO
 GO TO 370
 310   DO  j=1,nstrms
   xx1(j)=r(j,2)/r(nstrms,2)
 END DO
 GO TO 370
 330   DO  j=1,nstrms
   xx1(j)=xl(j,2)
 END DO
 GO TO 370
 350   DO  j=1,nstrms
   xx1(j)=xl(j,2)/xl(nstrms,2)
 END DO
 370   l1=is2(2)
 CALL alg01(datac(l1),data4(l1),ndata(2),xx1,xx1,x1,nstrms,nterp(2 ),0)
 DO  j=1,nstrms
   lamip1(j)=1.0-xx1(j)
 END DO
 390   CALL alg01(r(1,i),x(1,i),nstrms,r(1,i),x1,gama,nstrms,0,1)
 DO  j=1,nstrms
   gama(j)=ATAN(gama(j))
   sppg(j)=gama(j)+phi(j)
   cppg(j)=COS(sppg(j))
   sppg(j)=SIN(sppg(j))
   vv(j)=vm(j,i)
 END DO
 DO  j=1,itub
   dl(j)=xl(j+1,i)-xl(j,i)
   dsdl(j)=(s(j+1,i)-s(j,i))*g*ej/dl(j)
 END DO
 IF(i == 1.OR.nwork(i) >= 5)GO TO 430
 DO  j=1,itub
   dvmdvm(j)=0.0
   fx1(j)=(vw(j+1,i)+vw(j,i))/(r(j+1,i)+r(j,i))*(r(j+1,i)*vw(j+1,i)-r  &
       (j,i)*vw(j,i))/dl(j)
   fx2(j)=(h(j+1,i)-h(j,i))/dl(j)*g*ej
 END DO
 GO TO 450
 430   DO  j=1,itub
   fx1(j)=(tbeta(j+1,i)+tbeta(j,i))/(r(j+1,i)+r(j,i))*(r(j+1,i)*tbeta  &
       (j+1,i)-r(j,i)*tbeta(j,i))/dl(j)
   fx2(j)=(xi(j+1)-xi(j))/dl(j)*g*ej
 END DO
 450   vmax=0.0
 vmin=2500.0
 iter=0
 460   iter=iter+1
 ifail=0
 iconf1=0
 DO  j=1,nstrms
   vvold(j)=vv(j)
 END DO
 IF(i == 1.OR.nwork(i) >= 5)GO TO 810
 DO  j=1,itub
   x1=(h(j,i)+h(j+1,i))/2.0-(((vvold(j)+vvold(j+1))/2.0)**2+((vw(j,i)  &
       +vw(j+1,i))/2.0)**2)/(2.0*g*ej)
   IF(x1 >= hmin)GO TO 520
   IF(ipass <= nforce)GO TO 510
   IF(lnct < npage)GO TO 480
   WRITE(log2,500)
   lnct=1
   480   lnct=lnct+1
   WRITE(log2,490)ipass,i,iter,j,x1
   490   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamtu  &
       be,i3,53H  static enthalpy below limit in momentum equation at,e13 .5)
   500   FORMAT(1H1)
   510   ifail=1
   x1=hmin
   520   x2=(s(j,i)+s(j+1,i))/2.0
   x7=alg7(x1,x2)
   x2=(cppg(j)+cppg(j+1))*0.5
   x3=(sppg(j)+sppg(j+1))*0.5
   afun(j)=-2.0*x3*(dvmdm(j)+dvmdm(j+1))/(vvold(j)+vvold(j+1))-x2*(cr  &
       (j)+cr(j+1))
   bfun(j)=2.0*(fx2(j)-x7*dsdl(j)-fx1(j))
   IF(ipass == 1.OR.i == nstns)GO TO 580
   IF(ndata(i) == 0.OR.neqn == 3.OR.(nwork(i) == 0.AND.nwork(i+1) ==  &
       0))GO TO 560
   IF(nwork(i) == 0)GO TO 540
   x4=(tbeta(j,i)+tbeta(j+1,i))*0.5
   x5=(taneps(j)+taneps(j+1))*0.5
   530   bfun(j)=bfun(j)+x7*(dsdm(j)+dsdm(j+1))*(x3/(1.0+x4*x4)-x5*x4/(1.0+  &
       x4*x4)*0.5)-x5*(drvwdm(j)+drvwdm(j+1))*(vvold(j)+vvold(j+1))*0.5
   GO TO 580
   540   x4=(tbip1(j)+tbip1(j+1))*0.5
   x5 = (teip1(j)+teip1(j+1))*0.5
   GO TO 530
   560   bfun(j)=bfun(j)+x7*(dsdm(j)+dsdm(j+1))*x3
   580   vv(imid)=vvold(imid)**2
 END DO
 j=imid
 jinc=1
 590   jold=j
 j=j+jinc
 jj=jold
 IF(jinc == -1)jj=j
 IF(ABS(afun(j)) <= 1.0E-5) GO TO 660
 x1=-afun(jj)*(xl(j,i)-xl(jold,i))
 IF(ABS(x1) <= 1.0E-10)GO TO 660
 IF(x1 <= 88.0)GO TO 630
 IF(ipass <= nforce)GO TO 620
 IF(lnct < npage)GO TO 600
 WRITE(log2,500)
 lnct=1
 600   lnct=lnct+1
 WRITE(log2,610)ipass,i,iter,jj,x1
 610   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamtu  &
     be,i3,43H  momentum equation exponent above limit at,e13.5)
 620   ifail=1
 x1=88.0
 630   x1=EXP(x1)
 vv(j)=vv(jold)*x1+(1.0-x1)*bfun(jj)/afun(jj)
 640   IF(j == k)GO TO 670
 IF(j == nstrms)GO TO 650
 GO TO 590
 650   j=imid
 jinc=-1
 GO TO 590
 660   vv(j)=vv(jold)+bfun(jj)*(xl(j,i)-xl(jold,i))
 GO TO 640
 670   DO  j=k,nstrms
   IF(vv(j) <= 4.0*vvold(imid)**2)GO TO 676
   ifail=1
   IF(ipass <= nforce)GO TO 674
   CALL alg03(lnct,1)
   WRITE(log2,672)ipass,i,iter,j
   672   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
       NE,i3,50H  meridional velocity greater than twice mid value)
   674   vv(j)=4.0*vvold(imid)**2
   676   IF(vv(j) >= 1.0)GO TO 702
   IF(ipass <= nforce)GO TO 700
   IF(lnct < npage)GO TO 680
   WRITE(log2,500)
   lnct=1
   680   lnct=lnct+1
   WRITE(log2,690)ipass,i,iter,j,vv(j)
   690   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
       NE,i3,46H  (meridional velocity) squared below limit at,e13.5)
   700   vv(j)=1.0
   ifail=1
   CYCLE
   702   vv(j)=SQRT(vv(j))
 END DO
 x1=0.0
 DO  j=k,itub
   x1=x1+(xl(j+1,i)-xl(j,i))*ABS((vv(j+1)+vv(j))/(vvold(j+1)+vvold(j) )-1.0)
 END DO
 x1=x1/(xl(nstrms,i)-xl(k,i))
 x2=0.1
 IF(x1 < 0.2)x2=EXP(-11.52*x1)
 DO  j=k,nstrms
   vv(j)=vvold(j)+x2*(vv(j)-vvold(j))
 END DO
 IF(nloss(i) == 1.AND.nl2(i) == 0)CALL alg07
 DO  j=1,itub
   hs(j)=(h(j,i)+h(j+1,i))/2.0-(((vv(j)+vv(j+1))/2.0)**2+((vw(j,i)+vw  &
       (j+1,i))/2.0)**2)/(2.0*g*ej)
   IF(hs(j) >= hmin)GO TO 800
   IF(ipass <= nforce)GO TO 790
   IF(lnct < npage)GO TO 770
   WRITE(log2,500)
   lnct=1
   770   lnct=lnct+1
   WRITE(log2,780)ipass,i,iter,j,hs(j)
   780   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamtu  &
       be,i3,55H  static enthalpy below limit in continuity equation at,e 13.5)
   790   ifail=1
   hs(j)=hmin
   800   xm2(j)=alg9(hs(j),(s(j,i)+s(j+1,i))/2.0,((vv(j)+vv(j+1))/2.0)**2)
 END DO
 GO TO 1100
 810   j=imid
 jinc=1
 820   loop=1
 jold=j
 j=j+jinc
 jj=jold
 IF(jinc == -1)jj=j
 830   vold=vv(j)
 vav=(vold+vv(jold))/2.0
 ifaie=0
 iconf2=0
 x2=(tbeta(j,i)+tbeta(jold,i))/2.0
 x1=(xi(j)+xi(jold))/2.0+((xn*(r(j,i)+r(jold,i))/2.0)**2-vav**2*(1.  &
     0+x2*x2))/(2.0*g*ej)
 IF(x1 >= hmin)GO TO 870
 IF(ipass <= nforce)GO TO 860
 IF(lnct < npage)GO TO 840
 WRITE(log2,500)
 lnct=1
 840   lnct=lnct+1
 WRITE(log2,850)ipass,i,iter,jj,loop,x1
 850   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamtu  &
     be,i3,6H  loop,i3,43H  static h in momentum equn. below limit at,e 13.5)
 860   ifaie=1
 iconf2 = 1
 x1=hmin
 870   x3=(s(j,i)+s(jold,i))/2.0
 x7=alg7(x1,x3)
 x4=(sppg(j)+sppg(jold))*0.5
 x5=(cppg(j)+cppg(jold))*0.5
 x1=x5*(cr(j)+cr(jold))*0.5-fx1(jj)
 x12=1.0/(1.0+x2*x2)
 x8=(taneps(j)+taneps(jold))*0.5
 x11=fx2(jj)-x7*dsdl(jj)
 x6=x4*(dvmdm(j)+dvmdm(jold))*0.5-2.0*xn*x2*COS((gama(j)+gama(jold) )*0.5)
 IF(ipass == 1.OR.i == 1.OR.i == nstns)GO TO 920
 IF(neqn == 3)GO TO 900
 x11=x11+x7*(dsdm(j)+dsdm(jold))*0.5*(x4*x12-x8*x2*x12)
 x6=x6-x8*(drvwdm(j)+drvwdm(jold))*0.5
 GO TO 920
 900   x11=x11+x7*(dsdm(j)+dsdm(jold))*0.5*x4
 920   dv2dl=2.0*x12*(vav*(x6+vav*x1)+x11)
 dvmdvm(jj)=x12*(x1-x11/vav**2)
 x1=vv(jold)**2+dv2dl*(xl(j,i)-xl(jold,i))
 IF(x1 <= 9.0*vvold(imid)**2)GO TO 938
 iconf2=1
 ifaie=1
 IF(ipass <= nforce)GO TO 936
 CALL alg03(lnct,1)
 x1=SQRT(x1)
 x2=3.0*vvold(imid)
 WRITE(log2,934)ipass,i,iter,j,loop,x1,x2
 934   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
     NE,i3,6H  loop,i3,33H  meridional velocity above limit,e13.5,9H  l  &
     imit =,e13.5)
 936   x1=9.0*vvold(imid)**2
 938   IF(x1 >= 1.0)GO TO 950
 IF(ipass <= nforce)GO TO 944
 IF(lnct < npage)GO TO 930
 WRITE(log2,500)
 lnct=1
 930   lnct=lnct+1
 WRITE(log2,940)ipass,i,iter,j ,loop,x1
 940   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
     NE,i3,6H  loop,i3,46H  (meridional velocity) squared below limit a t,e13.5)
 944   x1=1.0
 ifaie=1
 iconf2=1
 950   vv(j)=SQRT(x1)
 IF(ABS(vv(j)/vold-1.0) <= tolnce/5.0)GO TO 990
 IF(loop >= lpmax)GO TO 960
 loop=loop+1
 GO TO 830
 960   iconf2=1
 IF(ipass <= nforce)GO TO 990
 IF(lnct < npage)GO TO 970
 WRITE(log2,500)
 lnct=1
 970   lnct=lnct+1
 WRITE(log2,980)ipass,i,iter,j,vv(j),vold
 980   FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
     NE,i3,38H  meridional velocity unconverged  vm=,e13.6,9H vm(OLD)=, e13.6)
 990   IF(ifaie == 1)ifail=1
 IF(iconf2 == 1)iconf1=1
 IF(j == nstrms)GO TO 1000
 IF(j == 1)GO TO 1010
 GO TO 820
 1000  j=imid
 jinc=-1
 GO TO 820
 1010  IF(i == 1)GO TO 1014
 IF(nloss(i) == 2.OR.(nloss(i) == 1.AND.nl2(i) == 0))CALL alg07
 1014  DO  j=1,itub
   x1=((vv(j)+vv(j+1))/2.0)**2*(1.0+((tbeta(j,i)+tbeta(j+1,i))/2.0)** 2)
   hs(j)=(xi(j)+xi(j+1))/2.0+((xn*(r(j,i)+r(j+1,i))/2.0)**2-x1)/(2.0* g*ej)
   IF(hs(j) >= hmin)GO TO 1080
   IF(ipass <= nforce)GO TO 1070
   IF(lnct < npage)GO TO 1060
   WRITE(log2,500)
   lnct=1
   1060  lnct=lnct+1
   WRITE(log2,780)ipass,i,iter,j,hs(j)
   1070  ifail=1
   hs(j)=hmin
   1080  xm2(j)=alg9(hs(j),(s(j,i)+s(j+1,i))/2.0,x1)
   IF(i == 1.OR.nloss(i) /= 1.OR.nl2(i) /= 0)CYCLE
   x1=(s(j,i)+s(j+1,i))/2.0
   x2=alg4(hs(j),x1)
   x4=alg8(hs(j),x1)
   x3=(xi(j)+xi(j))/2.0+(xn*((r(j,i)+r(j+1,i))/2.0))**2/(2.0*g*ej)
   x3=alg4(x3,x1)
   xm2(j)=xm2(j)*(1.0+x4*(loss(j)+loss(j+1))/2.0*x2/(x3*(1.0+(loss(j)  &
       +loss(j+1))/2.0*(1.0-x2/x3))))
 END DO
 1100  delw(1)=0.0
 dwdv=0.0
 x2=bblock(i)*bdist(i)
 x3=bblock(i)*(1.0-bdist(i))/xl(nstrms,i)
 DO  j=1,itub
   x1=dl(j)*(r(j+1,i)+r(j,i))*alg5(hs(j),(s(j,i)+s(j+1,i))/2.0)*(vv(j  &
       )+vv(j+1))*(cppg(j)+cppg(j+1))*pi/(4.0*sclfac**2)
   x1=x1*((lami(j)+lami(j+1))/2.0-wwbl(i)-x2-x3*(xl(j,i)+xl(j+1,i)))
   delw(j+1)=delw(j)+x1
   x4=0.0
   IF(j >= imid)GO TO 1130
   l1=j
   1110  x4=x4+dvmdvm(l1)
   IF(l1 >= imid-1)GO TO 1120
   l1=l1+1
   GO TO 1110
   1120  x4=x4/FLOAT(imid-j)
   GO TO 1200
   1130  l1=imid+1
   1140  x4=x4+dvmdvm(l1)
   IF(l1 >= j)GO TO 1150
   l1=l1+1
   GO TO 1140
   1150  x4=x4/FLOAT(j-imid+1)
   1200  dwdv=dwdv+x1*(1.0-xm2(j))*2.0/((vv(j)+vv(j+1))*(1.0-((xl(j,i)+xl(j  &
       +1,i))*0.5-xl(imid,i))*x4))
 END DO
 w=delw(nstrms)
 fm2=dwdv/w*vv(imid)
 DO  j=2,nstrms
   delw(j)=delw(j)/w
 END DO
 IF(dwdv <= 0.0)GO TO 1280
 IF(nmach(i) == 1)GO TO 1330
 IF(w < flow(icase).AND.iconf1 == 0)vmax=vv(imid)
 1220  dv=(flow(icase)-w)/dwdv
 IF(dv < -0.1*vv(imid))dv=-0.1*vv(imid)
 IF(dv > 0.1*vv(imid))dv= 0.1*vv(imid)
 1230  IF(ipass == 1.OR.(i /= 1.AND.nwork(i) <= 4))GO TO 1234
 IF(vv(imid)+dv < vmin)GO TO 1232
 dv=(vmin-vv(imid))*0.5
 1232  IF(vv(imid)+dv > vmax)GO TO 1234
 dv=(vmax-vv(imid))*0.5
 1234  DO  j=k,nstrms
   vv(j)=vv(j)+dv
   IF(vv(j) >= 1.0)CYCLE
   IF(ipass <= nforce)GO TO 1260
   IF(lnct < npage)GO TO 1240
   WRITE(log2,500)
   lnct=1
   1240  lnct=lnct+1
   WRITE(log2,1250)ipass,i,iter,j,vv(j)
   1250  FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,12H  streamli  &
       NE,i3,50H  meridional velocity below limit in continuity at,e13.5)
   1260  vv(j)=1.0
   ifail=1
 END DO
 GO TO 1340
 1280  IF(nmach(i) == 0)GO TO 1290
 IF(w < flow(icase).AND.iconf1 == 0)vmin=vv(imid)
 GO TO 1220
 1290  IF(vv(imid) < vmin.AND.iconf1 == 0)vmin=vv(imid)
 dv=-.1*vv(imid)
 1300  ifail=1
 IF(ipass <= nforce)GO TO 1230
 IF(lnct < npage)GO TO 1310
 WRITE(log2,500)
 lnct=1
 1310  lnct=lnct+1
 WRITE(log2,1320)ipass,i,iter
 1320  FORMAT(5X,4HPASS,i3,9H  station,i3,11H  iteration,i3,43H  other co  &
     ntinuity equation branch required)
 GO TO 1230
 1330  IF(vv(imid) > vmax.AND.iconf1 == 0)vmax=vv(imid)
 dv=0.1*vv(imid)
 GO TO 1300
 1340  x1=tolnce/5.0
 IF(neval(i) > 0)x1=x1/2.0
 IF(ABS(w/flow(icase)-1.0) > x1)GO TO 1354
 DO  j=k,nstrms
   IF(ABS(vv(j)/vvold(j)-1.0) > x1)GO TO 1354
 END DO
 GO TO 1390
 1354  IF(iter >= itmax)GO TO 1360
 IF(i == 1)GO TO 460
 IF((nloss(i) == 1.AND.nl2(i) == 0).OR.(nwork(i) >= 5.AND.nloss(i).  &
     EQ.2))CALL alg07
 GO TO 460
 1360  IF(ipass <= nforce)GO TO 1390
 IF(lnct < npage)GO TO 1370
 WRITE(log2,500)
 lnct=1
 1370  lnct=lnct+1
 x1=w/flow(icase)
 x2=vv(k)/vvold(k)
 x3=vv(imid)/vvold(imid)
 x4=vv(nstrms)/vvold(nstrms)
 WRITE(log2,1380)ipass,i,x1,x2,x3,x4
 1380  FORMAT(5X,4HPASS,i3,9H  station,i3,49H  momentum AND/OR continuity  &
     unconverged w/wspec=,f8.5,16H vm/vm(OLD) hub=,f8.5,5H mid=,f8.5,5  &
     h tip=,f8.5)
 1390  IF(ifail /= 0.AND.ifailo == 0)ifailo=i
 DO  j=1,nstrms
   vm(j,i)=vv(j)
 END DO
 IF(i /= 1)GO TO 1420
 DO  j=1,nstrms
   vw(j,1)=vv(j)*tbeta(j,1)
 END DO
 GO TO 1480
 1420  IF(nmix /= 1)GO TO 1440
 DO  j=1,nstrms
   s(j,i-1)=skeep(j)
   h(j,i-1)=hkeep(j)
   vw(j,i-1)=vwkeep(j)
 END DO
 1440  IF(nwork(i) >= 5)GO TO 1460
 tbeta(1,i)=0.0
 DO  j=k,nstrms
   tbeta(j,i)=(vw(j,i)-xn*r(j,i))/vv(j)
 END DO
 GO TO 1480
 1460  DO  j=1,nstrms
   vw(j,i)=vv(j)*tbeta(j,i)+xn*r(j,i)
   h(j,i)=xi(j)+xn*r(j,i)*vw(j,i)/(g*ej)
 END DO
 1480  CONTINUE
 RETURN
END SUBROUTINE alg26
