SUBROUTINE alg09
     
 REAL :: loss,lami,lamip1,lamim1
 
 DIMENSION xx1(21),xx2(21),xx3(21),xx4(21),xx5(21),xx6(21),sol(21),  &
     wpara(21),wd(21),dif(21),ws(21),pm1(21),xinc(21),beta1(21),talph1(  &
     21),ang(21),highm(21),wt(21),xmr(21)
 
 COMMON /ud3prt/ iprtc
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
 
 wmax=0.7
 l1=i+nl1(i)
 xn=speed(i)*spdfac(icase)*pi/(30.0*sclfac)
 IF(iprint == 0)GO TO 116
 l2=ABS(FLOAT(neval(i)))
 CALL alg03(lnct,7+nstrms)
 lnct=lnct-3
 IF(neval(i) > 0.AND.iprtc == 1) WRITE(log2,100) l1,i,l2
 IF(neval(i) < 0.AND.iprtc == 1) WRITE(log2,110) l1,i,l2
 100   FORMAT(2X,/,8X,57HLOSS coefficient determination for blade between  &
     &stations,i3,4H AND,i3,47H - as incorporated in above results  bla&
     &de TYPE,i2,/,8X,116(1H*),/,2X)
 110   FORMAT(2X,/,8X,57HLOSS coefficient determination for blade between &
     &stations,i3,4H AND,i3,47H - for purposes of comparison only   bla&
     &de TYPE,i2,/,8X,116(1H*),/,2X)
 116   l2=ndimen(i)+1
 SELECT CASE ( l2 )
   CASE (    1)
     GO TO 120
   CASE (    2)
     GO TO 140
   CASE (    3)
     GO TO 160
   CASE (    4)
     GO TO 180
 END SELECT
 120   DO  j=1,nstrms
   xx2(j)=r(j,l1)
   xx6(j)=r(j,i)
 END DO
 GO TO 200
 140   DO  j=1,nstrms
   xx2(j)=r(j,l1)/r(nstrms,l1)
   xx6(j)=r(j,i)/r(nstrms,i)
 END DO
 GO TO 200
 160   DO  j=1,nstrms
   xx2(j)=xl(j,l1)
   xx6(j)=xl(j,i)
 END DO
 GO TO 200
 180   DO  j=1,nstrms
   xx2(j)=xl(j,l1)/xl(nstrms,l1)
   xx6(j)=xl(j,i)/xl(nstrms,i)
 END DO
 200   l2=is2(i)
 CALL alg01(datac(l2),data5(l2),ndata(i),xx6,sol,x1,nstrms,nterp(i ),0)
 q=1.0
 IF(speed(i) < 0.0)GO TO 208
 IF(speed(i) > 0.0)GO TO 206
 IF(i < 3)GO TO 208
 ii=i-1
 204   IF(speed(ii) /= 0.0)GO TO 205
 IF(ii == 2)GO TO 208
 ii=ii-1
 GO TO 204
 205   IF(speed(ii) < 0.0)q=-1.0
 GO TO 208
 206   q=-1.0
 208   DO  j=1,nstrms
   talph1(j)=(vw(j,l1)-xn*r(j,l1))/vm(j,l1)
   dif(j)=1.0-vm(j,i)/vm(j,l1)*SQRT((1.0+tbeta(j,i)**2)/(1.0+talph1(j)**2))&
       +(vm(j,l1)*talph1(j)-vm(j,i)*tbeta(j,i))/(2.0*sol(j)*vm(j,l1) &
       *SQRT(1.0+talph1(j)**2))*q
 END DO
 l2=ABS(FLOAT(neval(i)))
 l3=ndiff(l2)
 CALL alg01(diff(1,l2),fdhub(1,l2),l3,dif,xx3,x1,nstrms,0,0)
 CALL alg01(diff(1,l2),fdmid(1,l2),l3,dif,xx4,x1,nstrms,0,0)
 CALL alg01(diff(1,l2),fdtip(1,l2),l3,dif,xx5,x1,nstrms,0,0)
 xx1(1)=0.1
 xx1(2)=0.5
 xx1(3)=0.9
 DO  j=1,nstrms
   xx1(4)=xx3(j)
   xx1(5)=xx4(j)
   xx1(6)=xx5(j)
   x1=(r(j,i)-r(1,i))/(r(nstrms,i)-r(1,i))
   CALL alg01(xx1,xx1(4),3,x1,wpara(j),x1,1,0,0)
 END DO
 DO  j=1,nstrms
   xmr(j)=0.0
   highm(j)=0.0
   ang(j)=0.0
   ws(j)=0.0
   xinc(j)=0.0
   beta1(j)=0.0
   wd(j)=wpara(j)*2.0*sol(j)*SQRT(1.0+tbeta(j,i)**2)
   wt(j)=wd(j)
 END DO
 IF(ndel(i) == 0)GO TO 384
 l2=is3(i)
 CALL alg01(delc(l2),delta(l2),ndel(i),xx2,pm1,x1,nstrms,1,0)
 IF(ndata(l1) == 0)GO TO 340
 CALL alg01(r(1,l1),x(1,l1),nstrms,r(1,l1),x1,xx1,nstrms,0,1)
 l2=ndimen(l1)+1
 SELECT CASE ( l2 )
   CASE (    1)
     GO TO 240
   CASE (    2)
     GO TO 260
   CASE (    3)
     GO TO 280
   CASE (    4)
     GO TO 300
 END SELECT
 240   DO  j=1,nstrms
   xx2(j)=r(j,l1)
 END DO
 GO TO 320
 260   DO  j=1,nstrms
   xx2(j)=r(j,l1)/r(j,nstrms)
 END DO
 GO TO 320
 280   DO  j=1,nstrms
   xx2(j)=xl(j,l1)
 END DO
 GO TO 320
 300   DO  j=1,nstrms
   xx2(j)=xl(j,l1)/xl(nstrms,l1)
 END DO
 320   l2=is2(l1)
 l3=ndata(l1)
 CALL alg01(datac(l2),data1(l2),l3,xx2,xx3,x1,nstrms,nterp(l1),0)
 CALL alg01(datac(l2),data3(l2),l3,xx2,xx4,x1,nstrms,nterp(l1),0)
 DO  j=1,nstrms
   x1=(ATAN((r(j,l1+1)-r(j,l1))/(x(j,l1+1)-x(j,l1)))+ATAN((r(j,l1)&
     -r(j,l1-1))/(x(j,l1)-x(j,l1-1))))/2.0
   beta1(j)=ATAN((TAN(xx3(j)/c1)*(1.0-xx1(j)*TAN(x1))&
           -TAN(x1)*TAN(xx4(j)/c1)*SQRT(1.0+xx1(j)**2))*COS(x1))
   xinc(j)=(ATAN(talph1(j))-beta1(j))*q
 END DO
 340   DO  j=1,nstrms
   ang(j)=xinc(j)+pm1(j)/c1
   x1=h(j,l1)-(vm(j,l1)**2+vw(j,l1)**2)/(2.0*g*ej)
   IF(x1 < hmin)x1=hmin
   x4=alg8(x1,s(j,l1))
   x2=(x4+1.0)/(x4-1.0)
   x3=SQRT(x2)
   x5=alg9(x1,s(j,l1),vm(j,l1)**2*(1.0+talph1(j)**2))
   xmr(j)=SQRT(x5)
   x6=x5
   IF(x6 < 1.0)x6=1.0
   x7=x3*ATAN(SQRT(x6-1.0)/x3)-ATAN(SQRT(x6-1.0))+ang(j)
   x10=0.0
   IF(x7 <= 0.0)GO TO 376
   x8=0.4*pi*(x3-1.0)
   IF(x7 > x8)GO TO 374
   x9 = 1.0
   k=1
   350   x10=x9-(x2+x9*x9)*(1.0+x9*x9)/(x9*x9*(x2-1.0))&
            *(x3*ATAN(x9/x3)-atan(x9)-x7)
   IF(ABS(x10-x9) <= 0.00001)GO TO 376
   IF(k > 20)GO TO 360
   k=k+1
   x9=x10
   GO TO 350
   360   IF(iprint == 0)GO TO 374
   CALL alg03(lnct,1)
   WRITE(log2,370)ipass,i,j
   370   FORMAT(5X,4HPASS,i3,9H  station,i3,12H  streamline,i3,58H  prandtl&
       &-meyer FUNCTION NOT converged - use inlet mach no)
   374   x10=SQRT(x6-1.0)
   376   highm(j)=SQRT(1.0+x10*x10)
   x1=(highm(j)+SQRT(x6))/2.0
   IF(x5 < 1.0)x1=x1*SQRT(x5)
   IF(x1 <= 1.0)GO TO 380
   x1=x1*x1
   ws(j)=(((x4+1.0)*x1/((x4-1.0)*x1+2.0))**(x4/(x4-1.0))*((x4+1.0)/(2  &
       .0*x4*x1-x4+1.0))**(1.0/(x4-1.0))-1.0)/((1.0+(x4-1.0)/2.0*x5)**(x4  &
       /(1.0-x4))-1.0)
   380   wt(j)=wd(j)+ws(j)
 END DO
 384   IF(iprint == 1)GO TO 400
 l2=is2(i)
 l3=nterp(i)
 l4=ndata(i)
 IF(nwork(i) >= 5) CALL alg01(datac(l2),data6(l2),l4,xx6,xx5,x1,nstrms,l3,0)
 CALL alg01(datac(l2),data1(l2),l4,xx6,xx1,x1,nstrms,l3,0)
 CALL alg01(datac(l2),data4(l2),l4,xx6,xx4,x1,nstrms,l3,0)
 CALL alg01(datac(l2),data3(l2),l4,xx6,xx3,x1,nstrms,l3,0)
 ndata(i)=nstrms
 l2=l2-1
 DO  j=1,nstrms
   k=l2+j
   datac(k)=xx6(j)
   IF(nwork(i) >= 5) data6(k)=xx5(j)
   data1(k)=xx1(j)
   IF(wt(j) > wmax)wt(j)=wmax
   data2(k)=wt(j)
   data3(k)=xx3(j)
   data4(k)=xx4(j)
   data5(k)=sol(j)
 END DO
 GO TO 450
 400   IF(lnct+3 <= npage)GO TO 420
 IF(iprtc /= 0) WRITE(log2,410)
 410   FORMAT(1H1)
 lnct=4+nstrms
 420 IF(iprtc == 1) WRITE(log2,430)
 430   FORMAT(5X,   'STREAM  INLET   OUTLET  CASCADE   DIFF       LOSS &
     &DIFFUSION  BLADE  INCIDENCE  EXPANSION INLET  EXPANDED SHOCK   TOT&
     &AL',/,5X,  '-LINE   RADIUS  RADIUS  SOLIDITY  FACTOR  PARAMETER &
     &LOSS     ANGLE    ANGLE      ANGLE    M.NO  MACH NO   LOSS   LOSS ',/,2X)
 lnct=lnct+3
 DO  j=1,nstrms
   x1=beta1(j)*c1*q
   x2=xinc(j)*c1
   x3=ang(j)*c1
   IF(iprtc == 1)  &
       WRITE(log2,460)j,r(j,l1),r(j,i),sol(j),dif(j),wpara(j),wd(j),x1,x2  &
       ,x3,xmr(j),highm(j),ws(j),wt(j)
 END DO
 450   CONTINUE
 460   FORMAT(i9,f10.3,f8.3,2F9.4,f10.5,f9.5,2F9.3,f10.3,f10.4,f8.4,f8.5,  &
     f9.5)
 RETURN
END FUNCTION NOT converged - use inlet mach no)
