SUBROUTINE alg05
     
 REAL :: loss,lami,lamip1,lamim1
 
 DIMENSION xx1(21),xx2(21),xx3(21),xx4(21),xx5(21)
 
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
 
 l1=ndimen(i)+1
 SELECT CASE ( l1 )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 140
   CASE (    4)
     GO TO 160
 END SELECT
 100   DO  j=1,nstrms
   xx5(j)=r(j,i)
 END DO
 GO TO 180
 120   DO  j=1,nstrms
   xx5(j)=r(j,i)/r(nstrms,i)
 END DO
 GO TO 180
 140   DO  j=1,nstrms
   xx5(j)=xl(j,i)
 END DO
 GO TO 180
 160   DO  j=1,nstrms
   xx5(j)=xl(j,i)/xl(nstrms,i)
 END DO
 180   l2=is2(i)
 l3=ndata(i)
 l4=nterp(i)
 CALL alg01(datac(l2),data1(l2),l3,xx5,work  ,x1,nstrms,l4,0)
 CALL alg01(datac(l2),data3(l2),l3,xx5,taneps,x1,nstrms,l4,0)
 DO  j=1,nstrms
   taneps(j)=TAN(taneps(j)/c1)
 END DO
 iw=nwork(i)
 il=nloss(i)
 IF(iw == 7.OR.il <= 3)  &
     CALL alg01(datac(l2),data2(l2),l3,xx5,loss  ,x1,nstrms,l4,0)
 IF(iw >= 5) CALL alg01(datac(l2),data6(l2),l3,xx5,xx1,x1,nstrms,l4,0)
 IF(il /= 4)GO TO 350
 DO  ii=i,nstns
   IF(nloss(ii) == 1)EXIT
 END DO
 210   l2=is2(ii)
 l3=ndata(ii)
 l4=nterp(ii)
 l1=ndimen(ii)+1
 SELECT CASE ( l1 )
   CASE (    1)
     GO TO 220
   CASE (    2)
     GO TO 240
   CASE (    3)
     GO TO 260
   CASE (    4)
     GO TO 280
 END SELECT
 220   DO  j=1,nstrms
   xx5(j)=r(j,ii)
 END DO
 GO TO 300
 240   DO  j=1,nstrms
   xx5(j)=r(j,ii)/r(nstrms,ii)
 END DO
 GO TO 300
 260   DO  j=1,nstrms
   xx5(j)=xl(j,ii)
 END DO
 GO TO 300
 280   DO  j=1,nstrms
   xx5(j)=xl(j,ii)/xl(nstrms,ii)
 END DO
 300   CALL alg01(datac(l2),data2(l2),l3,xx5,loss,x1,nstrms,l4,0)
 iii=i+nl1(i)+1
 DO  j=1,nstrms
   xx2(j)=0.0
   DO  ik=iii,ii
     xx2(j)=xx2(j)+SQRT((x(j,ik)-x(j,ik-1))**2+(r(j,ik)-r(j,ik-1))**2)
     IF(ik == i)xx3(j)=xx2(j)
   END DO
   xx3(j)=xx3(j)/xx2(j)
 END DO
 l1=ncurve(i)
 l2=nm(l1)
 l3=nrad(l1)
 DO  j=1,nstrms
   DO  k=1,l3
     CALL alg01(dm(1,k,l1),wfrac(1,k,l1),l2,xx3(j),xx2(k),x1,1,0,0)
   END DO
   x2=(r(j,ii)-r(1,ii))/(r(nstrms,ii)-r(1,ii))
   CALL alg01(terad(1,l1),xx2,l3,x2,x1,x1,1,0,0)
   loss(j)=loss(j)*x1
 END DO
 350   IF(iw < 5)GO TO 420
 IF(iw /= 5)GO TO 370
 DO  j=1,nstrms
   tbeta(j,i)=TAN((work(j)+xx1(j))/c1)
 END DO
 GO TO 420
 370   IF(iw == 7)GO TO 400
 DO  j=1,nstrms
   xx2(j)=TAN((ATAN((r(j,i+1)-r(j,i))/(x(j,i+1)-x(j,i)))+ATAN((r(j,i)  &
       -r(j,i-1))/(x(j,i)-x(j,i-1))))/2.0)
 END DO
 l1=is1(i)
 CALL alg01(rstn(l1),xstn(l1),nspec(i),r(1,i),x1,xx3,nstrms,0,1)
 DO  j=1,nstrms
   tbeta(j,i)=TAN(ATAN((TAN(work(j)/c1)*(1.0-xx3(j)*xx2(j))-xx2(j)*ta  &
       neps(j)*SQRT(1.0+xx3(j)**2))/SQRT(1.0+xx2(j)**2))+xx1(j)/c1)
 END DO
 GO TO 420
 400   xn=speed(i)*spdfac(icase)*pi/(30.0*sclfac)
 CALL alg01(datac(l2),data7(l2),l3,xx5,xx2,x1,nstrms,l4,0)
 CALL alg01(datac(l2),data8(l2),l3,xx5,xx3,x1,nstrms,l4,0)
 CALL alg01(datac(l2),data9(l2),l3,xx5,xx4,x1,nstrms,l4,0)
 ii=i+nl1(i)
 DO  j=1,nstrms
   x1=c1*ATAN((vw(j,ii)-xn*r(j,ii))/vm(j,ii))
   x2=xx3(j)
   IF(x1 < xx1(j))x2=xx4(j)
   loss(j)=loss(j)*(1.0+((x1-xx1(j))/(x2-xx1(j)))**2)
   IF(loss(j) > 0.5)loss(j)=0.5
   tbeta(j,i)=TAN((work(j)+(x1-xx1(j))*xx2(j))/c1)
 END DO
 420   RETURN
END SUBROUTINE alg05
