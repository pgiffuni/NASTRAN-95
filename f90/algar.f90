SUBROUTINE algar
     
    REAL :: loss,lami,lamip1,lamim1
 
    DIMENSION xx1(21),xx2(21),xx3(21),xx4(21),vmold(21),vmlold(21)
    DIMENSION deltar(59,30),pass(59)
 
    COMMON /ud3prt/ iprtc
    COMMON /contrl/ nanal,naero,narbit,loq1,loq2,loq3,loq4,loq5,loq6
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
 
    IF = 0
    log1=loq1
    log2=loq2
    log3=loq3
    log5=loq5
    log6=loq6
    IF (iprtc == 1) WRITE(log2,1)
1   FORMAT(1HT)
    pi=3.141592653589
    c1=180.0/pi
    hmin=50.0
    vmin = 25.0
    IF (iprtc == 1) WRITE(log2,50)
50  FORMAT(1H1,37X, 53HPROGRAM alg - compressor design - aerodynamic &
        &section,/,38X,53(1H*))
    lnct=2
    CALL alg02
    icase=1
100 IF (iprtc == 1) WRITE(log2,104) icase
104 FORMAT(1H1,9X,20HOUTPUT for point no.,i2,/,10X,22(1H*))
    lnct=2
    DO  i=1,30
        DO  j=1,59
            deltar(j,i)=0.0
        END DO
    END DO
    IF((icase == 1.AND.nread == 1).OR.(icase > 1.AND.ifailk == 0))GO TO 254
    IF(nsplit == 1)GO TO 170
    l1=nspec(1)
    xx1(1)=0.0
    DO  k=2,l1
        xx1(k)=xx1(k-1)+SQRT((rstn(k)-rstn(k-1))**2+(xstn(k)-xstn(k-1))**2 )
    END DO
    x1=1.0/xx1(l1)
    DO  k=2,l1
        xx1(k)=xx1(k)*x1
    END DO
    DO  k=1,11
        xx2(k)=FLOAT(k-1)*0.1
    END DO
    CALL alg01(xx1,xstn,l1,xx2,xx3,x1,11,0,0)
    CALL alg01(xx1,rstn,l1,xx2,xx4,x1,11,0,0)
    DO  k=2,11
        xx1(k)=xx1(k-1)+SQRT((xx3(k)-xx3(k-1))**2+(xx4(k)-xx4(k-1))**2)
        xx3(k-1)=(xx1(k)+xx1(k-1))*0.5
    END DO
    l2=is1(2)
    xx2(1)=ATAN2(rstn(l2)-rstn(1),xstn(l2)-xstn(1))
    l2=l2+nspec(2)-1
    xx2(2)=ATAN2(rstn(l2)-rstn(l1),xstn(l2)-xstn(l1))
    xi(1)=0.0
    xi(2)=xx1(11)
    CALL alg01(xi,xx2,2,xx3,phi,x1,10,1,0)
    CALL alg01(rstn,xstn,l1,xx3,x1,gama,10,0,1)
    xx3(1)=0.0
    DO  k=2,11
        xx3(k)=xx3(k-1)+COS(phi(k-1)+ATAN(gama(k-1)))*(xx4(k)+xx4(k-1))*&
            (xx1(k)-xx1(k-1))
    END DO
    x1=1.0/xx3(11)
    x2=1.0/xx1(11)
    DO  k=2,11
        xx1(k)=xx1(k)*x2
        xx3(k)=xx3(k)*x1
    END DO
    x1=1.0/FLOAT(itub)
    DO  k=1,nstrms
        xx2(k)=FLOAT(k-1)*x1
    END DO
    CALL alg01(xx1,xx3,11,xx2,delf,x1,nstrms,1,0)
    170   DO  i=1,nstns
        l1=is1(i)
        l2=nspec(i)
        xx1(1)=0.0
        vv(1)=0.0
        DO  k=2,l2
            l3=l1+k-1
            vv(k)=vv(k-1)+SQRT((rstn(l3)-rstn(l3-1))**2+(xstn(l3)-xstn(l3-1))**2)
        END DO
        x1=1.0/vv(l2)
        DO  k=2,l2
            xx1(k)=vv(k)*x1
        END DO
        DO  k=1,11
            xx2(k)=FLOAT(k-1)*0.1
        END DO
        CALL alg01(xx1,xstn(l1),l2,xx2,xx3,x1,11,0,0)
        CALL alg01(xx1,rstn(l1),l2,xx2,xx4,x1,11,0,0)
        DO  k=2,11
            xx1(k)=xx1(k-1)+SQRT((xx3(k)-xx3(k-1))**2+(xx4(k)-xx4(k-1))**2)
            gama(k-1)=(xx4(k)+xx4(k-1))*0.5
            xx3(k-1)=(xx1(k)+xx1(k-1))*0.5
        END DO
        IF(i == 1.OR.i == nstns)GO TO 234
        l3=is1(i+1)
        l4=is1(i-1)
        l5=l1
        xx2(1)=(ATAN2(rstn(l3)-rstn(l5),xstn(l3)-xstn(l5))+ATAN2(rstn(l5)-  &
            rstn(l4),xstn(l5)-xstn(l4)))*0.5
        l3=l3+nspec(i+1)-1
        l4=l4+nspec(i-1)-1
        l5=l5+l2-1
        xx2(2)=(ATAN2(rstn(l3)-rstn(l5),xstn(l3)-xstn(l5))+ATAN2(rstn(l5)-  &
            rstn(l4),xstn(l5)-xstn(l4)))*0.5
        GO TO 238
234     IF(i == nstns)GO TO 236
        l3=is1(2)
        xx2(1)=ATAN2(rstn(l3)-rstn(1),xstn(l3)-xstn(1))
        l4=nspec(1)
        l3=l3+nspec(2)-1
        xx2(2)=ATAN2(rstn(l3)-rstn(l4),xstn(l3)-xstn(l4))
        GO TO 238
236     l4=is1(i-1)
        xx2(1)=ATAN2(rstn(l1)-rstn(l4),xstn(l1)-xstn(l4))
        l4=l4+nspec(i-1)-1
        l3=l1+l2-1
        xx2(2)=ATAN2(rstn(l3)-rstn(l4),xstn(l3)-xstn(l4))
238     xi(1)=0.0
        xi(2)=xx1(11)
        CALL alg01(xi,xx2,2,xx3,phi,x1,10,1,0)
        CALL alg01(rstn(l1),xstn(l1),l2,gama,x1,gama,10,0,1)
        xx3(1)=0.0
        DO  k=2,11
            xx3(k)=xx3(k-1)+COS(phi(k-1)+ATAN(gama(k-1)))*(xx4(k)+xx4(k-1))*&
                (xx1(k)-xx1(k-1))
        END DO
        x1=1.0/xx3(11)
        DO  k=2,11
            xx3(k)=xx3(k)*x1
        END DO
        CALL alg01(xx3,xx1,11,delf,xl(1,i),x1,nstrms,1,0)
        x1=vv(l2)/xx1(11)
        DO  j=2,nstrms
            xl(j,i)=xl(j,i)*x1
        END DO
        CALL alg01(vv,xstn(l1),l2,xl(1,i),x(1,i),x1,nstrms,0,0)
        CALL alg01(vv,rstn(l1),l2,xl(1,i),r(1,i),x1,nstrms,0,0)
    END DO
254 IF(icase > 1)GO TO 270
    x1=(x(imid,2)-x(imid,1))**2+(r(imid,2)-r(imid,1))**2
    drdm2(1)=((r(nstrms,1)-r(1,1))**2+(x(nstrms,1)-x(1,1))**2)/x1
    l1=nstns-1
    DO  i=2,l1
        x2=(x(imid,i+1)-x(imid,i))**2+(r(imid,i+1)-r(imid,i))**2
        x3=x2
        IF(x1 < x3)x3=x1
        drdm2(i)=((r(nstrms,i)-r(1,i))**2+(x(nstrms,i)-x(1,i))**2)/x3
        x1=x2
    END DO
    drdm2(nstns)=((r(nstrms,nstns)-r(1,nstns))**2+(x(nstrms,nstns)-&
        x(1,nstns))**2)/x2
    270   DO  i=1,nstns
        wwbl(i)=wblock(i)
    END DO
    ipass=1
290 i=1
    IF((ipass > 1.OR.icase > 1).AND.ndata(1) == 1)GO TO 400
    l1=ndimen(1)+1
    SELECT CASE ( l1 )
        CASE (    1)
            GO TO 300
        CASE (    2)
            GO TO 320
        CASE (    3)
            GO TO 340
        CASE (    4)
            GO TO 360
    END SELECT
    300   DO  j=1,nstrms
        xx1(j)=r(j,1)
    END DO
    GO TO 380
    320   DO  j=1,nstrms
        xx1(j)=r(j,1)/r(nstrms,1)
    END DO
    GO TO 380
    340   DO  j=1,nstrms
        xx1(j)=xl(j,1)
    END DO
    GO TO 380
    360   DO  j=1,nstrms
        xx1(j)=xl(j,1)/xl(nstrms,1)
    END DO
380 l1=nterp(1)
    l2=ndata(1)
    CALL alg01(datac,data1,l2,xx1,s    ,x1,nstrms,l1,0)
    CALL alg01(datac,data2,l2,xx1,h    ,x1,nstrms,l1,0)
    CALL alg01(datac,data3,l2,xx1,tbeta,x1,nstrms,l1,0)
    DO  j=1,nstrms
        h(j,1)=alg6(s(j,1),h(j,1))
        s(j,1)=alg3(s(j,1),h(j,1))
        tbeta(j,1)=TAN(tbeta(j,1)/c1)
    END DO
400 IF(ipass > 1.OR.icase > 1)GO TO 420
    x1=flow(1)/(alg5(h,s)*pi*(r(nstrms,1)+r(1,1))*xl(nstrms,1))*sclfac **2
    DO  j=1,nstrms
        vm(j,1)=x1
    END DO
    IF(istag == 1)vm(1,1)=0.0
420 ifailo=0
    iffail=0
    ivfail=0
    DO  j=1,nstrms
        vmold(j)=vm(j,1)
    END DO
    GO TO 500
440 IF(ipass > 1)GO TO 460
    DO  j=1,nstrms
        vm(j,i)=vm(j,i-1)
    END DO
    IF(i-1 == istag)vm(1,i)=vm(2,i)
    IF(i == istag)vm(1,i)=0.0
460 iloss=1
    DO  j=1,nstrms
        vmold(j)=vm(j,i)
    END DO
    470   DO  j=1,nstrms
        vwkeep(j)=vw(j,i-1)
        skeep(j)=s(j,i-1)
        hkeep(j)=h(j,i-1)
    END DO
    x1=h(imid,i-1)-(vm(imid,i-1)**2+vw(imid,i-1)**2)/(2.0*g*ej)
    IF(x1 < hmin)x1=hmin
    psmid=alg4(x1,s(imid,i-1))
    IF(nmix == 1)CALL alg04(h(1,i-1),s(1,i-1),vw(1,i-1),r(1,i-1),r(1,  &
        i),x(1,i-1),x(1,i),vm(1,i-1),conmx,sclfac,g,ej,hmin,vmin,psmid,&
        nstrms,log2,lnct,IF)
    IF(IF == 0)GO TO 478
    ifailo=i-1
    GO TO 640
478 IF(nwork(i) == 0)GO TO 480
    CALL alg05
    IF(ntrans == 1.AND.ipass > 1)CALL alg06(r(1,i-1),r(1,i),x(1,i-1)  &
        ,x(1,i),h(1,i),s(1,i),vm(1,i),tbeta(1,i-1),tbeta(1,i),loss,contr,&
        sclfac,speed(i),spdfac(icase),g,ej,hmin,nstrms,pi)
    iter=0
    CALL alg07
    GO TO 500
    480   DO  j=1,nstrms
        h(j,i)=h(j,i-1)
        s(j,i)=s(j,i-1)
        vw(j,i)=0.0
        IF(i > istag.OR.j /= 1)vw(j,i)=vw(j,i-1)*rim1(j)/r(j,i)
    END DO
    500   DO  j=1,nstrms
        vmlold(j)=vm(j,i)
    END DO
    IF(neqn >= 2)GO TO 514
    CALL alg08
    GO TO 516
514 CALL alg26
516 IF(neval(i) <= 0)GO TO 590
    iprint=0
    CALL alg09
    IF(ifailo /= 0.AND.ipass > nforce)GO TO 550
    DO  j=1,nstrms
        IF(ABS(vm(j,i)/vmlold(j)-1.0) > tolnce/5.0)GO TO 530
    END DO
    GO TO 590
530 IF(iloss >= nliter(i))GO TO 550
    iloss=iloss+1
    DO  j=1,nstrms
        vmlold(j)=vm(j,i)
    END DO
    GO TO 470
550 IF(ipass <= nforce)GO TO 590
    IF(lnct+1 <= npage)GO TO 570
    IF (iprtc == 1) WRITE(log2,560)
560 FORMAT(1H1)
    lnct=1
570 lnct=lnct+1
    x1=vm(1,i)/vmlold(1)
    x2=vm(imid,i)/vmlold(imid)
    x3=vm(nstrms,i)/vmlold(nstrms)
    IF (iprtc == 1) WRITE(log2,580) ipass,i,x1,x2,x3
580 FORMAT(5X,4HPASS,i3,9H  station,i3,66H  vm profile NOT converged w&
        &ith loss recalc   vm NEW/vm prev  hub=,f9.6,6H  mid=,f9.6,7H case =,f9.6)
590 IF(nbl == 1.AND.(ifailo == 0.OR.ipass <= nforce))CALL alg10
    DO  j=1,nstrms
        xim1(j)=x(j,i)
        rim1(j)=r(j,i)
        IF(i == istag.AND.j == 1)CYCLE
        IF(ABS(vm(j,i)/vmold(j)-1.0) > tolnce)ivfail=ivfail+1
        IF(ABS(delw(j)-delf(j)) > tolnce)iffail=iffail+1
    END DO
    IF(nmax == 1.OR.(ipass == 1.AND.nread == 1))GO TO 624
    x1=fm2
    IF(x1 < 1.0-xmmax)x1=1.0-xmmax
    x2=1.0
    IF(i == 1.OR.nwork(i) >= 5)x2=1.0+tbeta(imid,i)**2
    x1=1.0/(1.0+x1*drdm2(i)/(rconst*x2))
    l3=nstrms-2
    CALL alg01(delw,xl(1,i),nstrms,delf(2),xx1(2),x1,l3,1,0)
    xx=xl(imid,i)
    DO  j=2,itub
        xl(j,i)=xl(j,i)+x1*(xx1(j)-xl(j,i))
    END DO
    l1=ipass
    IF(l1 <= 59)GO TO 618
    l1=59
    DO  k=1,58
        deltar(k,i)=deltar(k+1,i)
    END DO
618 deltar(l1,i)=xl(imid,i)-xx
    l1=is1(i)
    l2=nspec(i)
    xx1(1)=0.0
    DO  k=2,l2
        kk=l1-1+k
        xx1(k)=xx1(k-1)+SQRT((xstn(kk)-xstn(kk-1))**2+(rstn(kk)-rstn(kk-1) )**2)
    END DO
    CALL alg01(xx1,rstn(l1),l2,xl(2,i),r(2,i),x1,l3,0,0)
    CALL alg01(xx1,xstn(l1),l2,xl(2,i),x(2,i),x1,l3,0,0)
624 IF(ipass > nforce.AND.ifailo /= 0)GO TO 640
    IF(i == nstns)GO TO 630
    i=i+1
    GO TO 440
630 IF(ipass >= nmax)GO TO 640
    IF(ifailo /= 0)GO TO 635
    IF(ivfail == 0.AND.iffail == 0)GO TO 640
635 ipass=ipass+1
    GO TO 290
640 CALL alg11
    l1=nstns
    IF(ifailo /= 0)l1=ifailo
    iprint=1
    DO  i=2,l1
        IF(neval(i) /= 0)CALL alg09
    END DO
    IF(nplot /= 0)CALL alg12
    IF(ifailo /= 0)GO TO 750
    IF(npunch == 0)GO TO 680
    WRITE(log3,660)(delf(j),j=1,nstrms)
660 FORMAT(6F12.8)
    WRITE(log3,670)((r(j,i),x(j,i),xl(j,i),i,j,j=1,nstrms),i=1,nstns)
670 FORMAT(3F12.8,2I3)
    680   DO  i=1,nstns
        IF(nout1(i) == 0)CYCLE
        WRITE(log3,690)(r(j,i),j,i,j=1,nstrms)
690     FORMAT(f12.8,60X,2I4)
    END DO
    l1=log3
    IF(narbit /= 0)l1=log6
    DO  i=1,nstns
        IF(nout2(i) == 0)CYCLE
        l2=is1(i)
        l3=l2+nspec(i)-1
        WRITE(l1,710)nspec(i),(xstn(k),rstn(k),k=l2,l3)
710     FORMAT(i3,/,(2F12.7))
        xn=speed(i)
        IF(i == nstns)GO TO 714
        IF(speed(i) /= speed(i+1).AND.nwork(i+1) /= 0)xn=speed(i+1)
714     xn=xn*spdfac(icase)*pi/(30.0*sclfac)
        DO  j=1,nstrms
            xx1(j)=ATAN((vw(j,i)-xn*r(j,i))/vm(j,i))*c1
        END DO
        WRITE(l1,730)(r(j,i),xx1(j),j,i,j=1,nstrms)
730     FORMAT(2F12.8,48X,2I4)
    END DO
750 IF(nstplt == 0)GO TO 759
    l1=ipass
    IF(l1 > 59)l1=59
    DO  k=1,l1
        pass(k)=FLOAT(k)
    END DO
    DO  k=1,nstns
        IF (iprtc == 1) WRITE(log2,756) k
756     FORMAT(1H1,53X,19HDELTA l for station,i3,/,2X)
        CALL alg25(l1,ipass,log2,pass,deltar(1,k))
    END DO
759 IF(icase >= ncase)GO TO 760
    icase=icase+1
    ifailk=ifailo
    GO TO 100
760 IF (iprtc == 1) WRITE(log2,770)
770 FORMAT(1HS)
 
    RETURN
END SUBROUTINE algar
