SUBROUTINE alg07
     
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
 
    l1=i+nl1(i)
    l2=i+nl2(i)
    iw=nwork(i)
    il=nloss(i)
    xn=speed(i)*spdfac(icase)*pi/(30.0*sclfac)
    SELECT CASE ( iw )
        CASE (    1)
            GO TO 100
        CASE (    2)
            GO TO 250
        CASE (    3)
            GO TO 270
        CASE (    4)
            GO TO 290
        CASE (    5)
            GO TO 440
        CASE (    6)
            GO TO 440
        CASE (    7)
            GO TO 440
    END SELECT
100 SELECT CASE ( il )
        CASE (    1)
            GO TO 110
        CASE (    2)
            GO TO 190
        CASE (    3)
            GO TO 210
        CASE (    4)
            GO TO 110
    END SELECT
110 IF(l2 /= i)GO TO 150
    DO  j=1,nstrms
        IF(ipass == 1.AND.iter == 0)GO TO 120
        IF(iter == 0)vv(j)=vm(j,i)
        x1=h(j,i)-(vv(j)**2+vw(j,i)**2)/(2.0*g*ej)
        x2=h(j,i)-(vw(j,i)**2-(vw(j,i)-xn*r(j,i))**2)/(2.0*g*ej)
        IF(x1 < hmin)x1=hmin
        IF(x2 < hmin)x2=hmin
        x3=1.0/(1.0+loss(j)*(1.0-alg4(x1,s(j,i))/alg4(x2,s(j,i))))
        GO TO 130
120     x3=1.0
130     h(j,i)=alg2(s(j,l1),work(j)/x3)
        s(j,i)=alg3(work(j),h(j,i))
    END DO
    GO TO 230
    150   DO  j=1,nstrms
        IF(ipass == 1.AND.l2 > i)GO TO 160
        x1=h(j,l1)-(vw(j,l1)**2-(vw(j,l1)-xn*r(j,l1))**2)/(2.0*g*ej)+xn**2  &
            *(r(j,i)**2-r(j,l1)**2)/(2.0*g*ej)
        IF(x1 < hmin)x1=hmin
        x2=h(j,l2)-(vm(j,l2)**2+vw(j,l2)**2)/(2.0*g*ej)
        x3=h(j,l2)-(vw(j,l2)**2-(vw(j,l2)-xn*r(j,l2))**2)/(2.0*g*ej)
        IF(x2 < hmin)x2=hmin
        IF(x3 < hmin)x3=hmin
        x4=1.0-loss(j)/alg4(x1,s(j,l1))*(alg4(x3,s(j,l2))-alg4(x2,s(j,l2)) )
        GO TO 170
160     x4=1.0
170     h(j,i)=alg2(s(j,l1),work(j)/x4)
        s(j,i)=alg3(work(j),h(j,i))
    END DO
    GO TO 230
    190   DO  j=1,nstrms
        h(j,i)=h(j,l1)+(alg2(s(j,l1),work(j))-h(j,l1))/loss(j)
        s(j,i)=alg3(work(j),h(j,i))
    END DO
    GO TO 230
    210   DO  j=1,nstrms
        s(j,i)=s(j,l1)+loss(j)
        h(j,i)=alg2(s(j,i),work(j))
    END DO
    230   DO  j=1,nstrms
        vw(j,i)=(xn*rim1(j)*vw(j,i-1)+(h(j,i)-h(j,i-1))*g*ej)/(xn*r(j,i))
    END DO
    GO TO 570
    250   DO  j=1,nstrms
        h(j,i)=work(j)
        vw(j,i)=(xn*rim1(j)*vw(j,i-1)+(h(j,i)-h(j,i-1))*g*ej)/(xn*r(j,i))
    END DO
    GO TO 330
    270   DO  j=1,nstrms
        vw(j,i)=work(j)/r(j,i)
    END DO
    GO TO 310
    290   DO  j=1,nstrms
        vw(j,i)=work(j)
    END DO
    310   DO  j=1,nstrms
        h(j,i)=h(j,i-1)+xn*(r(j,i)*vw(j,i)-rim1(j)*vw(j,i-1))/(g*ej)
    END DO
330 SELECT CASE ( il )
        CASE (    1)
            GO TO 340
        CASE (    2)
            GO TO 400
        CASE (    3)
            GO TO 420
        CASE (    4)
            GO TO 340
    END SELECT
340 IF(l2 /= i)GO TO 370
    DO  j=1,nstrms
        IF(ipass == 1.AND.iter == 0)GO TO 350
        IF(iter == 0)vv(j)=vm(j,i)
        x1=h(j,i)-(vv(j)**2+vw(j,i)**2)/(2.0*g*ej)
        x2=h(j,i)-(vw(j,i)**2-(vw(j,i)-xn*r(j,i))**2)/(2.0*g*ej)
        IF(x1 < hmin)x1=hmin
        IF(x2 < hmin)x2=hmin
        x3=1.0/(1.0+loss(j)*(1.0-alg4(x1,s(j,i))/alg4(x2,s(j,i))))
        GO TO 360
350     x3=1.0
360     s(j,i)=alg3(x3*alg4(h(j,i),s(j,l1)),h(j,i))
    END DO
    GO TO 570
    370   DO  j=1,nstrms
        IF(ipass == 1.AND.l2 > i)GO TO 380
        x1=h(j,l1)-(vw(j,l1)**2-(vw(j,l1)-xn*r(j,l1))**2)/(2.0*g*ej)+xn**2  &
            *(r(j,i)**2-r(j,l1)**2)/(2.0*g*ej)
        IF(x1 < hmin)x1=hmin
        x2=h(j,l2)-(vm(j,l2)**2+vw(j,l2)**2)/(2.0*g*ej)
        x3=h(j,l2)-(vw(j,l2)**2-(vw(j,l2)-xn*r(j,l2))**2)/(2.0*g*ej)
        IF(x2 < hmin)x2=hmin
        IF(x3 < hmin)x3=hmin
        x4=1.0-loss(j)/alg4(x1,s(j,l1))*(alg4(x3,s(j,l2))-alg4(x2,s(j,l2)) )
        GO TO 390
380     x4=1.0
390     s(j,i)=alg3(x4*alg4(h(j,i),s(j,l1)),h(j,i))
    END DO
    GO TO 570
    400   DO  j=1,nstrms
        s(j,i)=alg3(alg4(h(j,l1)+loss(j)*(h(j,i)-h(j,l1)),s(j,l1)),h(j,i))
    END DO
    GO TO 570
    420   DO  j=1,nstrms
        s(j,i)=s(j,l1)+loss(j)
    END DO
    GO TO 570
    440   DO  j=1,nstrms
        xi(j)=h(j,i-1)-xn*rim1(j)*vw(j,i-1)/(g*ej)
    END DO
    SELECT CASE ( il )
        CASE (    1)
            GO TO 460
        CASE (    2)
            GO TO 510
        CASE (    3)
            GO TO 550
        CASE (    4)
            GO TO 460
    END SELECT
460 IF(l2 /= i)GO TO 490
    DO  j=1,nstrms
        x2=xi(j)+(xn*r(j,i))**2/(2.0*g*ej)
        IF(ipass == 1.AND.iter == 0) GO TO 470
        IF(iter == 0) vv(j) = vm(j,i)
        x1=x2-vv(j)**2*(1.0+tbeta(j,i)**2)/(2.0*g*ej)
        IF(x1 < hmin)x1=hmin
        IF(x2 < hmin)x2=hmin
        x3=1.0/(1.0+loss(j)*(1.0-alg4(x1,s(j,i))/alg4(x2,s(j,i))))
        GO TO 480
470     x3=1.0
480     s(j,i)=alg3(x3*alg4(x2,s(j,l1)),x2)
    END DO
    GO TO 570
    490   DO  j=1,nstrms
        x4=xi(j)+(xn*r(j,i))**2/(2.0*g*ej)
        IF(x4 < hmin)x4=hmin
        x1=alg4(x4,s(j,l1))
        IF(ipass == 1.AND.l2 > i)GO TO 500
        x2=xi(j)+(xn*r(j,l2))**2/(2.0*g*ej)
        x3=h(j,l2)-(vm(j,l2)**2+vw(j,l2)**2)/(2.0*g*ej)
        IF(x2 < hmin)x2=hmin
        IF(x3 < hmin)x3=hmin
        x1=x1-loss(j)*(alg4(x2,s(j,l2))-alg4(x3,s(j,l2)))
500     s(j,i)=alg3(x1,x4)
    END DO
    GO TO 570
510 IF(ipass == 1.AND.iter == 0)GO TO 530
    DO  j=1,nstrms
        IF(iter == 0)vv(j)=vm(j,i)
        x1=h(j,i-1)+xn*(vv(j)*(tbeta(j,i)+xn*r(j,i)/vv(j))*r(j,i)-rim1(j)*  &
            vw(j,i-1))/(g*ej)
        IF(x1 < hmin)x1=hmin
        x2=alg4(h(j,l1)+(x1-h(j,l1))*loss(j),s(j,l1))
        s(j,i)=alg3(x2,x1)
    END DO
    GO TO 570
    530   DO  j=1,nstrms
        s(j,i)=s(j,l1)
    END DO
    GO TO 570
    550   DO  j=1,nstrms
        s(j,i)=s(j,l1)+loss(j)
    END DO

570 RETURN
END SUBROUTINE alg07
