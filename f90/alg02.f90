SUBROUTINE alg02
     
    LOGICAL :: debug
    REAL :: loss,lami,lamip1,lamim1
    DIMENSION       ii(21,30),jj(21,30),idata(24),rdata(6),NAME(2)
    COMMON /ud3prt/ iprtc
    COMMON /udsign/ nsign
    COMMON /upage / limit,lq
    COMMON /ud300c/ nstns,nstrms,nmax,nforce,nbl,ncase,nsplit,nread,  &
        npunch,npage,nset1,nset2,istag,icase,ifailo,ipass,  &
        i,ivfail,iffail,nmix,ntrans,nplot,iloss,lnct,itub,  &
        imid,ifail,iter,log1,log2,log3,log4,log5,log6,  &
        iprint,nmany,nstplt,neqn,nspec(30),nwork(30),  &
        nloss(30),ndata(30),nterp(30),nmach(30),nl1(30),  &
        nl2(30),ndimen(30),is1(30),is2(30),is3(30),  &
        neval(30),ndiff(4),ndel(30),nliter(30),nm(2),  &
        nrad(2),ncurve(30),nwhich(30),nout1(30),nout2(30),  &
        nout3(30),nblade(30),dm(11,5,2),wfrac(11,5,2),  &
        r(21,30),xl(21,30),x(21,30),h(21,30),s(21,30),  &
        vm(21,30),vw(21,30),tbeta(21,30),diff(15,4),  &
        fdhub(15,4),fdmid(15,4),fdtip(15,4),terad(5,2),  &
        datac(100),data1(100),data2(100),data3(100),  &
        data4(100),data5(100),data6(100),data7(100),  &
        data8(100),data9(100),flow(10),speed(30),  &
        spdfac(10),bblock(30),bdist(30),wblock(30),  &
        wwbl(30),xstn(150),rstn(150),delf(30),delc(100),  &
        delta(100),title(18),drdm2(30),rim1(30),xim1(30)
    COMMON /ud300c/ work(21),loss(21),taneps(21),xi(21),vv(21),  &
        delw(21),lami(21),lamim1(21),lamip1(21),phi(21),  &
        cr(21),gama(21),sppg(21),cppg(21),hkeep(21),  &
        skeep(21),vwkeep(21),delh(30),delt(30),visk,shape,  &
        sclfac,ej,g,tolnce,xscale,pscale,plow,rlow,xmmax,  &
        rconst,fm2,hmin,c1,pi,contr,conmx
    EQUIVALENCE     (h(1,1),ii(1,1)),(s(1,1),jj(1,1))
    DATA    NAME  / 4HALG0, 4H2     /
 
    debug = .false.
    CALL sswtch (20,j)
    IF (j == 1) debug =.true.
    neval(1) = 0
    CALL fread (log1,title,18,1)
    IF (iprtc == 1) WRITE (log2,110) title
110 FORMAT (10X,10HINPUT DATA, /10X,10(1H*), //10X,5HTITLE,34X,2H= , 18A4)
    lnct = lnct + 4
    CALL alg1 (lnct)
    CALL fread (log1,idata,21,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',111,IDATA,21)
    nstns  = idata( 1)
    nstrms = idata( 2)
    nmax   = idata( 3)
    nforce = idata( 4)
    nbl    = idata( 5)
    ncase  = idata( 6)
    nsplit = idata( 7)
    nset1  = idata( 8)
    nset2  = idata( 9)
    nread  = idata(10)
    npunch = idata(11)
    nplot  = idata(12)
    npage  = idata(13)
    ntrans = idata(14)
    nmix   = idata(15)
    nmany  = idata(16)
    nstplt = idata(17)
    neqn   = idata(18)
    nle    = idata(19)
    nte    = idata(20)
    nsign  = idata(21)
    IF (nstrms == 0) nstrms = 11
    IF (nmax   == 0) nmax   = 40
    IF (nforce == 0) nforce = 10
    IF (ncase  == 0) ncase  = 1
    IF (npage  == 0) npage  = 60
    lq    = log2
    limit = npage
    CALL alg03 (lnct,19)
    IF (iprtc == 1) WRITE (log2,130) nstns,nstrms,nmax,nforce,nbl,  &
        ncase,nsplit,nset1,nset2,nread,npunch,nplot,npage,ntrans,  &
        nmix,nmany,nstplt,neqn,nle,nte,nsign
130 FORMAT (//10X,'NUMBER OF STATIONS',21X,1H=,i3, /10X,'NUMBER OF ',  &
        'STREAMLINES',18X,1H=,i3, /10X,20HMAX NUMBER of passes,19X,  &
        1H=,i3, /10X,30HMAX NUMBER of arbitrary passes,9X,1H=,i3,  &
        /10X,29HBOUNDARY layer calc indicator,10X,1H=,i3, /10X,  &
        24HNUMBER of running points,15X,1H=,i3, /10X,  &
        33HSTREAMLINE distribution indicator,6X,1H=,i3, /10X,  &
        34HNUMBER of loss/d-factor curve sets,5X,1H=,i3, /10X,  &
        34HNUMBER of loss/t.e.loss curve sets,5X,1H=,i3, /10X,  &
        26HSTREAMLINE INPUT indicator,13X,1H=,i3, /10X,  &
        27HSTREAMLINE output indicator,12X,1H=,i3, /10X,  &
        24HPRECISION plot indicator,15X,1H=,i3, /10X,  &
        24HMAX NUMBER of lines/page,15X,1H=,i3, /10X,  &
        29HWAKE transport calc indicator,10X,1H=,i3, /10X,  &
        32HMAINSTREAM mixing calc indicator,7X,1H=,i3, /10X,  &
        33HNO of stations from analytic secn,6X,1H=,i3, /10X,  &
        27HLINE-printer plot indicator,12X,1H=,i3, /10X,  &
        32HMOMENTUM equation FORM indicator,7X,1H=,i3, /10X,  &
        30HSTATION NUMBER at leading edge,9X,1H=,i3, /10X,  &
        31HSTATION NUMBER at trailing edge,8X,1H=,i3, /10X,  &
        37HCOMPRESSOR dir. of rotation indicator,2X,1H=,i3)
    itub = nstrms - 1
    imid = nstrms/2 + 1
    IF (nmany == 0) GO TO 136
    CALL fread (log1,nwhich,nmany,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',132,NWHICH,NMANY)
    CALL alg03 (lnct,2)
    IF (iprtc == 1) WRITE (log2,134) (nwhich(i),i=1,nmany)
134 FORMAT (//10X,'GEOMETRY COMES FROM ANALYTIC SECTION FOR STATIONS', 23I3)
136 CALL alg03 (lnct,7)
    CALL fread (log1,rdata,6,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',136,RDATA,6)
    g      = rdata(1)
    ej     = rdata(2)
    sclfac = rdata(3)
    tolnce = rdata(4)
    visk   = rdata(5)
    shape  = rdata(6)
    IF (g   ==   0.0) g  = 32.174
    IF (ej  ==   0.0) ej = 778.16
    IF (sclfac == 0.) sclfac = 12.0
    IF (tolnce == 0.) tolnce = 0.001
    IF (visk ==  0.0) visk  = 0.00018
    IF (shape ==  0.0) shape = 0.7
    IF (iprtc == 1) WRITE (log2,150) g,ej,sclfac,tolnce,visk,shape
150 FORMAT (//10X,22HGRAVITATIONAL constant,17X,1H=,f8.4, /10X,  &
        17HJOULES equivalent,22X,1H=,f8.3, /10X,  &
        29HLINEAR DIMENSION scale factor,10X,1H=,f8.4, /10X,  &
        15HBASIC tolerance,24X,1H=,f8.5, /10X,  &
        19HKINEMATIC viscosity,20X,1H=,f8.5, /10X, 17HB.l. shape factor,22X,1H=,f8.5)
    CALL alg03 (lnct,7)
    CALL fread (log1,rdata,6,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',151,RDATA,6)
    xscale = rdata(1)
    pscale = rdata(2)
    rlow   = rdata(3)
    plow   = rdata(4)
    xmmax  = rdata(5)
    rconst = rdata(6)
    IF (xmmax == 0.0) xmmax  = 0.6
    IF (rconst == 0.0) rconst = 6.0
    IF (iprtc == 1) WRITE (log2,160) xscale,pscale,rlow,plow,xmmax, rconst
160 FORMAT (//10X,29HPLOTTING scale for dimensions,10X,1H=,f7.3, /10X,  &
        28HPLOTTING scale for pressures,11X,1H=,f7.3, /10X,  &
        22HMINIMUM radius on plot,17X,1H=,f7.3, /10X,  &
        24HMINIMUM pressure on plot,15X,1H=,f7.3, /10X,  &
        40HMAXIMUM m-squared in relaxation factor =,f8.4, /10X,  &
        29HCONSTANT in relaxation factor,10X,1H=,f8.4)
    CALL alg03 (lnct,3)
    CALL fread (log1,rdata,2,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',162,RDATA,2)
    contr = rdata(1)
    conmx = rdata(2)
    IF (iprtc == 1) WRITE (log2,164) contr,conmx
164 FORMAT (//10X,22HWAKE transfer constant,17X,1H=,f8.5, /10X,  &
        25HTURBULENT mixing constant,14X,1H=,f8.5)
    CALL alg03 (lnct,5+ncase)
    DO  k = 1,ncase
        CALL fread (log1,flow(k),1,0)
        CALL fread (log1,spdfac(k),1,1)
    END DO
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',171,FLOW,NCASE)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',172,SPDFAC,NCASE)
    IF (iprtc == 1) WRITE(log2,180) (k,flow(k),spdfac(k),k=1,ncase)
180 FORMAT (//10X,21HPOINTS TO be computed,  //10X,2HNO,6X,8HFLOWRATE,  &
        4X,12HSPEED factor, //,(10X,i2,f13.3,f14.3))
    CALL fread (log1,l1,1,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',180,L1,1)
    DO  k = 1,l1
        CALL fread (log1,xstn(k),1,0)
        CALL fread (log1,rstn(k),1,1)
    END DO
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',191,XSTN,L1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',192,RSTN,L1)
    istag = 0
    IF (rstn(1) == 0.0) istag = 1
    nspec(1) = l1
    CALL alg03 (lnct,7+l1)
    IF (iprtc == 1) WRITE (log2,200) l1,(xstn(k),rstn(k),k=1,l1)
200 FORMAT (//10X,'ANNULUS / COMPUTING STATION GEOMETRY', //10X,  &
        24HSTATION  1  specified by,i3,7H points, //17X,4HXSTN,8X,  &
        4HRSTN,//,(f22.4,f12.4))
    is1(1) = 1
    last   = l1
    DO  i = 2,nstns
        CALL fread (log1,l1,1,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',210,L1,1)
        next = last + 1
        last = last + l1
        IF (last > 150) GO TO 550
        DO  k = next,last
            CALL fread (log1,xstn(k),1,0)
            CALL fread (log1,rstn(k),1,1)
        END DO
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',215,XSTN(NEXT),LAST-NEXT+1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',216,RSTN(NEXT),LAST-NEXT+1)
        IF (rstn(next) == 0.0) istag = i
        CALL alg03 (lnct,5+l1)
        is1(i) = next
        nspec(i) = l1
        IF (iprtc == 1) WRITE (log2,230) i,l1,(xstn(k),rstn(k), k=next,last)
    END DO
230 FORMAT (//10X,7HSTATION,i3,14H  specified by,i3,7H points, //17X,  &
        4HXSTN,8X,4HRSTN, //,(f22.4,f12.4))
    speed(1)  = 0.0
    CALL fread (log1,idata,4,1)
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',233,IDATA,4)
    l1        = idata(1)
    nterp(1)  = idata(2)
    ndimen(1) = idata(3)
    nmach(1)  = idata(4)
    DO  k = 1,l1
        CALL fread (log1,rdata,4,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',234,RDATA,4)
        datac(k) = rdata(1)
        data1(k) = rdata(2)
        data2(k) = rdata(3)
        data3(k) = rdata(4)
    END DO
    CALL alg03 (lnct,7+l1)
    is2(1)   = 1
    ndata(1) = l1
    last = l1
    IF (iprtc == 1) WRITE (log2,250) l1,nterp(1),ndimen(1),nmach(1),  &
        (datac(k),data1(k),data2(k),data3(k),k=1,l1)
250 FORMAT (//10X,24HSTATION calculation DATA,   //7X,  &
        18HSTATION  1  ndata=,i3,7H nterp=,i2,8H ndimen=,i2,  &
        7H nmach=,i2, //11X,5HDATAC,6X,14HTOTAL pressure,4X,  &
        17HTOTAL temperature,4X,11HWHIRL angle, //, (5X,f12.4,f15.4,f19.3,f18.3))
    DO  k = 1,l1
        data1(k) = data1(k)*sclfac**2
    END DO
    lastd    = 0
    nout1(1) = 0
    nout2(1) = 0
    DO  i = 2,nstns
        logn = log1
        IF (nmany == 0) GO TO 258
        DO  l1 = 1,nmany
            IF (nwhich(l1) == i) GO TO 256
        END DO
        GO TO 258
256     logn = log5
258     CALL fread (logn,idata,16,1)
        !WKBD IF (DEBUG .AND. LOGN.EQ.LOG1) CALL BUG1 ('ALG02   ',258,IDATA,16)
        ndata(i)  = idata(1)
        nterp(i)  = idata(2)
        ndimen(i) = idata(3)
        nmach(i)  = idata(4)
        nwork(i)  = idata(5)
        nloss(i)  = idata(6)
        nl1(i)    = idata(7)
        nl2(i)    = idata(8)
        neval(i)  = idata(9)
        ncurve(i) = idata(10)
        nliter(i) = idata(11)
        ndel(i)   = idata(12)
        nout1(i)  = idata(13)
        nout2(i)  = idata(14)
        nout3(i)  = idata(15)
        nblade(i) = idata(16)
        l1 = 3
        IF (ndata(i) /= 0) l1 = l1 + 5 + ndata(i)
        IF (ndel(i)  /= 0) l1 = l1 + 3 + ndel(i)
        CALL alg03 (lnct,l1)
        IF (iprtc == 1) WRITE (log2,270) i,ndata(i),nterp(i),ndimen(i),  &
            nmach(i),nwork(i),nloss(i),nl1(i),nl2(i),neval(i),ncurve(i)  &
            ,      nliter(i),ndel(i),nout1(i),nout2(i),nout3(i),nblade(i)
270     FORMAT (//7X,7HSTATION,i3, 8H  ndata=,i3,7H nterp=,i2,8H ndimen=,  &
            i2,7H nmach=,i2,7H nwork=,i2,7H nloss=,i2,5H nl1=,i3,  &
            5H nl2=,i3,7H neval=,i2,8H ncurve=,i2,8H nliter=,i3,  &
            6H ndel=,i3, /19X,6HNOUT1=,i2,7H nout2=,i2,7H nout3=,i2, 8H nblade=,i3)
        speed(i) = 0.0
        IF (ndata(i) == 0) CYCLE
        next   = last + 1
        last   = last + ndata(i)
        is2(i) = next
        IF (last > 100) GO TO 550
        CALL fread (logn,speed(i),1,1)
        !WKBD IF (DEBUG .AND.LOGN.EQ.LOG1) CALL BUG1 ('ALG02   ',271,SPEED(I),1)
        DO  k = next,last
            CALL fread (logn,rdata,6,1)
            !WKBD IF (DEBUG .AND. LOGN.EQ.LOG1) CALL BUG1 ('ALG02   ',272,RDATA,6)
            datac(k) = rdata(1)
            data1(k) = rdata(2)
            data2(k) = rdata(3)
            data3(k) = rdata(4)
            data4(k) = rdata(5)
            data5(k) = rdata(6)
            CALL fread (logn,rdata,4,1)
            !WKBD IF (DEBUG .AND. LOGN.EQ.LOG1) CALL BUG1 ('ALG02   ',273,RDATA,4)
            data6(k) = rdata(1)
            data7(k) = rdata(2)
            data8(k) = rdata(3)
            data9(k) = rdata(4)
        END DO
        IF (iprtc == 1) WRITE (log2,290) speed(i),(datac(k),data1(k),  &
            data2(k),data3(k),data4(k),data5(k),data6(k),data7(k),  &
            data8(k),data9(k),k=next,last)
290     FORMAT (//10X,7HSPEED =,f9.2, //13X,5HDATAC,7X,5HDATA1,7X,5HDATA2,  &
            7X,5HDATA3,7X,5HDATA4,7X,5HDATA5,7X,5HDATA6,7X,5HDATA7,7X,  &
            5HDATA8,7X,5HDATA9, //, (10X,f9.4,f12.3,f13.6,f11.4,f12.5,f12.5,4F12.4))
        IF (nwork(i) /= 1) GO TO 296
        DO  k = next,last
            data1(k) = data1(k)*sclfac**2
        END DO
296     IF (neval(i) > 0 .AND. nstrms > ndata(i)) last = last + nstrms -  &
            ndata(i)
        IF (ndel(i) == 0) CYCLE
        next   = lastd + 1
        lastd  = lastd + ndel(i)
        is3(i) = next
        IF (lastd > 100) GO TO 550
        DO  k = next,lastd
            CALL fread (log1,delc(k), 1,0)
            CALL fread (log1,delta(k),1,1)
        END DO
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',298,DELC(NEXT),LASTD-NEXT+1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',299,DELTA(NEXT),LASTD-NEXT+1)
        IF (iprtc == 1) WRITE(log2,310)(delc(k),delta(k),k=next,lastd)
310     FORMAT (//13X,4HDELC,8X,5HDELTA, //,(10X,f9.4,f12.4))
    END DO
    CALL alg03 (lnct,5+nstns)
    DO  i = 1,nstns
        CALL fread (log1,rdata,3,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',323,RDATA,3)
        wblock(i) = rdata(1)
        bblock(i) = rdata(2)
        bdist(i)  = rdata(3)
    END DO
    IF (iprtc == 1) WRITE (log2,340) (i,wblock(i),bblock(i), bdist(i),i=1,nstns)
340 FORMAT (//10X,'BLOCKAGE FACTOR SPECIFICATIONS', //10X,'STATION  ',  &
        ' WALL BLOCKAGE   WAKE BLOCKAGE   WAKE DISTRIBUTION FACTOR',  &
        //,(10X,i4,f16.5,f16.5,f19.3))
    IF (nset1 == 0) GO TO 380
    DO  k = 1,nset1
        CALL fread (log1,l1,1,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',342,L1,1)
        DO  j = 1,l1
            CALL fread (log1,rdata,4,1)
            !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',343,RDATA,4)
            diff(j,k)  = rdata(1)
            fdhub(j,k) = rdata(2)
            fdmid(j,k) = rdata(3)
            fdtip(j,k) = rdata(4)
        END DO
        CALL alg03 (lnct,6+l1)
        IF (iprtc == 1) WRITE (log2,360) k,l1,(diff(j,k),fdhub(j,k),  &
            fdmid(j,k),fdtip(j,k),j=1,l1)
360     FORMAT (//10X,'LOSS PARAMETER / DIFFUSION FACTOR CURVES FOR BLADE'  &
            ,      ' TYPE',i2,i5,' D-FACTORS GIVEN', //15X,9HDIFFUSION,5X,  &
            'L O S S   P A R A M E T E R S', /16X,7HFACTORS,8X,3HHUB,  &
            9X,3HMID,8X,3HTIP,//,(15X,f8.3,f13.5,f12.5,f11.5))
        ndiff(k) = l1
    END DO
380 IF (nset2 == 0) GO TO 450
    DO  k = 1,nset2
        CALL fread (log1,idata,2,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',385,IDATA,2)
        l1 = idata(1)
        l2 = idata(2)
        CALL alg03 (lnct,7+l1)
        nm(k)   = l1
        nrad(k) = l2
        CALL fread (log1,terad(1,k),1,1)
        !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',391,TERAD(1,K),1)
        DO  j = 1,l1
            CALL fread (log1,rdata,2,1)
            !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',398,RDATA,2)
            dm(j,1,k)    = rdata(1)
            wfrac(j,1,k) = rdata(2)
        END DO
        IF (iprtc == 1) WRITE (log2,410) k,l1,l2,terad(1,k),(dm(j,1,k),  &
            wfrac(j,1,k),j=1,l1)
410     FORMAT (//10X,'FRACTIONAL LOSS DISTRIBUTION CURVES FOR BLADE ',  &
            'CLASS',i2,i5,' POINTS GIVEN AT',i3,' RADIAL LOCATIONS', //  &
            10X,'FRACTION OF COMPUTING STATION LENGTH AT BLADE EXIT =',  &
            f7.4, //10X,'FRACTION OF MERIDIONAL CHORD',4X,  &
            'LOSS/LOSS AT TRAILING EDGE', //,(15X,f11.4,20X,f11.4))
        IF (l2 == 1) CYCLE
        DO  l = 2,l2
            CALL alg03 (lnct,5+l1)
            CALL fread (log1,terad(l,k),1,1)
            !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',411,TERAD(L,K),1)
            DO  j = 1,l1
                CALL fread (log1,rdata,2,1)
                !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',412,RDATA,2)
                dm(j,l,k)    = rdata(1)
                wfrac(j,l,k) = rdata(2)
            END DO
            IF (iprtc == 1) WRITE (log2,430) terad(l,k),(dm(j,l,k),  &
                wfrac(j,l,k),j=1,l1)
        END DO
430     FORMAT (//10X,'FRACTION OF COMPUTING STATION LENGTH AT BLADE ',  &
            'EXIT =',f7.4, //10X,'FRACTION OF MERIDIONAL CHORD',4X,  &
            'LOSS/LOSSAT TRAILING EDGE', //,(15X,f11.4,20X,f11.4))
    END DO
450 IF (nsplit == 0 .AND. nread == 0) GO TO 570
    DO  j = 1,nstrms,6
        CALL fread (log1,delf(j),6,1)
    END DO
    !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',455,DELF,NSTRMS)
    l1 = 5
    IF (nstrms >= 16) l1 = 8
    CALL alg03 (lnct,l1)
    IF (iprtc == 1) WRITE (log2,470)
    l1 = nstrms
    IF (nstrms > 15) l1 = 15
    IF (iprtc == 1) WRITE (log2,480) (j,j=1,l1)
480 FORMAT (//10X,'STREAMLINE',i5,14I7)
470 FORMAT (//10X,'PROPORTIONS OF TOTAL FLOW BETWEEN HUB AND EACH ',  &
        'STREAMLINE ARE TO BE AS FOLLOWS')
    IF (iprtc == 1) WRITE(log2,490) (delf(j),j=1,l1)
490 FORMAT (10X,4HFLOW,7X,15F7.4)
    IF (nstrms <= 15) GO TO 500
    l1 = l1 + 1
    IF (iprtc == 1) WRITE (log2,480) (j,j=l1,nstrms)
    IF (iprtc == 1) WRITE (log2,490) (delf(j),j=l1,nstrms)
500 IF (nread == 0) GO TO 570
    DO  i = 1,nstns
        DO  j = 1,nstrms
            CALL fread (log1,rdata,3,0)
            !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',501,RDATA,3)
            r(j,i)  = rdata(1)
            x(j,i)  = rdata(2)
            xl(j,i) = rdata(3)
            CALL fread (log1,idata,2,1)
            !WKBD IF (DEBUG) CALL BUG1 ('ALG02   ',502,IDATA,2)
            ii(j,i) = idata(1)
            jj(j,i) = idata(2)
        END DO
    END DO
    CALL alg03 (lnct,5+nstrms)
    IF (iprtc == 1) WRITE (log2,520)
520 FORMAT (//10X,'ESTIMATED STREAMLINE COORDINATES')
    DO  i = 1,nstns
        IF (i > 1) CALL alg03 (lnct,3+nstrms)
        IF (iprtc == 1) WRITE (log2,540) (i,j,r(j,i),x(j,i),xl(j,i),  &
            ii(j,i),jj(j,i),j=1,nstrms)
    END DO
540 FORMAT (//10X,'STATION  STREAMLINE   RADIUS  AXIAL COORDINATE  ',  &
        'L -COORDINATE    CHECKS-  I    J', //, (3X,2I11,f14.4,f12.4,f16.4,i17,i5))
    GO TO 570
550 WRITE  (log2,560)
560 FORMAT (////10X,'JOB STOPPED - TOO MUCH INPUT DATA')
    CALL mesage (-37,0,NAME)
570 RETURN
END SUBROUTINE alg02
