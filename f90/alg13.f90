SUBROUTINE alg13 (ibl,ys,yp,xs,xp,ysemi,xsemi,log1,log2,n,iprint,  &
                  beta1,beta2,p,q,yzero,t,yone,xdel,ydel,z,axialc, &
                  lnct,ifcord,sq,sb,isecn,xsemj,ysemj,istak,xhere, &
                  x,ss,nstns,r,dx,y,dy,ss1,bx,sigma,ccord,isplit,  &
                  yzeros,ts,yones,zspmxt,perspj,inast,irle,irte,   &
                  tharr)
     
 
    INTEGER, INTENT(IN OUT)                  :: ibl
    REAL, INTENT(OUT)                        :: ys(21,80)
    REAL, INTENT(OUT)                        :: yp(21,80)
    REAL, INTENT(OUT)                        :: xs(21,80)
    REAL, INTENT(OUT)                        :: xp(21,80)
    REAL, INTENT(OUT)                        :: ysemi(21,31)
    REAL, INTENT(OUT)                        :: xsemi(21,31)
    INTEGER, INTENT(IN OUT)                  :: log1
    INTEGER, INTENT(IN OUT)                  :: log2
    INTEGER, INTENT(OUT)                     :: n
    INTEGER, INTENT(IN OUT)                  :: iprint
    REAL, INTENT(IN OUT)                     :: beta1
    REAL, INTENT(IN OUT)                     :: beta2
    REAL, INTENT(OUT)                        :: p
    REAL, INTENT(OUT)                        :: q
    REAL, INTENT(OUT)                        :: yzero
    REAL, INTENT(OUT)                        :: t
    REAL, INTENT(OUT)                        :: yone
    REAL, INTENT(IN)                         :: xdel
    REAL, INTENT(IN)                         :: ydel
    REAL, INTENT(IN)                         :: z
    REAL, INTENT(OUT)                        :: axialc
    INTEGER, INTENT(OUT)                     :: lnct
    INTEGER, INTENT(IN OUT)                  :: ifcord
    REAL, INTENT(IN OUT)                     :: sq
    REAL, INTENT(OUT)                        :: sb
    INTEGER, INTENT(IN OUT)                  :: isecn
    REAL, INTENT(OUT)                        :: xsemj(21,31)
    REAL, INTENT(OUT)                        :: ysemj(21,31)
    INTEGER, INTENT(IN)                      :: istak
    REAL, INTENT(IN OUT)                     :: xhere(100)
    REAL, INTENT(OUT)                        :: x(100)
    REAL, INTENT(OUT)                        :: ss(100)
    INTEGER, INTENT(IN OUT)                  :: nstns
    REAL, INTENT(IN OUT)                     :: r(10,21)
    REAL, INTENT(OUT)                        :: dx(100)
    REAL, INTENT(OUT)                        :: y(100)
    REAL, INTENT(OUT)                        :: dy(100)
    REAL, INTENT(IN OUT)                     :: ss1(80,4)
    REAL, INTENT(IN)                         :: bx
    REAL, INTENT(OUT)                        :: sigma(100)
    REAL, INTENT(OUT)                        :: ccord(1)
    REAL, INTENT(IN OUT)                     :: yzeros
    REAL, INTENT(OUT)                        :: ts
    REAL, INTENT(OUT)                        :: yones
    REAL, INTENT(IN)                         :: zspmxt
    REAL, INTENT(IN)                         :: perspj

    INTEGER, INTENT(IN)                      :: isplit
    INTEGER, INTENT(IN OUT)                  :: inast
    INTEGER, INTENT(IN)                      :: irle
    INTEGER, INTENT(IN)                      :: irte

    REAL, INTENT(OUT)                        :: tharr(21,10)
    REAL :: ix,iy,ixy,ipx,ipy,ixd,iyd,ixn,iyn,ixyn

    DIMENSION  NAME(2),xxm(81),snadum(10)
    DIMENSION  s(80),phi(11), thick2(80),xm(81),ym(80),am(80)
    DIMENSION  xspltm(45),yspltm(45),sspltm(45),xsplts(45),  &
               ysplts(45),xspltp(45),yspltp(45),thick(45)

    COMMON /udstr2/ nbldes,stag(21),chordd(21)

    DATA    NAME  / 4HALG1, 4H3    /
 
    f1(a) = a*EXP(1.0-a*sq)*sq
    f2(a) = (sq-1.0)*a*EXP(1.0+a*(1.0-sq))
    f3(a,b,c,d) = b/a**3*EXP(a*xd)*(a*xd-2.0) + c*(xd+sq) + d
    f4(a,b) = ABS(a-b)/(a-b)
    f5(a,b,c) = b/a**2*EXP(a*xd)*(a*xd-1.0) + c
    f6(xab) = SQRT(rdius**2-(xab-x1)**2) + y1
    f7(xab) =-SQRT(rdius**2-(xab-x1)**2) + y1
    f8(xab) =-1./SQRT(rdius**2-(xab-x1)**2)*(xab-x1)
 
10  FORMAT (1H1)
    a    = 0.
    d    = 0.
    bta1 = beta1
    bta2 = beta2
    beta3= 0.0
    pi   = 3.1415926535
    c1   = 180.0/pi
    IF (iprint >= 2) GO TO 40
 
    WRITE  (log2,20) ibl,p,q,beta1,beta2,yzero,t,yone,z,axialc
20  FORMAT (1H1,44X,'STREAMSURFACE geometry on streamline NUMBER',i3, &
        /45X,46(1H*), //20X,'P',15X,'=',f7.4,6X,&
        '(d2ydx2 of meanline at leading edge as a fraction of its maximum value.)',&
        /20X,1HQ,15X,1H=,f7.4,6X,&
        '(d2ydx2 of meanline at trailing edge as a fraction of its maximum value.)',&
        /20x,'BETA1',11X,'=',f7.3,6X,'(blade inlet angle.)', /20X,'BETA2',11X,&
        '=',f7.3,6X,'(blade outlet angle.)', /20X,'YZERO',11X,'=',f8.5,  &
        5X,'(blade leading edge radius as a fraction of chord.)', /20X,'t',15X,'='&
        ,f8.5,5X,'(blade maximum thickness as a fraction of chord.)', /20X,&
        'YONE',12X,'=',f8.5,5X,'(blade trailing edge half-thickness as a fraction of chord.)',&
        /20X,'Z',15X,'=',f7.4,6X,'(location of maximum thickness as a fraction of mean line.)',&
        /20X,'CORD',12X,1H=,f7.4,6X,'(chord OR meridional chord of section.)')
    IF (isecn == 1 .OR. isecn == 3) WRITE (log2,30) sq,sb
30  FORMAT (20X,1HS,15X,1H=,f7.4,6X,'(inflection point as a fractionof meridional chord.)',&
        /20X,'BETA3',11X,'=',f7.3,6X,'(change inangle from leading edge.)')
40  IF (iprint == 3) GO TO 55
    lnct = lnct + 2
    IF (lnct <= 60)GO TO 45
    lnct = 3
    WRITE  (log2,10)
45  WRITE  (log2,50) ibl,p,q,beta1,beta2,yzero,t,yone,z,axialc
50  FORMAT (2X, /5X,'LINE',i3,'  p=',f7.4,'  q=',f7.4,'  beta1=',f7.3,  &
        '  beta2=',f7.3,'  yzero=',f7.5,'  t/c=',f7.5,'  yone=',f7.5,'  z=',&
        f7.4,6H  axc=,f7.3)
55  IF (isecn == 1) GO TO 60
    IF (isecn == 3) GO TO 130
    IF (isecn == 2) GO TO 150
    h  = 1.0/(1.0+SQRT((1.0-q)/(1.0-p)))
    hh = h*h
    oa = 4.0*(TAN(beta1/c1)-TAN(beta2/c1))/(p/(1.0-p)*hh+h-1.0/3.0)
    oa48 = oa/48.0
    xk2  =-hh/(8.0*(1.0-p))*oa
    b    = hh*h/12.0*oa+TAN(beta1/c1)
    c    =-hh*hh*oa48
    xmlc = SQRT(1.0+(oa48*(1.0-h)**4+xk2+b+c)**2)
    GO TO 160
60  nq = 1
    sb = beta1+sb
    g1 = 1.0/sq
    r1 = f1(g1)
    g2 = g1+5.0
    r2 = f1(g2)
    s2 = f4(r2,p)
70  g3 = g2+(p-r2)*(g2-g1)/(r2-r1)
    r3 = f1(g3)
    s3 = f4(r3,p)
    IF (ABS(r3-p) <= 0.0001) GO TO 90
    IF (nq > 50) GO TO 1290
    nq = nq + 1
    IF (ABS(s2-s3) <= 0.0001) GO TO 80
    g1 = g3
    r1 = r3
    GO TO 70
80  g2 = g3
    r2 = r3
    s2 = s3
    GO TO 70
90  a1 = g3
    nq = 1
    g1 = 1.0/(sq-1.0)
    r1 = f2(g1)
    g2 = g1-5.0
    r2 = f2(g2)
    s2 = f4(r2,q)
100 g3 = g2+(q-r2)*(g2-g1)/(r2-r1)
    r3 = f2(g3)
    s3 = f4(r3,q)
    IF (ABS(r3-q) <= 0.0001) GO TO 120
    IF (nq > 50) GO TO 1290
    nq = nq + 1
    IF (ABS(s2-s3) <= 0.0001) GO TO 110
    g1 = g3
    r1 = r3
    GO TO 100
110 g2 = g3
    r2 = r3
    s2 = s3
    GO TO 100
120 a2 = g3
    b1 = a1**2*(TAN(beta1/c1)-TAN(sb/c1))/ (1.0-(a1*sq+1.0)*EXP(-a1*sq))
    cc1= TAN(sb/c1)+b1/a1**2
    e1 = (a1*sq+2.0)*b1/a1**3*EXP(-a1*sq)
    b2 = a2**2*(TAN(beta2/c1)-TAN(sb/c1))/  &
        (1.0+(a2*(1.0-sq)-1.0)*EXP(a2*(1.0-sq)))
    cc2= TAN(sb/c1)+b2/a2**2
    d2 = 2.*(b2/a2**3-b1/a1**3)+sq*(cc1-cc2) + e1
    xd = 1.0-sq
    r2 = f3(a2,b2,cc2,d2)
    xmlc = SQRT(1.0+r2**2)
    GO TO 160
130 i1 = 1
    beta3 = beta1+sb
    s0 = 0.
    x0 = 0.
    y0 = 0.
    y21= 0.0
    i2 = FLOAT(n)*sq
    IF (i2 <= 1) sq = 0.0
    IF (i2 <= 1) beta3 = beta1
    IF (i2 <= 1) GO TO 140
    xrnge = sq
    fact  = sq
    CALL alg18 (beta1,beta3,i1,i2,fact,x0,y0,s0,xrnge,y11,x11,y21, rdius1,s,c1)
    i1 = i2
    x0 = sq
    y0 = y21
    s0 = s(i1)
140 i2 = n
    fact  = 1.-sq
    xrnge = fact
    CALL alg18 (beta3,beta2,i1,i2,fact,x0,y0,s0,xrnge,y12,x12,y22, rdius2,s,c1)
    xmlc = SQRT(1.0+y22**2)
    GO TO 160
150 CALL alg18 (beta1,beta2,1,n,1.0,0.0,0.0,0.0,1.0,y1,x1,y2,rdius, s,c1)
    xmlc  = SQRT(1.+y2**2)
    chord = xmlc/(1.-2.*yzero*(1.-xmlc))
    fcslmn= 1.0 - chord*2.*yzero
    GO TO 170
160 chord = xmlc/(1.0-yzero+xmlc*(yzero+ABS(yone*SIN(beta2/c1))))
    fcslmn= 1.0 - chord*(yzero+ABS(yone*SIN(beta2/c1)))
170 IF (ifcord == 1) axialc = axialc/chord
    yzero = yzero*chord/fcslmn
    yone  = yone*chord/fcslmn
    t   = t*chord/fcslmn
    s(1)= 0.0
    xx  = 0.0
    xn  = n
    IF (isecn == 2) GO TO 240
    at = (yzero-t/2.0)/(2.0*z**3)
    ct = (t/2.0-yzero)*3.0/(2.0*z)
    dt = yzero
    et = (yone-t/2.0)/(1.0-z)**3 - 1.5*(yzero-t/2.0)/(z**2*(1.0-z))
    ft = 1.5*(yzero-t/2.0)/z**2
    ht = t/2.0
    IF (isecn == 3) GO TO 240
    delx = 1.0/(10.0*(xn-1.0))
    ASSIGN 190 TO isec1
    ASSIGN 290 TO isec2
    IF (isecn == 0) GO TO 180
    ASSIGN 200 TO isec1
    ASSIGN 300 TO isec2
    180   DO  j  = 2,n
        DO  jj = 1,11
            GO TO isec1, (190,200)
190         phi(jj) = SQRT(1.0+(oa/12.0*(xx-h)**3+xk2*2.0*xx+b)**2)
            GO TO 220
200         xd = xx - sq
            IF (xd > 0.0) GO TO 210
            phi(jj) = SQRT(1.0+(f5(a1,b1,cc1))**2)
            GO TO 220
210         phi(jj) = SQRT(1.+(f5(a2,b2,cc2))**2)
220         xx = xx + delx
        END DO
        xx = xx - delx
        s(j) = s(j-1)+ (phi(1)+ phi(11)+ 4.0*(phi(2)+phi(4)+phi(6)+phi(8)+  &
            phi(10))+2.0*(phi(3)+phi(5)+phi(7)+phi(9)))/(30.0*(xn-1.0))
    END DO
240 delx = 1.0/(xn-1.0)
    IF (isecn /= 2) GO TO 250
    t2 = t/2.
    tprim2 = t2 - yzero
    c2 = 2.*c1
    aform = (tprim2+rdius*(1.-COS((beta1-beta2)/c2)))/xmlc*2.
    phis  = ACOS((1.-aform**2)/(1.+aform**2))
    rs    = yzero + xmlc/2./SIN(phis)
    yss   = rdius - rs + t2
    bform = (rdius*(1.-COS((beta1-beta2)/c2))-tprim2)/xmlc*2.
    phip  = ACOS((1.-bform**2)/(1.+bform**2))
    phi2  = ABS((beta1-beta2)/c1)
250 xm(1) = 0.0
    IF (isecn /= 3) GO TO 260
    ymm = 0.0
    xmm = 0.0
    i2  = sq*FLOAT(n)
    i3  = i2
    IF (i2 <=  1) i2 = n + 1
    delx= sq/FLOAT(i2-1)
    IF (i3 /= i2) i3 = 1
    delxx = (1.-sq)/FLOAT(n-i3)
    IF (i2 == n+1) delx = delxx
    260   DO  j = 1,n
        sn = s(j)/s(n)
        IF (isecn == 2) GO TO 340
        IF (sn    > z) GO TO 270
        thick2(j) = (at*sn**2+ct)*sn + dt
        GO TO 280
270     sn = sn - z
        thick2(j) = (et*sn+ft)*sn**2 + ht
280     IF (isecn == 3) GO TO 320
        GO TO isec2, (290,300)
290     ym(j)  = oa48*(xm(j)-h)**4+xk2*xm(j)**2 + b*xm(j) + c
        yprime = oa/12.0*(xm(j)-h)**3 + xk2*2.0*xm(j) + b
        GO TO 370
300     xd = xm(j) - sq
        IF (xd > 0.0) GO TO 310
        ym(j)  = f3(a1,b1,cc1,e1)
        yprime = f5(a1,b1,cc1)
        GO TO 370
310     ym(j)  = f3(a2,b2,cc2,d2)
        yprime = f5(a2,b2,cc2)
        GO TO 370
320     IF (xm(j)-sq > 0.0 .OR. xm(j) == 0.0 .AND. sq == 0.0) GO TO 330
        IF (beta1 == beta3) GO TO 360
        bta1  = beta1
        bta2  = beta3
        rdius = rdius1
        y1 = y11
        x1 = x11
        GO TO 350
330     IF (beta2 == beta3) GO TO 360
        rdius = rdius2
        x1 = x12
        y1 = y12
        bta1 = beta3
        bta2 = beta2
        GO TO 350
340     phix = (sn-0.5)*phi2
        thick2(j) = yss*COS(phix) + SQRT(rs**2-yss**2*SIN(phix)**2) -rdius
350     ym(j)  = f6(xm(j))
        yprime = f8(xm(j))
        IF (bta1-bta2 < 0.0) yprime = -yprime
        IF (bta1-bta2 < 0.0) ym(j) = f7(xm(j))
        IF (isecn == 2) GO TO 370
        IF (j == i3) delx = delxx
        GO TO 370
360     yprime = TAN(beta3/c1)
        IF (j /= 1) xmm = xm(j-1)/fcslmn - yzero
        IF (j /= 1) ymm = ym(j-1)/fcslmn
        ym(j) = yprime*(xm(j)-xmm) + ymm
        IF (j == i3) delx = delxx
370     xm(j+1) = xm(j) + delx
        fypr = 1.0/SQRT(1.0+yprime**2)
        xs(ibl,j) = (xm(j)-thick2(j)*yprime*fypr+yzero)*fcslmn
        ys(ibl,j) = (ym(j)+thick2(j)*fypr)*fcslmn
        xp(ibl,j) = (xm(j)+thick2(j)*yprime*fypr+yzero)*fcslmn
        yp(ibl,j) = (ym(j)-thick2(j)*fypr)*fcslmn
        am(j)  = ATAN(yprime)*c1
        xxm(j) = xm(j)
        IF(j == n) stager = ATAN(ym(n)/xm(n))*c1
        IF(j == n) stag(ibl) = stager
        xm(j) = (xm(j)+yzero)*fcslmn
        ym(j) = ym(j)*fcslmn
        thick2(j) = thick2(j)*fcslmn
        s(j) = s(j)*fcslmn
    END DO
    IF (isplit == 0) GO TO 530
    xspltm(1) = 1. - perspj
    k1 = 25
    xspltm(k1) = 1.
    k11 = k1 - 1
    delxx = perspj/FLOAT(k11)
    DO  j = 2,k11
        xspltm(j) = xspltm(j-1) + delxx
    END DO
    CALL alg15 (xm,ym,n,xspltm,yspltm,k1,1)
    yle = yspltm(1)
    CALL alg15 (xm,s,n,xspltm,sspltm,k1,1)
    CALL alg15 (xm,am,n,xspltm,ss1(1,3),k1,1)
    sspls = sspltm(1)
    DO  j = 1,k1
        sspltm(j) = sspltm(j) - sspls
    END DO
    SELECT CASE ( isplit )
        CASE (    1)
            GO TO 410
        CASE (    2)
            GO TO 420
    END SELECT
410 xnorms = SQRT((xspltm(k1)-xspltm(1))**2+(yspltm(k1)-yspltm(1))**2)
    chords = xnorms/(1.-yzeros+xnorms*(yzeros+ABS(yones*SIN(beta2/ c1))))
    fcslms = (perspj-chords*(yzeros+ABS(yones*SIN(beta2/c1))))/perspj
    yzeros = yzeros*chords/fcslms
    yones  = yones *chords/fcslms
    ts  = ts*chords/fcslms
    at  = (yzeros-ts/2.)/(2.*zspmxt**3)
    ct  = (ts/2.-yzeros)*3./(2.*zspmxt)
    dt  = yzeros
    et  = (yones-ts/2.)/(1.-zspmxt)**3-1.5*(yzeros-ts/2.)/  &
        (zspmxt**2*(1.-zspmxt))
    ft  = 1.5*(yzeros-ts/2.)/zspmxt**2
    ht  = ts/2.
    GO TO 450
420 yzs = yzeros
    ts1 = ts
    beta1 = ss1(1,3)
    y1  =-COS(beta1/c1)/(SIN(beta1/c1)-SIN(beta2/c1))
    x1  = SIN(beta1/c1)/(SIN(beta1/c1)-SIN(beta2/c1))
    rdius = ABS(1./(SIN(beta1/c1)-SIN(beta2/c1)))
    y2  = TAN((beta1+beta2)/(2.*c1))
    xmlcs  = SQRT(1.+y2**2)
    chords = xmlcs/(1.0-2.*yzeros*(1.0-xmlcs))
    fcslms = 1.0-chords*2.*yzeros
    yzeros = yzeros*chords/fcslms
    ts  = ts*chords/fcslms
    ss1(1,1) = 0.
    delx  = 1./(xn-1.)
    t2    = ts/2.
    tprim2= t2-yzeros
    c2    = 2.*c1
    aform = (tprim2+rdius*(1.-COS((beta1-beta2)/c2)))/xmlcs*2.
    phis  = ACOS((1.-aform**2)/(1.+aform**2))
    rs    = yzeros + xmlcs/2./SIN(phis)
    yss   = rdius - rs + t2
    bform = (rdius*(1.-COS((beta1-beta2)/c2))-tprim2)/xmlcs*2.
    phip  = ACOS((1.-bform**2)/(1.+bform**2))
    rp    = xmlcs/2./SIN(phip)-yzeros
    ypp   = rdius - rp - t2
    xx = 0.
    DO  j = 2,n
        xx = xx + delx
        phi1 = ATAN(-1./SQRT(rdius**2-(xx-x1)**2)*(xx-x1))
        IF (beta1 < 0.) phi1 = -phi1
        phi2 = ABS(beta1/c1-phi1)
        ss1(j,1) = rdius*phi2
    END DO
    DO  j = 1,n
        ss1(j,1) = ss1(j,1)/ss1(n,1)
        phix = (ss1(j,1)-.5)*phi2
        ss1(j,2) =(yss*COS(phix)+SQRT(rs**2-yss**2*SIN(phix)**2)-rdius)/t2
    END DO
    CALL alg14 (xspltm,yspltm,k1,xspltm,xdum,ss1(1,3),k1,1)
    xnorms = SQRT(perspj**2+(yspltm(k1)-yspltm(1))**2)
    chords = xnorms/(1.-2.*yzs*(1.-xnorms))
    fcslms = (perspj-chords*2.*yzs)/perspj
    ts     = ts1*chords/fcslms
    yzeros = yzs*chords/fcslms
    450   DO  j = 1,k1
        sn = sspltm(j)/sspltm(k1)
        IF (isplit  > 1) GO TO 480
        IF (sn > zspmxt) GO TO 460
        thick(j) = (at*sn**2+ct)*sn + dt
        GO TO 470
460     sn = sn - zspmxt
        thick(j) = (et*sn+ft)*sn**2 + ht
470     fypr   = 1./SQRT(1.+TAN(ss1(j,3)/c1)**2)
        yprime = TAN(ss1(j,3)/c1)
        GO TO 490
480     CALL alg15 (ss1,ss1(1,2),n,sn,thick(j),1,1)
        thick(j) = thick(j)*ts/2.
        fypr   = 1.0/SQRT(1.0+ss1(j,3)**2)
        yprime = ss1(j,3)
490     xspltp(j) = (xspltm(j)-(1.-perspj)+thick(j)*yprime*fypr+yzeros)*  &
            fcslms+(1.-perspj)
        xsplts(j) = (xspltm(j)-(1.-perspj)-thick(j)*yprime*fypr+yzeros)*  &
            fcslms+(1.-perspj)
        yspltp(j) = (yspltm(j)-yle-thick(j)*fypr)*fcslms+yle
        ysplts(j) = (yspltm(j)-yle+thick(j)*fypr)*fcslms+yle
        xspltm(j) = (xspltm(j)-(1.-perspj)+yzeros)*fcslms+(1.-perspj)
        yspltm(j) = (yspltm(j)-yle)*fcslms+yle
        thick(j)  = thick(j)*fcslms
        sspltm(j) = sspltm(j)*fcslms
    END DO
    IF (isplit > 1) ss1(1,3) = ATAN(ss1(1,3))*c1
    yzeros = yzeros*fcslms
    areas  = pi/2.*yzeros**2
    area2  = areas
    yint   =-4./(3.*pi)*yzeros*areas*SIN(ss1(1,3)/c1)
    xint   = yzeros*(1.-COS(ss1(1,3)/c1)*4./(3.*pi))*areas
    DO  j = 2,k1
        dela  = (thick(j)+thick(j-1))*(sspltm(j)-sspltm(j-1))
        areas = areas + dela
        xint  = xint + dela*(xspltm(j)+xspltm(j-1))/2.
        yint  = yint + dela*(yspltm(j)+yspltm(j-1))/2.
    END DO
    IF (isplit < 2) GO TO 520
    xint  = xint + area2*(xspltm(k1)+4.*yzeros/(3.*pi)*COS(beta2/c1))
    yint  = yint + area2*(yspltm(k1)+4.*yzeros/(3.*pi)*SIN(beta2/c1))
    areas = areas + area2
520 xbars = xint/areas
    ybars = yint/areas
530 CONTINUE
    yzero = yzero*fcslmn
    IF (inast == 0) GO TO 550
    nasnum = irte - irle + 1
    CALL alg15 (x,ss,100,xhere(irle),snadum(irle),nasnum,1)
    sndum1 = snadum(irle)
    sndum2 = snadum(irte)
    DO  j = irle,irte
        snadum(j) = (snadum(j)-sndum1)/(sndum2-sndum1)
        CALL alg15 (xxm,thick2,n,snadum(j),tharr(ibl,j),1,1)
        tharr(ibl,j) = tharr(ibl,j)*2.*axialc
    END DO
550 CONTINUE
    area = pi/2.0*yzero**2
    xint = yzero*(1.0-COS(beta1/c1)*4.0/(3.0*pi))*area
    yint =-4.0/(3.0*pi)*yzero*area*SIN(beta1/c1)
    DO  j = 2,n
        dela = (thick2(j)+thick2(j-1))*(s(j)-s(j-1))
        area = area + dela
        xint = xint + dela*(xm(j)+xm(j-1))/2.0
        yint = yint + dela*(ym(j)+ym(j-1))/2.0
    END DO
    IF (isecn /= 2) GO TO 570
    area2= pi/2.*yzero**2
    xint = xint + area2*(xm(n)+4.*yzero/(3.*pi)*COS(beta2/c1))
    yint = yint + area2*(ym(n)+4.*yzero/(3.*pi)*SIN(beta2/c1))
    area = area + area2
570 xbar = xint/area
    ybar = yint/area
    xbarb= xbar
    ybarb= ybar
    ybar = ybar + ydel/axialc
    xbar = xbar + xdel/axialc
    ax   = 1./99.
    dx(1)= 0.
    DO  ik = 2,100
        dx(ik) = dx(ik-1) + ax
    END DO
    ymm = 0.0
    xmm = 0.0
    DO  ik = 1,100
        xab = dx(ik)
        IF (isecn == 0) GO TO 590
        IF (isecn == 1) GO TO 600
        IF (isecn == 2) GO TO 640
        IF (isecn == 3) GO TO 620
590     y(ik) = (oa48*(xab-h)**4+xab**2*xk2+b*xab+c)*fcslmn
        ss1(ik,1) = oa/12.*(xab-h)**3+xk2*2.*xab+b
        GO TO 660
600     xd = xab - sq
        IF (xd > 0.) GO TO 610
        y(ik) = f3(a1,b1,cc1,e1)*fcslmn
        ss1(ik,1) = f5(a1,b1,cc1)
        GO TO 660
610     y(ik) = f3(a2,b2,cc2,d2)*fcslmn
        ss1(ik,1) = f5(a2,b2,cc2)
        GO TO 660
620     IF (xab-sq > 0.0 .OR. xab == 0.0 .AND. sq == 0.0) GO TO 630
        IF (beta1 == beta3) GO TO 650
        rdius = rdius1
        x1 = x11
        y1 = y11
        bta1 = beta1
        bta2 = beta3
        GO TO 640
630     IF (beta2 == beta3) GO TO 650
        rdius = rdius2
        x1 = x12
        y1 = y12
        bta1 = beta3
        bta2 = beta2
640     y(ik) = f6(xab)*fcslmn
        ss1(ik,1) = f8(xab)
        IF (bta1-bta2 < 0.0) ss1(ik,1) = -ss1(ik,1)
        IF (bta1-bta2 < 0.0) y(ik) = f7(xab)*fcslmn
        GO TO 660
650     ss1(ik,1) = TAN(beta3/c1)
        IF (ik /= 1) ymm = y(ik-1)/fcslmn
        IF (ik /= 1) xmm = dx(ik-1)
        y(ik) = (ss1(ik,1)*(xab-xmm)+ymm)*fcslmn
660     sigma(ik) = dx(ik)*fcslmn + yzero
    END DO
    CALL alg15 (sigma,y,100,dx,dy,100,1)
    CALL alg15 (sigma,ss1(1,1),100,dx,y,100,1)
    CALL alg15 (dx,dy,100,xbar,xab,1,1)
    CALL alg15 (dx,y,100,xbar,xbc,1,1)
    xbar = xbarb
    ybar = ybarb
    ix   = 0.0
    iy   = 0.0
    ixy  = 0.0
    DO  j = 2,n
        dela = (thick2(j)+thick2(j-1))*(s(j)-s(j-1))
        ixd  = (thick2(j)+thick2(j-1))**3*(s(j)-s(j-1))/12.0
        iyd  = (thick2(j)+thick2(j-1))*(s(j)-s(j-1))**3/12.0
        cosang = COS((am(j)+am(j-1))/c1)
        ixn  = (ixd+iyd+(ixd-iyd)*cosang)/2.0
        iyn  = (ixd+iyd-(ixd-iyd)*cosang)/2.0
        ixyn = 0.0
        IF (am(j)+am(j-1) /= 0.0) ixyn = ((ixn-iyn)*cosang-ixd+iyd)/  &
            (2.0*SIN((am(j)+am(j-1))/c1))
        ix  = ix + ixn + dela*((ym(j)+ym(j-1))/2.0-ybar)**2
        iy  = iy + iyn + dela*((xm(j)+xm(j-1))/2.0-xbar)**2
        ixy = ixy+ ixyn+ dela*(ybar-(ym(j)+ym(j-1))/2.0)*(xbar-(xm(j)+  &
            xm(j-1))/2.0)
    END DO
    ang = ATAN(2.0*ixy/(iy-ix))
    ipx = (ix+iy)/2.0+(ix-iy)/2.0*COS(ang)-ixy*SIN(ang)
    ipy = (ix+iy)/2.0-(ix-iy)/2.0*COS(ang)+ixy*SIN(ang)
    ang = ang/2.0*c1
    xml = xm(n)
    yml = ym(n)
    camber = beta1 - beta2
    IF (iprint >= 2) GO TO 790
    lnct = 47
    IF (isecn == 1 .OR. isecn == 3) lnct = 49
    WRITE (log2,680) chord,stager,camber,area,xbar,ybar,ix,iy,ixy,ang,  &
        ipx,ang,ipy,ang
680 FORMAT ( /16X,'NORMALISED results - all the following refer to ablade&
                     &having a meridional chord projection of unity'    ,       &
        /16X,100(1H*),//20X,'BLADE chord',4X,'=',f7.4, //20X,'STAGGER&
                     &angle  ='    ,f7.3 , //20X,'CAMBER angle   =',f7.3, //20X,&
        'SECTION area   =',f7.5,//20X,'LOCATION of centroid&
                     & relative to leading edge'    , //30X,'XBar =',f8.5, /30X,&
        'YBAR =',f8.5, //20X,'SECOND moments of area about &
                     &centroid'    , //30X,'IX   =',f8.5, /30X,'IY   =',f8.5, &
        /30X,'IXY  =',f8.5, //20X,'ANGLE of inclination of (one)&
                     & principal axis to x  axis ='    ,f7.3, //20X,'PRINCIPAL &
                     &second moments of area about centroid'    , //30X,'IPX  =',&
        f7.5,6X,'(at',f7.3,' with  x  axis)', /30X,'IPY  =',f7.5,&
        6X,'(at',f7.3,' with  y  axis)', //)

690 FORMAT (27X,'POINT',8X,'M e a n l i n e  d a t a',13X,'SURFACE coordinate&
               &DATA'    , /27X,'NUMBER',5X,'X',7X,'Y',5X,'ANGLE thickness',&
        9X,'X1',6X,'Y1',6X,'X2',6X,'Y2', //)
    WRITE (log2,690)
    DO  j = 1,n
        IF (lnct /= 60) GO TO 700
        WRITE (log2,10)
        WRITE (log2,690)
        lnct = 4
700     lnct = lnct + 1
        tm   = thick2(j)*2.0
        WRITE (log2,720) j,xm(j),ym(j),am(j),tm,xs(ibl,j),ys(ibl,j),  &
            xp(ibl,j),yp(ibl,j)
    END DO
720 FORMAT (27X,i3,f13.5,f8.5,f7.3,f8.5,f16.5,3F8.5)
    IF (isplit == 0) GO TO 760
    IF (lnct  <= 40) GO TO 730
    WRITE (log2,10)
    lnct = 1
730 WRITE  (log2,740)
740 FORMAT (//10X,'SPLITTER coordinates', /10X,21(1H*), //2X)
    WRITE (log2,690)
    lnct = lnct + 11
    n  = k1
    DO  j = 1,n
        tm = thick(j)*2.
        xs(ibl,j) = xsplts(j)
        xp(ibl,j) = xspltp(j)
        yp(ibl,j) = yspltp(j)
        ys(ibl,j) = ysplts(j)
        WRITE (log2,720) j,xspltm(j),yspltm(j),ss1(j,3),tm,xs(ibl,j),  &
            ys(ibl,j),xp(ibl,j),yp(ibl,j)
        lnct = lnct + 1
        IF (lnct <= 60) CYCLE
        WRITE (log2,10)
        WRITE (log2,690)
        lnct = 4
    END DO
760 CONTINUE
    DO  j = 1,n
        xm(j) = xs(ibl,j)
        ym(j) = ys(ibl,j)
        am(j) = xp(ibl,j)
        thick2(j) = yp(ibl,j)
    END DO
    WRITE  (log2,780) ibl
780 FORMAT (1H1,45X,'NORMALISED plot of section NUMBER',i3, /2X)
    CALL alg16 (n,log2,xm,ym,am,thick2)
790 a2  = axialc**2
    a4  = a2**2
    ix  = ix*a4
    iy  = iy*a4
    ixy = ixy*a4
    ipx = ipx*a4
    ipy = ipy*a4
    IF (istak > 1) GO TO 800
    xbar= istak
    IF (istak == 0) ybar = 0.
    IF (istak == 1) ybar = yml
800 rle = yzero*axialc
    IF (isplit /= 0) GO TO 810
    chord = chord*axialc
    ccord(ibl) = chord
    area = area*a2
    xc  = rle - xbar*axialc - xdel
    yc  =-ybar*axialc - ydel
    xtc = (xml-xbar)*axialc - xdel
    ytc = (yml-ybar)*axialc - ydel
    GO TO 860
810 rle = yzeros*axialc
    chord = chords*axialc
    areas = areas*axialc**2
    xc  = (xspltm(1)-xbar)*axialc - xdel
    yc  = (yspltm(1)-ybar)*axialc - ydel
    xtc = (xspltm(k1)-xbar)*axialc - xdel
    ytc = (yspltm(k1)-ybar)*axialc - ydel
    xbars = (xbars-xbar)*axialc - xdel
    ybars = (ybars-ybar)*axialc - ydel
    IF (iprint >= 2) GO TO 940
    SELECT CASE ( isplit )
        CASE (    1)
            GO TO 820
        CASE (    2)
            GO TO 840
    END SELECT
820 WRITE  (log2,830) chord,rle,xc,yc,xbars,ybars,areas
830 FORMAT (1H1,31X,'DIMENSIONAL results - all results refer TO a blade of &
               &specified chord'    , /32X,69(1H*), //20X,'BLADE chord',4X,'=',  &
        1P,e12.5,//20X,'END radius',5X,'=',1P,e12.5,8X,'CENTERED at x=',  &
        1P,e12.5,3H y=,1P,e13.5, //20X,'LOCATION of centroid at x=',  &
        1P,e12.5,' AND y=',1P,e12.5, //20X,'SECTION area   =',1P,e12.5, //2X)
    GO TO 900
840 WRITE  (log2,850) chord,rle,xc,yc,xtc,ytc,xbars,ybars,areas
850 FORMAT (1H1,31X,'DIMENSIONAL results - all results refer to a blade of &
               &specified chord'    , /32X,69(1H*), //20X,'BLADE chord',4X,1H=,  &
        1P,e12.5,//20X,'END radius',5X,1H=,1P,e12.5,8X,'CENTERED at x=',  &
        1P,e12.5,' y=',1P,e13.5, /64X,'AND x=',1P,e12.5,' y=',1P,e13.5,  &
        /20X,'LOCATION of centroid at x=',1P,e12.5,' AND y=',1P,e12.5,  &
        //20X,'SECTION area   =',1P,e12.5, //2X)
    GO TO 900
860 CONTINUE
    IF (iprint >= 2) GO TO 940
    IF (isecn  == 2) GO TO 880
    WRITE  (log2,870) chord,rle,xc,yc,area,ix,iy,ixy,ipx,ang,ipy,ang
870 FORMAT (1H1,31X,'DIMENSIONAL results - all results refer to a blade of &
             &specified chord'    , /32X,69(1H*),//20X,'BLADE chord',4X,'=',  &
        1P,e12.5,//20X,'L.e.radius',5X,'=',1P,e12.5,8X,'CENTERED at x=',  &
        1P,e13.5,3H y=,1P,e13.5, //20X,'SECTION area   =',1P,e12.5,//20X,  &
        'SECOND moments of area about centroid', //30X,'IX   =',1P,e12.5,  &
        /30X,'IY   =',1P,e12.5, /30X,'IXY  =',1P,e12.5, //20X,'PRINCIPAL &
             &second moments of area about centroid'    , //30X,'IPX  =',1P,e12.5,  &
        '  (at',0P,f7.3,' with  x  axis)', /30X,'IPY  =',1P,e12.5,  &
        '  (at',0P,f7.3,' with  y  axis)', //)
    GO TO 910
880 CONTINUE
    WRITE (log2,890) chord,rle,xc,yc,xtc,ytc,area,ix,iy,ixy,ipx,ang, ipy,ang
890 FORMAT (1H1,31X,'DIMENSIONAL results - all results refer to a blade of &
             &specified chord'    , /32X,69(1H*),//20X,'BLADE chord',4X,'=',  &
        1P,e12.5, //20X,'END radii',6X,'=',1P,e12.5,8X,'CENTERED at x=',  &
        1P,e13.5,' y=',1P,e13.5, /64X,'AND x=',1P,e13.5,' y=',1P,e13.5,  &
        /20X,'SECTION area   =',1P,e12.5, //20X,'SECOND moments of are &
             &a about centroid'    , //30X,'IX   =',1P,e12.5, /30X,'IY   =',1P,e12.5,&
        /30X,'IXY  =',1P,e12.5, //20X,'PRINCIPAL second moments of area &
             &about centroid'    , //30X,'IPX  =',1P,e12.5,'  (at',0P,f7.3,  &
        ' with  x  axis)', /30X,'IPY  =',1P,e12.5,'  (at',0P,f7.3,  &
        ' with  y  axis)', //)
900 CONTINUE
910 WRITE (log2,920)
    WRITE (log2,930)
920 FORMAT(4X,'PT',5X,'SURFACE',10(1H-),'ONE',8X,'SURFACE',10(1H-),'Two',&
        10X,'PT',5X,'SURFACE',10(1H-),'ONE',8X,'SURFACE',10(1H-),'TWO')
930 FORMAT (4X,'NO',8X,'X',13X,'Y',13X,'X',13X,'Y',12X,'NO',8X,'X',13X,  &
        'Y',13X,'X',13X,'Y', //)
    lnct = 24
    940   DO  j = 1,n
        xs(ibl,j) = (xs(ibl,j) - xbar)*axialc - xdel
        ys(ibl,j) = (ys(ibl,j) - ybar)*axialc - ydel
        xp(ibl,j) = (xp(ibl,j) - xbar)*axialc - xdel
        yp(ibl,j) = (yp(ibl,j) - ybar)*axialc - ydel
        IF (iprint  >= 2) CYCLE
        IF ((j/2)*2 /= j) CYCLE
        IF (lnct   /= 60) GO TO 950
        lnct = 4
        WRITE (log2,10)
        WRITE (log2,920)
        WRITE (log2,930)
950     lnct = lnct + 1
        jm1  = j - 1
        WRITE (log2,960) jm1,xs(ibl,jm1),ys(ibl,jm1),xp(ibl,jm1),  &
            yp(ibl,jm1),j,xs(ibl,j),ys(ibl,j),xp(ibl,j), yp(ibl,j)
960     FORMAT (3X,i3,4(2X,1P,e12.5),6X,i3,4(2X,1P,e12.5))
    END DO
    chordd(ibl) = chord
    IF (isplit > 1) isecn = isplit
    IF (iprint >= 2) GO TO 1000
    IF (lnct  > 24) WRITE (log2,980)
980 FORMAT (1H1)
    IF (lnct > 24) lnct = 2
    lnct = lnct + 5
    IF (isecn == 2) GO TO 1030
    WRITE  (log2,990)
990 FORMAT (//48X,'POINTS describing leading edge radius', //48X,  &
        'POINT no.',6X,'X',13X,'Y', /2X)
1000 eps = beta1 + 180.0
    IF (isecn == 2) GO TO 1030
    DO  j = 1,31
        xsemi(ibl,j) = xc - rle*SIN(eps/c1)
        ysemi(ibl,j) = yc + rle*COS(eps/c1)
        eps = eps - 6.0
        IF (iprint >= 2) CYCLE
        WRITE (log2,1010) j,xsemi(ibl,j),ysemi(ibl,j)
        lnct = lnct + 1
1010    FORMAT (48X,i5,1P,e17.5,1P,e14.5)
    END DO
    GO TO 1090
1030 phiss = phis - ABS((beta1-beta2)/c2)
    phipp = ABS((beta1-beta2))/c2 - phip
    eps   = beta1 + 180.0
    eps2  = beta2 + 90.
    delep = (180.-(phiss+phipp)*c1)/28.
    DO  j = 1,31
        IF (j /= 1) GO TO 1040
        xsemi(ibl,j) = xp(ibl,1)
        ysemi(ibl,j) = yp(ibl,1)
        xsemj(ibl,j) = xs(ibl,n)
        ysemj(ibl,j) = ys(ibl,n)
        eps  = eps  - phipp*c1
        eps2 = eps2 - phiss*c1
        CYCLE
1040    IF (j /= 31) GO TO 1050
        xsemi(ibl,j) = xs(ibl,1)
        ysemi(ibl,j) = ys(ibl,1)
        ysemj(ibl,j) = yp(ibl,n)
        xsemj(ibl,j) = xp(ibl,n)
        CYCLE
1050    xsemi(ibl,j) = xc  - rle*SIN(eps/c1)
        ysemi(ibl,j) = yc  + rle*COS(eps/c1)
        xsemj(ibl,j) = xtc + rle*COS(eps2/c1)
        ysemj(ibl,j) = ytc + rle*SIN(eps2/c1)
        eps  = eps  - delep
        eps2 = eps2 - delep
    END DO
    IF (iprint >= 2) GO TO 1090
    WRITE  (log2,1070)
1070 FORMAT (//39X,'POINTS describing leading AND trailing edges',  &
        /25X,'LEADING edge',22X,'TRAILING edge', /2X,'POINT no.',4X,8X,  &
        'X',14X,'Y',12X,8X,'X',14X,'Y', /2X)
    WRITE (log2,1080) (j,xsemi(ibl,j),ysemi(ibl,j),xsemj(ibl,j),  &
        ysemj(ibl,j),j=1,31)
    lnct = lnct + 31
1080 FORMAT (6X,i2,7X,1P,e17.5,1P,e14.5,2X,1P,e17.5,1P,e14.5)
1090 ssurf = axialc
    ss2   =  bx - axialc*xbar  - xdel
    sbar  = ss2 + axialc*xbarb + xdel
    DO  ik = 1,100
        ss(ik) = ss(ik) - sbar
    END DO
    CALL alg15 (ss,x,100,0.0,sbar,1,1)
    CALL alg15 (xhere,r(1,ibl),nstns,sbar,rxbar,1,0)
    xbarc = xbar
    ybarc = ybar
    xbar  = xbarb + xdel/axialc
    ybar  = ybarb + ydel/axialc
    ss1(1,1) = ss(1)
    s23 = axialc/99.
    ss(1) = ss(1) + ss2
    DO  ik = 2,100
        ss1(ik,1) = ss(ik)
        ss(ik) = ss(ik-1) + s23
    END DO
    sigmao = (xab-ybar)/rxbar*axialc
    DO  ik = 2,100
        IF (xbar == dx(ik)) GO TO 1140
        IF (xbar > dx(ik-1) .AND. xbar < dx(ik)) GO TO 1150
    END DO
    WRITE  (log2,1130)
1130 FORMAT (1H1,23H xbar cannot be located)
1140 sigma(ik) = sigmao
    kl = ik + 1
    GO TO 1160
1150 kl = ik
    sigma(ik-1) = sigmao
1160 ssdum = ss(kl-1)
    ss(kl-1) = 0.
    yp1 = xbc
    rx1 = rxbar
    DO  ik = kl,100
        xsurf = ss2 + dx(ik)*ssurf + ss1(1,1)
        CALL alg15 (ss1(1,1),x,100,xsurf,xdum,1,1)
        CALL alg15 (xhere,r(1,ibl),nstns,xdum,rx2,1,0)
        sigma(ik) = sigma(ik-1) + (y(ik)/rx2+yp1/rx1)/2.*(ss(ik)-ss(ik-1))
        yp1 = y(ik)
        rx1 = rx2
    END DO
    ss(kl-1) = ssdum
    ssdum  = ss(kl)
    sigdum = sigma(kl)
    sigma(kl) = sigmao
    ss(kl) = 0.
    rx1 = rxbar
    yp1 = xbc
    km  = kl - 1
    DO  ik = 1,km
        kj  = kl - ik
        xsurf = ss2 + dx(kj)*ssurf + ss1(1,1)
        CALL alg15 (ss1(1,1),x,100,xsurf,xdum,1,1)
        CALL alg15 (xhere,r(1,ibl),nstns,xdum,rx2,1,0)
        sigma(kj) = sigma(kj+1)-(y(kj)/rx2+yp1/rx1)/2.*(ss(kj+1)-ss(kj))
        yp1 = y(kj)
        rx1 = rx2
    END DO
    sigma(kl) = sigdum
    ss(kl) = ssdum
    DO  ik = 1,100
        ss(ik) = ss1(ik,1)
    END DO
    xbar = xbarc
    ybar = ybarc
    DO  ik = 1,n
        ss1(ik,1) = ss2 + ((xs(ibl,ik)+xdel)/axialc+xbar)*ssurf+ss(1)
        ss1(ik,2) = ss2 + ((xp(ibl,ik)+xdel)/axialc+xbar)*ssurf+ss(1)
    END DO
    DO  ik =1,31
        ss1(ik,3) = ss2 + ((xsemi(ibl,ik)+xdel)/axialc+xbar)*ssurf+ss(1)
    END DO
    IF (isecn /= 2) GO TO 1230
    DO  ik = 1,31
        ss1(ik,4) = ss2 + ((xsemj(ibl,ik)+xdel)/axialc+xbar)*ssurf+ss(1)
    END DO
    CALL alg15 (ss,x,100,ss1(1,4),ss1(1,4),31,1)
1230 CALL alg15 (ss,x,100,ss1(1,1),ss1(1,1),n,1)
    CALL alg15 (ss,x,100,ss1(1,2),ss1(1,2),n,1)
    CALL alg15 (ss,x,100,ss1(1,3),ss1(1,3),31,1)
    IF (istak > 1) GO TO 1250
    IF (istak == 1) sigmao = sigma(100)
    IF (istak == 0) sigmao = sigma(1)
    DO  ik = 1,100
        sigma(ik) = sigma(ik) - sigmao
    END DO
    1250  DO  ik = 1,100
        dx(ik) = (dx(ik)-xbar)*axialc - xdel
        dy(ik) = (dy(ik)-ybar)*axialc - ydel
    END DO
    DO  mk = 1,4
        IF (isecn /= 2 .AND. mk == 4) CYCLE
        IF (mk == 4 .OR. mk == 3) nnn = 31
        IF (mk == 1 .OR. mk == 2) nnn = n
        DO  ik = 1,nnn
            IF (mk == 1) yp1 = ys(ibl,ik)
            IF (mk == 2) yp1 = yp(ibl,ik)
            IF (mk == 3) yp1 = ysemi(ibl,ik)
            IF (mk == 4) yp1 = ysemj(ibl,ik)
            IF (mk == 1) rx1 = xs(ibl,ik)
            IF (mk == 2) rx1 = xp(ibl,ik)
            IF (mk == 3) rx1 = xsemi(ibl,ik)
            IF (mk == 4) rx1 = xsemj(ibl,ik)
            CALL alg15 (dx,dy,100,rx1,rxbar,1,1)
            delly = yp1 - rxbar
            CALL alg15 (xhere,r(1,ibl),nstns,ss1(ik,mk),rab,1,0)
            delsig = delly/rab
            CALL alg15 (dx,sigma,100,rx1,xab,1,1)
            ss1(ik,mk) = xab + delsig
        END DO
    END DO
    RETURN

1290 WRITE  (log2,1300)
1300 FORMAT (1H1,10X,'ITERATIVE solution for constant fails - case abandoned')
    CALL mesage (-37,0,NAME)

END SUBROUTINE alg13
