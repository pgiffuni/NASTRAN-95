SUBROUTINE alg11
     
 REAL :: loss,lami,lamip1,lamim1
 DIMENSION       xm(21),ppg(21),v(21),pt(21),ps(21),wt(21),pn(21),  &
     p1(21),deltp(21,30),ts(21),solid(21),deltb(21),  &
     tr(21,30),rmdv(21,6),idata(6),rdata(6),name1(2), name2(2)
 COMMON /system/ ksystm(90),lpunch
 COMMON /ud3prt/ iprtc,istrml,ipgeom
 COMMON /algino/ iscr
 COMMON /udstr2/ nbldes,stag(21),chordd(21)
 COMMON /udsign/ nsign
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
 EQUIVALENCE    (idata(1),rdata(1))
 DATA    name1, name2 /4HPLOA,4HD2  ,4HTEMP,4H    /
 
 opr   = 0.0
 oeff  = 1.0
 pfac  = 550.0
 ilast = nstns
 
!     LOCATE COMPUTING STATION NUMBER AT THE BLADE LEADING EDGE AND
!     AT THE BLADE TRAILING EDGE.
 
 ledgeb = 0
 itrleb = 0
 DO  ible = 1,nstns
   nout3s = nout3(ible)/10
   IF (nout3(ible) == 1 .OR. nout3s == 1) ledgeb = ible
   IF (nout3(ible) == 2 .OR. nout3s == 2) itrleb = ible
 END DO
 IF (ifailo /= 0) ilast = ifailo
 DO  i = 1,ilast
   CALL alg03 (lnct,7+nstrms)
   IF (iprtc == 1) WRITE(log2,100) i
   100 FORMAT (//10X,'STATION',i3,'  FLOW-FIELD DESCRIPTION', /10X,  &
       34(1H*), //,'  STREAM      -----MESH-POINT COORDS------',  &
       3X,16(1H-),'V E L O C I T I E S,16(1H-)    RADIUS OF  ',  &
       'STREAMLINE   STATION',/,'  -LINE       RADIUS    X-COORD'  &
       ,       '    L-COORD   MERIDIONAL TANGENTIAL   AXIAL',6X,'RADIAL',  &
       4X,'TOTAL    CURVATURE SLOPE ANGLE LEAN ANGLE',/)
   CALL alg01 (r(1,i),x(1,i),nstrms,r(1,i),x1,gama,nstrms,0,1)
   IF (i /= 1 .AND. i /= nstns) GO TO 130
   l1 = 1
   l2 = 2
   IF (i == 1) GO TO 110
   l2 = nstns
   l1 = l2 - 1
   110 DO  j = 1,nstrms
     cr(j)  = 0.0
     phi(j) = ATAN2(r(j,l2)-r(j,l1),x(j,l2)-x(j,l1))
   END DO
   GO TO 150
   130 DO  j = 1,nstrms
     x1 = SQRT((r(j,i+1)-r(j,i))**2+(x(j,i+1)-x(j,i))**2)
     x2 = SQRT((r(j,i)-r(j,i-1))**2+(x(j,i)-x(j,i-1))**2)
     x3 = ATAN2(r(j,i+1)-r(j,i),x(j,i+1)-x(j,i))
     x4 = ATAN2(r(j,i)-r(j,i-1),x(j,i)-x(j,i-1))
     cr(j) = (x3-x4)/(x1+x2)*2.0
     IF (cr(j) /= 0.0) cr(j) = 1.0/cr(j)
     phi(j) = (x3+x4)/2.0
   END DO
   150 DO  j = 1,nstrms
     va = vm(j,i)*COS(phi(j))
     vr = vm(j,i)*SIN(phi(j))
     fi = phi(j)*c1
     ga = ATAN(gama(j))*c1
     ppg(j) = fi + ga
     v(j) = SQRT(vm(j,i)**2+vw(j,i)**2)
     
!     STORE RADIUS AT BLADE LEADING AND TRAILING EDGES, ALL STREAMLINES
     
     IF (icase == 1 .AND. i == ledgeb) rmdv(j,5) = r(j,i)
     IF (icase == 1 .AND. i == itrleb) rmdv(j,6) = r(j,i)
     IF (iprtc == 1) WRITE (log2,170) j,r(j,i),x(j,i),xl(j,i),  &
         vm(j,i),vw(j,i),va,vr,v(j),cr(j),fi,ga
   END DO
   170 FORMAT (i6,f14.4,2F11.4,5F11.2,1X,f10.2,2F11.3)
   CALL alg03 (lnct,nstrms+4)
   IF (iprtc == 1) WRITE (log2,180)
   180 FORMAT (/8H  stream,7X,4HMACH,6X,4(1H-),9HPRESSURES,4(1H-),5X,  &
       17H---temperatures--,4X,8HSPECIFIC,4X,17H---enthalpies----,  &
       4X,7HENTROPY,6X,4HFLOW,3X,11H(phi+gamma), /,7H  -line,7X,  &
       6HNUMBER,5X,5HTOTAL,6X,6HSTATIC,5X,5HTOTAL,6X,6HSTATIC,5X,  &
       6HWEIGHT,5X,5HTOTAL,6X,6HSTATIC,16X,5HANGLE,/)
   DO  j = 1,nstrms
     deltb(j) = 0.0
     hs = h(j,i) - v(j)**2/(2.0*g*ej)
     IF (hs < hmin) hs = hmin
     xm(j) = SQRT(alg9(hs,s(j,i),v(j)**2))
     pt(j) = alg4(h(j,i),s(j,i))
     ptins = pt(j)/sclfac**2
     ps(j) = alg4(hs,s(j,i))
     psins = ps(j)/sclfac**2
     tt    = alg7(h(j,i),s(j,i))
     ts(j) = alg7(hs,s(j,i))
     wt(j) = alg5(hs,s(j,i))
     alpha = 0.0
     IF (i /= istag .OR. j /= 1) alpha = c1*ATAN(vw(j,i)/vm(j,i))
     
!     STORE DENSITY AT BLADE LEADING EDGE FOR ALL STREAMLINES
     
     IF (icase == 1 .AND. i == ledgeb) rmdv(j,2) = wt(j)
     IF (iprtc == 1) WRITE (log2,200) j,xm(j),ptins,psins,tt,ts(j),  &
         wt(j),h(j,i),hs,s(j,i),alpha,ppg(j)
   END DO
   200 FORMAT (i6,f14.4,2F11.4,2F11.3,f12.6,f10.3,f11.3,f12.6,f10.3, f11.3)
   IF (i /= 1) GO TO 220
   p1bar = 0.0
   h1bar = 0.0
   p1(1) = pt(1)
   pn(1) = pt(1)
   DO  j = 1,itub
     p1(j+1) = pt(j+1)
     pn(j+1) = pt(j+1)
     x1    = (delf(j+1)-delf(j))/2.0
     p1bar = p1bar + x1*(pt(j)+pt(j+1))
     h1bar = h1bar + x1*(h(j,1)+h(j+1,1))
   END DO
   hbar  = h1bar
   s1bar = alg3(p1bar,h1bar)
   pnbar = p1bar
   hnbar = h1bar
   snbar = s1bar
   l1keep= 1
   CYCLE
   220 ifle  = 0
   ifte  = 0
   IF (nwork(i) == 0) GO TO 230
   ifte  = 1
   IF (i == nstns .OR. nwork(i+1) == 0 .OR. speed(i) == speed(i+1)) GO TO 240
   ifle  = 1
   GO TO 240
   230 IF (i == nstns .OR. nwork(i+1) == 0) GO TO 240
   ifle  = 1
   240 IF (ifte == 0) GO TO 500
   CALL alg03 (lnct,nstrms+8)
   xn    = speed(i)*spdfac(icase)
   xblade= 10.0
   IF (nblade(i) /= 0) xblade = ABS(FLOAT(nblade(i)))
   l1    = xblade
   IF (iprtc == 1) WRITE (log2,250) i,xn,l1
   250 FORMAT (/10X,'STATION',i3,' IS WITHIN OR AT THE TRAILING EDGE OF',  &
       ' A BLADE ROTATING AT',f8.1,' RPM  NUMBER OF BLADES IN ',  &
       'ROW =',i3, /10X,109(1H*), //,'  STREAM      BLADE     ',  &
       'RELATIVE    RELATIVE   RELATIVE  DEVIATION    BLADE    ',  &
       '  LEAN    PRESS DIFF    LOSS     DIFFUSION   DELTA P',  &
       /,'  -LINE',7X,'SPEED     VELOCITY    MACH NO.  FLOW',  &
       ' ANGLE   ANGLE      ANGLE      ANGLE   ACROSS BLADE  ',  &
       'COEFF      FACTOR     ON Q',/)
   q = 1.0
   IF (speed(i) < 0.0) GO TO 290
   IF (speed(i) > 0.0) GO TO 280
   IF (i < 3) GO TO 290
   ii = i - 1
   260 IF (speed(ii) /= 0.0) GO TO 270
   IF (ii == 2) GO TO 290
   ii = ii - 1
   GO TO 260
   270 IF (speed(ii) < 0.0) q = -1.0
   GO TO 290
   280 q  = -1.0
   290 l1 = ndimen(i) + 1
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
   300 DO  j = 1,nstrms
     taneps(j) = r(j,i)
   END DO
   GO TO 380
   320 DO  j = 1,nstrms
     taneps(j) = r(j,i)/r(nstrms,i)
   END DO
   GO TO 380
   340 DO  j = 1,nstrms
     taneps(j) = xl(j,i)
   END DO
   GO TO 380
   360 DO  j = 1,nstrms
     taneps(j) = xl(j,i)/xl(nstrms,i)
   END DO
   380 l1 = is2(i)
   IF (nwork(i) == 5 .OR. nwork(i) == 6) CALL alg01 (datac(l1),  &
       data6(l1),ndata(i),taneps,deltb,x1,nstrms,nterp(i),0)
   CALL alg01 (datac(l1),data5(l1),ndata(i),taneps,solid,x1,nstrms,  &
       nterp(i),0)
   CALL alg01 (datac(l1),data3(l1),ndata(i),taneps,taneps,x1,nstrms,  &
       nterp(i),0)
   l1 = i + nl1(i)
   l2 = l1
   IF (nloss(i) == 1 .OR. nloss(i) == 4 .OR. nwork(i) == 7) l2 = i + nl2(i)
   xn = xn*pi/(30.0*sclfac)
   DO  j = 1,nstrms
     u  = xn*r(j,i)
     vr = SQRT(vm(j,i)**2+(vw(j,i)-u)**2)
     xmr  = xm(j)*vr/v(j)
     beta = ATAN(tbeta(j,i))*c1
     bbeta= 0.0
     IF (nwork(i) == 5 .OR. nwork(i) == 6) bbeta = beta - deltb(j)
     deltb(j) = deltb(j)*q
     delp = 0.0
     IF (i == nstns .OR. nwork(i+1) == 0 .OR. speed(i) /= speed(i+1))  &
         GO TO 390
     x1 = SQRT((r(j,i+1)-r(j,i))**2+(x(j,i+1)-x(j,i))**2)
     x2 = SQRT((r(j,i)-r(j,i-1))**2+(x(j,i)-x(j,i-1))**2)
     x3 = xblade
     delp = pi*r(j,i)*wt(j)/(sclfac**2*x3*g)*(tbeta(j,i)/  &
         (1.0+tbeta(j,i)**2)*ts(j)*g*ej*((s(j,i+1)-s(j,i))/x1 +  &
         (s(j,i)-s(j,i-1))/x2)+vm(j,i)/r(j,i)*((r(j,i+1)*vw(j,i+1) -  &
         r(j,i)*vw(j,i))/x1+(r(j,i)*vw(j,i)-r(j,i-1)*vw(j,i-1))/x2))
     deltp(j,i) = delp
     390 hri  = h(j,i) - (v(j)**2-vr**2)/(2.0*g*ej)
     prd  = alg4(hri,s(j,l1))
     pr   = alg4(hri,s(j,i))
     tr(j,i) = alg7(hri,s(j,i))
     prl2 = pr
     psl2 = ps(j)*sclfac**2
     IF (l2 == i) GO TO 400
     prl2 = h(j,l2) - (vw(j,l2)**2 - (vw(j,l2) - xn*r(j,l2))**2)/ (2.0*g*ej)
     prl2 = alg4(prl2,s(j,l2))
     psl2 = h(j,l2) - (vw(j,l2)**2+vm(j,l2)**2)/(2.0*g*ej)
     psl2 = alg4(psl2,s(j,l2))
     400 coef = (prd-pr)/(prl2-psl2)
     dif  = 0.0
     IF (solid(j) == 0.0) GO TO 410
     x2   = vw(j,l1) - xn*r(j,l1)
     x1   = SQRT(vm(j,l1)**2+x2**2)
     x3   = vw(j,i) - u
     dif  = 1.0 - vr/x1 + (x2-x3)/(2.0*x1*solid(j))*q
     410 prl1 = prl2
     psl1 = psl2
     IF (l2 == l1) GO TO 420
     psl1 = h(j,l1) - (vw(j,l1)**2 + vm(j,l1)**2)/(2.0*g*ej)
     prl1 = psl1 + (vm(j,l1)**2 + (vw(j,l1)-xn*r(j,l1))**2)/(2.0*g*ej)
     psl1 = alg4(psl1,s(j,l1))
     prl1 = alg4(prl1,s(j,l1))
     420 dpq  = (ps(j)-psl1)/(prl1-psl1)
     IF (iprtc == 1) WRITE (log2,434) j,u,vr,xmr,beta,deltb(j),bbeta,  &
         taneps(j),delp,coef,dif,dpq
   END DO
   434 FORMAT (i6,f14.2,f11.2,f11.4,4F11.3,f11.4,f11.5,f10.4,f11.4)
   CALL alg03 (lnct,nstrms+5)
   pbar = 0.0
   hbar = 0.0
   DO  j = 1,itub
     x1   = (delf(j+1)-delf(j))/2.0
     pbar = pbar + x1*(pt(j)+pt(j+1))
     hbar = hbar + x1*(h(j,i)+h(j+1,i))
   END DO
   rbar1= pbar/p1bar
   dh1  = (hbar-h1bar)/h1bar
   eff1 = 0.0
   IF (hbar /= h1bar) eff1 = (alg2(s1bar,pbar)-h1bar)/(hbar-h1bar)
   opr  = rbar1
   IF (eff1 /= 0.0) oeff = eff1
   IF (l1 == l1keep) GO TO 460
   l1keep= l1
   pnbar = 0.0
   hnbar = 0.0
   DO  j = 1,nstrms
     pn(j) = alg4(h(j,l1),s(j,l1))
   END DO
   DO  j = 1,itub
     x1    = (delf(j+1)-delf(j))/2.0
     pnbar = pnbar + x1*(pn(j)+pn(j+1))
     hnbar = hnbar + x1*(h(j,l1)+h(j+1,l1))
   END DO
   snbar = alg3(pnbar,hnbar)
   460 effn  = 0.0
   IF (hnbar /= hbar) effn = (alg2(snbar,pbar)-hnbar)/(hbar-hnbar)
   rbarn = pbar/pnbar
   dhn   = (hbar-hnbar)/hnbar
   IF (iprtc == 1) WRITE (log2,470) i,l1,i,i,l1,i,rbar1,rbarn,eff1,  &
       effn,dh1,dhn
   470 FORMAT (/,'  STREAM',7X,'INLET THROUGH STATION',i3,7X,'STATION',  &
       i3,' THROUGH STATION',i3,5X,'MEAN VALUES',6X,  &
       'INLET TO STA.',i2,'   STA.',i2,' TO STA.',i2, /,  &
       '  -LINE',6X,'PRESSURE  ISENTROPIC  DELTA H    PRESSURE  ',  &
       'ISENTROPIC  DELTA H     PRESSURE RATIO',f14.4,f19.4, /15X,  &
       'RATIO   EFFICIENCY  ON H1        RATIO   EFFICIENCY  ON ',  &
       'H1       ISEN EFFY',2F19.4, /80X,'DELTA H ON H1',f15.4, f19.4)
   DO  j = 1,nstrms
     rbar1 = pt(j)/p1(j)
     eff1  = 0.0
     IF (h(j,i) /= h(j,1)) eff1 = (alg2(s(j,1),pt(j))-h(j,1))/  &
         (h(j,i)-h(j,1))
     dh1   = (h(j,i)-h(j,1))/h(j,1)
     rbarn = pt(j)/pn(j)
     effn  = 0.0
     IF (h(j,i) /= h(j,l1)) effn = (alg2(s(j,l1),pt(j))-h(j,l1))/  &
         (h(j,i)-h(j,l1))
     dhn = (h(j,i)-h(j,l1))/h(j,l1)
     IF (iprtc == 1) WRITE (log2,490) j,rbar1,eff1,dh1,rbarn,effn,dhn
   END DO
   490 FORMAT (i6,f14.4,f10.4,f11.4,f12.4,f10.4,f11.4)
   500 IF (ifle == 0) CYCLE
   CALL alg03 (lnct,nstrms+8)
   xn = speed(i+1)*spdfac(icase)
   ip = i + 1
   xblade = 10.0
   IF (nblade(ip) /= 0) xblade = ABS(FLOAT(nblade(ip)))
   l1 = xblade
   IF (iprtc == 1) WRITE (log2,510) i,xn,l1
   510 FORMAT (/10X,'STATION',i3,' IS AT THE LEADING EDGE OF A BLADE ',  &
       'ROATING AT',f9.1,' RPM  NUMBER OF BLADES IN ROW =',i3,  &
       /10X,99(1H*), //,'  STREAM      BLADE     RELATIVE   ',  &
       'RELATIVE   RELATIVE  INCIDENCE    BLADE      LEAN    ',  &
       'PRESS DIFF', /,'  -LINE       SPEED     VELOCITY   MACH',  &
       ' NO.  FLOW ANGLE   ANGLE      ANGLE      ANGLE   ACROSS', ' BLADE',/)
   xn = xn*pi/(30.0*sclfac)
   q  = 1.0
   IF (speed(ip) < 0.0) GO TO 550
   IF (speed(ip) > 0.0) GO TO 540
   IF (ip < 3) GO TO 550
   ii = ip - 1
   520 IF (speed(ii) /= 0.0) GO TO 530
   IF (ii == 2) GO TO 550
   ii = ii - 1
   GO TO 520
   530 IF (speed(ii) < 0.0) q = -1.0
   GO TO 550
   540 q = -1.0
   550 DO  j = 1,nstrms
     cr(j) = 0.0
     taneps(j) = 0.0
   END DO
   IF (nwork(i) /= 0 .OR. ndata(i) == 0) GO TO 660
   l1 = ndimen(i) + 1
   SELECT CASE ( l1 )
     CASE (    1)
       GO TO 570
     CASE (    2)
       GO TO 590
     CASE (    3)
       GO TO 610
     CASE (    4)
       GO TO 630
   END SELECT
   570 DO  j = 1,nstrms
     taneps(j) = r(j,i)
   END DO
   GO TO 650
   590 DO  j = 1,nstrms
     taneps(j) = r(j,i)/r(nstrms,i)
   END DO
   GO TO 650
   610 DO  j = 1,nstrms
     taneps(j) = xl(j,i)
   END DO
   GO TO 650
   630 DO  j = 1,nstrms
     taneps(j) = xl(j,i)/xl(nstrms,i)
   END DO
   650 l1 = is2(i)
   CALL alg01 (datac(l1),data1(l1),ndata(i),taneps,cr,x1,nstrms, nterp(i),0)
   CALL alg01 (datac(l1),data3(l1),ndata(i),taneps,taneps,x1,nstrms,  &
       nterp(i),0)
   660 bbeta = 0.0
   DO  j = 1,nstrms
     u   = xn*r(j,i)
     vr  = SQRT(vm(j,i)**2 + (vw(j,i)-u)**2)
     xmr = xm(j)*vr/v(j)
     tr(j,i) = alg7(h(j,i)-(v(j)**2-vr**2)/(2.0*g*ej),s(j,i))
     beta = ATAN((vw(j,i)-u)/vm(j,i))*c1
     
!     STORE REL. MACH, REL. VEL AND REL. FLOW ANGLE FOR ALL STREAMLINES
!     AT THE BLADE LEADING EDGE
     
     IF (icase /= 1 .OR. i /= ledgeb) GO TO 675
     rmdv(j,1) = xmr
     rmdv(j,3) = vr
     rmdv(j,4) = beta
     675 CONTINUE
     deltb(j) = 0.0
     IF (nwork(i) /= 0 .OR. ndata(i) == 0) GO TO 670
     bbeta = ATAN((TAN(cr(j)/c1)*(1.0-gama(j)*TAN(phi(j))) -  &
         TAN(phi(j))*TAN(taneps(j)/c1)*SQRT(1.0+gama(j)**2))* COS(phi(j)))*c1
     deltb(j) = (beta-bbeta)*q
     670 x1   = SQRT((r(j,i+1)-r(j,i))**2+(x(j,i+1)-x(j,i))**2)
     delp = pi*r(j,i)*2.0*wt(j)/(sclfac**2*xblade*g)*(SIN(beta/c1)*  &
         COS(beta/c1)*g*ej*ts(j)*(s(j,i+1)-s(j,i))/x1+vm(j,i)/  &
         (r(j,i)*x1)*(r(j,i+1)*vw(j,i+1)-r(j,i)*vw(j,i)))
     deltp(j,i) = delp
     IF (iprtc == 1) WRITE (log2,690) j,u,vr,xmr,beta,deltb(j),bbeta,  &
         taneps(j),delp
   END DO
   690 FORMAT (i6,f14.2,f11.2,f11.4,4F11.3,f11.4)
 END DO
 IF (nbl == 0) GO TO 770
 l1 = (ilast-1)/10 + 1
 CALL alg03 (lnct,3+5*l1)
 IF (iprtc /= 1) GO TO 770
 WRITE  (log2,710)
 710 FORMAT (/10X,'ANNULUS WALL BOUNDARY LAYER CALCULATION RESULTS',  &
     /10X,47(1H*))
 DO  k = 1,l1
   l2 = 10*(k-1) + 1
   l3 = l2 + 9
   IF (l3 > ilast) l3 = ilast
   WRITE (log2,730) (i,i=l2,l3)
   WRITE (log2,740) (delh(i),i=l2,l3)
   WRITE (log2,750) (delt(i),i=l2,l3)
   WRITE (log2,760) (wwbl(i),i=l2,l3)
 END DO
 730 FORMAT (/,' STATION NUMBER',14X,10I10)
 740 FORMAT (' HUB DISPLACEMENT THICKNESS',4X,10F10.5)
 750 FORMAT (' CASE DISPLACEMENT THICKNESS',3X,10F10.5)
 760 FORMAT (' BLOCKAGE AREA FRACTION',8X,10F10.5)
 770 CALL alg03 (lnct,4)
 IF (iprtc == 1 .AND. ivfail == 0.AND.iffail == 0) WRITE (log2,780)  &
     icase,ipass
 IF (ifailo /= 0) WRITE (log2,790) icase,ipass,ifailo
 IF (ifailo == 0 .AND. (ivfail /= 0.OR.iffail /= 0))  &
     WRITE (log2,800) icase,ipass,ivfail,iffail
 780 FORMAT (/10X,'POINT NO',i3,'   PASS',i3,'   THE CALCULATION IS ',  &
     'CONVERGED', /10X,52(1H*))
 790 FORMAT (/10X,'POINT NO',i3,'   PASS',i3,'   THE CALCULATION FAIL',  &
     'ED AT STATION',i3, /10X,60(1H*))
 800 FORMAT (/10X,'POINT NO',i3,'   PASS',i3,'   THE CALCULATION IS ',  &
     'NOT FULLY CONVERGED  IVFAIL =',i3,'  IFFAIL =',i3, /10X, 88(1H*))
 power = flow(icase)*(hbar-h1bar)*ej/pfac
 IF (iprtc == 1) WRITE (log2,810) spdfac(icase),flow(icase),opr, oeff,power
 810 FORMAT (10X,'SPEED FACTOR =',f10.3,'  FLOW =',f8.3,'  TOTAL PRES',  &
     'SURE RATIO =',f7.3,'  ISENTROPIC EFFICIENCY =',f6.4, '  POWER =',e11.4)
 IF (iprtc == 0) WRITE (log2,815) icase,ipass,spdfac(icase),  &
     flow(icase),opr,oeff,power
 815 FORMAT  (18H     for point no.,i3,5H pass,i3,15H - speed factor,  &
     10X,1H=,f10.4 / 32X,4HFLOW,18X,1H=,f10.4, /  &
     32X,23HTOTAL pressure ratio  =,f10.4, /32X,'ISENTROPIC ',  &
     'EFFICIENCY =',f10.4, /32X,'POWER',17X,1H=,e10.4)
 IF (ifailo /= 0) GO TO 920
 l1 = 2
 820 DO  i = l1,nstns
   nout3s = nout3(i)/10
   IF (nout3s == 0) nout3s = nout3(i)
   IF (nout3s == 1 .OR. nout3s == 3) GO TO 840
 END DO
 GO TO 920
 840 l2 = i
 l3 = i + 1
 DO  i = l3,nstns
   nout3s = nout3(i)/10
   nout3t = nout3(i) - nout3s*10
   IF (nout3s == 0) nout3t = 1
   IF (nout3s == 0) nout3s = nout3(i)
   IF (nout3s == 2 .OR. nout3s == 3) EXIT
 END DO
 860 l3 = i
 CALL alg03 (lnct,10)
 IF (iprtc == 1) WRITE (log2,870) l2,l3
 870 FORMAT (/10X,'DATA FOR NASTRAN PROGRAM FOR BLADE BETWEEN STATIONS'  &
     ,      i3,' AND',i3, /10X,61(1H*),//)
 IF (nout3t == 2) GO TO 891
 IF (iprtc  == 1) WRITE (log2,871)
 871 FORMAT (' NAME   CODE    DELTA P   ELEMENT',7X,  &
     'MESHPOINTS -  J   I',9X,'J   I',9X,'J   I',/)
 lnct  = lnct - 4
 ielem = 0
 xsign =-FLOAT(nsign)
 l4    = l2 + 1
 idata(1) = name1(1)
 idata(2) = name1(2)
 idata(3) = 60
 DO  j = 1,itub
   DO  i = l4,l3
     CALL alg03 (lnct,2)
     ielem = ielem + 1
     l5 = i - 1
     l6 = j + 1
     IF (i == l3) GO TO 880
     pload = xsign*((deltp(j,l5)+deltp(l6,l5)+deltp(l6,i))/3.0)
     IF (nblade(i) < 0) pload = xsign*((deltp(j,l5)+deltp(j,i) +  &
         deltp(l6,l5)+deltp(l6,i))*0.25)
     IF (iprtc == 1) WRITE (log2,900) pload,ielem,l6,l5,l6,i,j,l5
     rdata(4) = pload
     idata(5) = ielem
     CALL WRITE (iscr,idata,5,1)
     ielem = ielem + 1
     IF (nblade(i) >= 0) pload = xsign*((deltp(j,l5)+deltp(l6,i)+  &
         deltp(j,i))/3.0)
     IF (iprtc == 1) WRITE (log2,900) pload,ielem,j,l5,l6,i,j,i
     rdata(4) = pload
     idata(5) = ielem
     CALL WRITE (iscr,idata,5,1)
     CYCLE
     880 pload    = xsign*((deltp(j,l5)+deltp(l6,l5))/3.0)
     IF (nblade(i) < 0) pload = pload*0.75
     IF (iprtc == 1) WRITE (log2,900) pload,ielem,j,l5,l6,l5,l6,i
     rdata(4) = pload
     idata(5) = ielem
     CALL WRITE (iscr,idata,5,1)
     ielem    = ielem + 1
     IF (nblade(i) >= 0) pload = xsign*(deltp(j,l5)/3.0)
     IF (iprtc == 1) WRITE (log2,900) pload,ielem,j,l5,l6,i,j,i
     rdata(4) = pload
     idata(5) = ielem
     CALL WRITE (iscr,idata,5,1)
   END DO
 END DO
 900 FORMAT (' PLOAD2   60',f12.5,i7,14X,3(i10,i4))
 l1 = l3
 891 IF (nout3t == 1) GO TO 820
 
!     OUTPUT RELATIVE TOTAL TEMPERATURES AT NODES ON *TEMP* CARDS
 
 CALL alg03 (lnct,10)
 lnct = lnct - 6
 IF (iprtc == 1) WRITE (log2,892)
 892 FORMAT (//,' NAME   CODE    DELTA T   NODE',10X,'MESHPOINTS -  ',  &
     'J   I   COORDINATES -   RADIAL       AXIAL',/)
 inode = 1
 idata(1) = name2(1)
 idata(2) = name2(2)
 idata(3) = 70
 DO  j = 1,nstrms
   DO  i = l2,l3
     CALL alg03(lnct,1)
     idata(4) = inode
     rdata(5) = tr(j,i)
     CALL WRITE (iscr,idata,5,1)
     IF (iprtc == 1) WRITE (log2,912) tr(j,i),inode,j,i,r(j,i),x(j,i)
     inode = inode + 1
   END DO
 END DO
 912 FORMAT (' TEMP     70',f12.5,i6,21X,2I4,16X,f10.4,2X,f10.4)
 GO TO 820
 920 CONTINUE
 
!     PUNCH STREAML2 BULK DATA CARDS FOR EACH STREAMLINE
!     CHANGE THE SIGN ON THE STAGGER AND FLOW ANGLES FOR STREAML2 CARDS.
!     THIS CHANGE IS NECESSARY BECAUSE OF THE AERODYNAMIC PROGRAMS IN
!     NASTRAN MODULE AMG THAT USE THESE ANGLES.
 
 IF (ledgeb*itrleb == 0) GO TO 940
 IF (istrml == -1 .OR. istrml == 1) GO TO 940
 WRITE (log2,931)
 nstnsx = itrleb - ledgeb + 1
 DO  ileb = 1,nstrms
   radius = (rmdv(ileb,5)+rmdv(ileb,6))/2.0
   bspace = (6.283185*radius)/FLOAT(nbldes)
   stag(ileb  ) = -1.0*stag(ileb  )
   rmdv(ileb,4) = -1.0*rmdv(ileb,4)
   WRITE (lpunch,932) ileb,nstnsx,stag(ileb),chordd(ileb),radius,  &
       bspace,rmdv(ileb,1),rmdv(ileb,2),ileb,ileb,rmdv(ileb,3), rmdv(ileb,4)
   WRITE (log2,933) ileb,nstnsx,stag(ileb),chordd(ileb),radius,  &
       bspace,rmdv(ileb,1),rmdv(ileb,2),rmdv(ileb,3),rmdv(ileb,4)
 END DO
 931 FORMAT (//10X,47HNASTRAN - streaml2 - compressor blade bulk DATA,  &
     /10X,49(1H*), /,'  SLN  NSTNS  STAGGER    CHORD    RADIUS',  &
     '    BSPACE     MACH       DEN       VEL      FLOWA',/)
 932 FORMAT (8HSTREAML2,2I8,f8.3,3F8.5,2F8.6,5H+strl,i2,5H+strl,i2,  &
     f8.1,f8.3 )
 933 FORMAT (i5,i6,2X,f8.3,3(2X,f8.5),2(2X,f8.6),2X,f8.1,2X,f8.3)
 940 CONTINUE
 RETURN
END SUBROUTINE alg11
