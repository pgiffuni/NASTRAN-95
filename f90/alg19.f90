SUBROUTINE alg19(log1,log2,log3,log5,nlines,nspec,kpts,rsta, &
                 xsta,r,zr,b1,b2,tc,pi,c1,nblade,ccord,block, &
                 alpb,epslon,ifangs,ipunch,naero)
  
 INTEGER, INTENT(IN OUT)                  :: log1
 INTEGER, INTENT(IN OUT)                  :: log2
 INTEGER, INTENT(IN)                      :: log3
 INTEGER, INTENT(IN)                      :: log5
 INTEGER, INTENT(IN)                      :: nlines
 INTEGER, INTENT(IN OUT)                  :: nspec
 INTEGER, INTENT(IN OUT)                  :: kpts(1)
 REAL, INTENT(IN OUT)                     :: rsta(21,10)
 REAL, INTENT(IN OUT)                     :: xsta(21,10)
 REAL, INTENT(IN)                         :: r(10,21)
 REAL, INTENT(IN OUT)                     :: zr(1)
 REAL, INTENT(IN OUT)                     :: b1(1)
 REAL, INTENT(IN OUT)                     :: b2(1)
 REAL, INTENT(IN OUT)                     :: tc(1)
 REAL, INTENT(IN OUT)                     :: pi
 REAL, INTENT(IN OUT)                     :: c1
 INTEGER, INTENT(IN OUT)                  :: nblade
 REAL, INTENT(IN)                         :: ccord(1)
 REAL, INTENT(IN OUT)                     :: block(10,21)
 REAL, INTENT(IN OUT)                     :: alpb(10,21)
 REAL, INTENT(IN OUT)                     :: epslon(10,21)
 INTEGER, INTENT(IN OUT)                  :: ifangs(1)
 INTEGER, INTENT(IN)                      :: ipunch
 INTEGER, INTENT(IN OUT)                  :: naero

 DIMENSION       idata(24),rdata(6),   &
      nr(10),nterp(10),nmach(10),nloss(10),nl1(10),  &
     nl2(10),neval(10),ncurve(10),nliter(10),ndel(10),  &
     rr(21,10),xloss(21,10),rte(5),dm(11,5),  &
     dvfrac(11,5),rdte(21),deltad(21),ac(21),f137b(8),   &
     f137s(5),f142tc(7),f161d(8,5),f195m(8,2),f164xb(8), &
     f172k(7),nout1(10),nout2(10),sol(21),dev(10,21),  &
     devv(10,5),dx(10),x(10),doo(5), nout3(10),nblad(10)

 COMMON /ud3prt/ iprtc

 DATA    f137b / 0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0/
 DATA    f137s / 0.4,0.8,1.2,1.6,2.0/
 DATA    f142tc/ 0.0,0.02,0.04,0.06,0.08,0.10,0.12/
 DATA    f161d / 0.0,0.009,0.17,0.29,0.42,0.59,0.79,1.05,0.0,0.12,  &
     0.30,0.51,0.75,1.05,1.47,2.07,0.0,0.16,0.33,0.61,  &
     0.95,1.42,2.12,3.07,0.0,0.17,0.40,0.72,1.11,1.71,  &
     2.62,3.95,0.0,0.2,0.44,0.78,1.21,1.90,3.01,4.75/
 DATA    f195m / 0.17,0.173,0.179,0.189,0.206,0.232,0.269,0.310,  &
     0.25,0.255,0.261,0.268,0.278,0.292,0.312,0.342 /
 DATA    f164xb/ 0.965,0.945,0.921,0.890,0.850,0.782,0.679,0.550/
 DATA    f172k / 0.0,0.160,0.331,0.521,0.74,1.0,1.300/
 
 lmax   = 60
 CALL fread (log1,idata,6,1)
 nrad   = idata(1)
 ndpts  = idata(2)
 ndatr  = idata(3)
 nswitc = idata(4)
 nle    = idata(5)
 nte    = idata(6)
 CALL fread (log1,rdata,2,1)
 xkshpe = rdata(1)
 speed  = rdata(2)
 CALL fread (log1,idata,3,1)
 nout1(nle) = idata(1)
 nout2(nle) = idata(2)
 nout3(nle) = idata(3)
 IF (iprtc == 1) WRITE (log2,20) nrad,ndpts,ndatr,nswitc,nle,nte,  &
     xkshpe,speed,nle,nout1(nle),nout2(nle),nout3(nle)
 20 FORMAT (1H1,9X,'DATA INTERFACING ROUTINE - DEVIATION CALCULATIONS'  &
     ,      ' AND DATA FORMATTING', /10X,69(1H*), /10X,5HINPUT, /10X,  &
     5(1H*), //10X,6HNRAD =,i3,9H  ndpts =,i3,9H  ndatr =,i3,  &
     11H  nswitch =,i2,7H  nle =,i2,7H  nte =,i3, //10X,  &
     8HXKSHPE =,f7.4,9H  speed =,f9.1, //10X,'AT LEADING EDGE ',  &
     '(STATION,I3,9H) NOUT1 =',i2,9H  nout2 =,i2,9H  nout3 =,i2)
 lnct = 10
 k    = nle + 1
 DO  i = k,nte
   CALL fread (log1,idata,14,1)
   nr(i)     = idata(1)
   nterp(i)  = idata(2)
   nmach(i)  = idata(3)
   nloss(i)  = idata(4)
   nl1(i)    = idata(5)
   nl2(i)    = idata(6)
   neval(i)  = idata(7)
   ncurve(i) = idata(8)
   nliter(i) = idata(9)
   ndel(i)   = idata(10)
   nout1(i)  = idata(11)
   nout2(i)  = idata(12)
   nout3(i)  = idata(13)
   nblad(i)  = idata(14)
   IF (lnct+6+nr(i) <= lmax) GO TO 50
   IF (iprtc /= 0) WRITE (log2,40)
   40 FORMAT (1H1)
   lnct = 1
   50 lnct = lnct + 6 + nr(i)
   IF (iprtc == 1) WRITE (log2,60) i,nr(i),nterp(i),nmach(i),  &
       nloss(i),nl1(i),nl2(i),neval(i),ncurve(i),nliter(i),  &
       ndel(i),nout1(i),nout2(i),nout3(i),nblad(i)
   60 FORMAT (/10X,7HSTATION,i3,7H   nr =,i3,9H  nterp =,i2,9H  nmach =,  &
       i2,9H  nloss =,i2,7H  nl1 =,i3,7H  nl2 =,i3,9H  neval =,i2,  &
       8HNCURVE =,i2,10H  nliter =,i3,8H  ndel =,i2, /22X,  &
       7HNOUT1 =,i2,9H  nout2 =,i2,9H  nout3 =,i2,9H  nblad =,i3)
   l1 = nr(i)
   DO  j = 1,l1
     CALL fread (log1,rdata,2,1)
     rr(j,i)    = rdata(1)
     xloss(j,i) = rdata(2)
   END DO
   IF (iprtc == 1) WRITE (log2,90) (rr(j,i),xloss(j,i),j=1,l1)
 END DO
 90 FORMAT (/14X,6HRADIUS,6X,15HLOSS descriptor,//,(f20.4,f17.6))
 IF (lnct+7+ndpts <= lmax) GO TO 100
 IF (iprtc /= 0) WRITE (log2,40)
 lnct = 1
 100 lnct = lnct + 2
 IF (iprtc == 1) WRITE (log2,110) nrad
 110 FORMAT (/10X,28HDEVIATION fraction curves at,i2,6H radii)
 DO  k = 1,nrad
   IF (lnct+5+ndpts <= lmax) GO TO 120
   IF (iprtc /= 0) WRITE (log2,40)
   lnct = 1
   120 lnct = lnct + 5 + ndpts
   CALL fread (log1,rte(k),1,1)
   DO  j = 1,ndpts
     CALL fread (log1,rdata,2,1)
     dm(j,k)     = rdata(1)
     dvfrac(j,k) = rdata(2)
   END DO
   IF (iprtc == 1) WRITE (log2,150) rte(k),(dm(j,k),dvfrac(j,k), j=1,ndpts)
 END DO
 150 FORMAT (/10X,5HRTE =,f8.4, //15X,2HDM,10X,6HDVFRAC, //, (f20.5,f13.5))
 DO  j = 1,ndatr
   CALL fread (log1,rdata,3,1)
   rdte(j)   = rdata(1)
   deltad(j) = rdata(2)
   ac(j)     = rdata(3)
 END DO
 IF (lnct+3+ndatr <= lmax) GO TO 170
 IF (iprtc /= 0) WRITE (log2,40)
 lnct = 1
 170 lnct = lnct + 3 + ndatr
 IF (iprtc == 1) WRITE (log2,180) (rdte(j),deltad(j),ac(j), j=1,ndatr)
 180 FORMAT (/15X,4HRDTE,6X,6HDELTAD,9X,2HAC,//,(f20.4,f11.3,f13.4))
 IF (lnct+6+nlines <= lmax) GO TO 190
 IF (iprtc /= 0) WRITE (log2,40)
 lnct = 1
 190 lnct = lnct + 6 + nlines
 IF (iprtc == 1) WRITE (log2,200)
 200 FORMAT (/10X,7HRESULTS, /,10X,7(1H*))
 IF (iprtc == 1) WRITE (log2,210)
 210 FORMAT (/5X,10HSTREAMLINE,5X,5HBETA1,6X,5HBETA2,5X,6HCAMBER,7X,  &
     3HT/c,8X,3HA/c,6X,8HSOLIDITY,4X,11HADDIT. devn,4X, 15HTOTAL deviation,/)
 DO  j = 1,nlines
   xj = j
   CALL alg15 (zr,b1,nspec,xj,beta1,1,0)
   CALL alg15 (zr,b2,nspec,xj,beta2,1,0)
   CALL alg15 (zr,tc,nspec,xj,thick,1,0)
   q = 1.0
   IF (speed > 0.0) q = -1.0
   camber = (beta1-beta2)*q
   solid  = ccord(j)*FLOAT(nblade)/(pi*(r(nle,j)+r(nte,j)))
   beta1  = beta1*q
   CALL alg15 (f137b,f195m(1,nswitc),8,beta1,xms,1,0)
   CALL alg15 (f137b,f164xb,8,beta1,xb,1,0)
   CALL alg15 (f142tc,f172k,7,thick,xkdt,1,0)
     DO 220 k = 1,5
 220 CALL alg15 (f137b,f161d(1,k),8,beta1,doo(K),1,0)
   CALL alg15 (f137s,doo,5,solid,do,1,1)
   CALL alg15 (rdte,deltad,ndatr,r(nte,j),dadd,1,0)
   CALL alg15 (rdte,ac,ndatr,r(nte,j),aonc,1,0)
   sol(j) = solid
   dev(nte,j) = (dadd+do*xkshpe*xkdt+camber*solid**(-xb)*  &
       (xms+0.5*(aonc-0.5)))*q
   beta2  = beta2*q
   DO  i = nle,nte
     CALL alg15 (rsta(1,i),xsta(1,i),kpts(i),r(i,j),x(i),1,0)
   END DO
   dx(nle) = 0.0
   k = nle + 1
   DO  i = k,nte
     dx(i) = dx(i-1) + SQRT((x(i)-x(i-1))**2+(r(i,j)-r(i-1,j))**2)
   END DO
   x1 = dx(nte)
   DO  i = k,nte
     dx(i) = dx(i)/x1
   END DO
   l2 = nte - nle - 1
   k  = nle + 1
   DO  l1 = 1,nrad
     CALL alg15 (dm(1,l1),dvfrac(1,l1),ndpts,dx(k),devv(k,l1),l2,0)
   END DO
   kk = nte - 1
   DO  i  = k,kk
     DO  l1 = 1,nrad
       doo(l1) = devv(i,l1)
     END DO
     CALL alg15 (rte,doo,nrad,r(nte,j),devfr,1,0)
     dev(i,j) = dev(nte,j)*devfr
   END DO
   IF (iprtc == 1) WRITE (log2,300) j,beta1,beta2,camber,thick,  &
       aonc,solid,dadd,dev(nte,j)
 END DO
 300 FORMAT (i11,f14.3,2F11.3,2F11.4,f12.5,f14.4,f17.4)
 IF (ifangs(nle) == 0) GO TO 340
 IF (naero == 0) GO TO 330
 idata(1)  = nlines
 idata(2)  = 0
 idata(3)  = 0
 idata(4)  = 0
 idata(5)  = 0
 idata(6)  = 0
 idata(7)  = 0
 idata(8)  = 0
 idata(9)  = 0
 idata(10) = 0
 idata(11) = 0
 idata(12) = 0
 idata(13) = nout1(nle)
 idata(14) = nout2(nle)
 idata(15) = nout3(nle)
 idata(16) = 0
 CALL WRITE (log5,idata,16,1)
 CALL WRITE (log5,0.0,1,1)
 310 FORMAT (i3,11(2X,1H0),3I3,3H  0, /,4H 0.0)
 DO  j = 1,nlines
   rdata(1) = r(nle,j)
   rdata(2) = alpb(nle,j)
   rdata(3) = 0.0
   rdata(4) = epslon(nle,j)
   rdata(5) = 0.0
   rdata(6) = 0.0
   CALL WRITE (log5,rdata,6,1)
   rdata(1) = 0.0
   rdata(2) = 0.0
   rdata(3) = 0.0
   rdata(4) = 0.0
   CALL WRITE (log5,rdata,4,1)
 END DO
 320 FORMAT (2F12.7,12X,f12.7,24X,/,4H 0.0,44X)
 IF (ipunch == 0) GO TO 340
 330 WRITE (log3,310) nlines,nout1(nle),nout2(nle),nout3(nle)
 WRITE (log3,320) (r(nle,j),alpb(nle,j),epslon(nle,j),j=1,nlines)
 340 DO  i = k,nte
   DO  j = 1,nlines
     rdte(j) = r(i,j)
   END DO
   CALL alg15 (rr(1,i),xloss(1,i),nr(i),rdte,deltad,nlines,0)
   nx = log5
   IF (naero == 0) nx = log3
   IF (nx == log3) GO TO 360
   idata(1)  = nlines
   idata(2)  = nterp(i)
   idata(3)  = 0
   idata(4)  = nmach(i)
   idata(5)  = 6
   idata(6)  = nloss(i)
   idata(7)  = nl1(i)
   idata(8)  = nl2(i)
   idata(9)  = neval(i)
   idata(10) = ncurve(i)
   idata(11) = nliter(i)
   idata(12) = ndel(i)
   idata(13) = nout1(i)
   idata(14) = nout2(i)
   idata(15) = nout3(i)
   idata(16) = nblad(i)
   CALL WRITE (log5,idata,16,1)
   CALL WRITE (log5,speed,1,1)
   DO  j = 1,nlines
     rdata(1) = r(i,j)
     rdata(2) = alpb(i,j)
     rdata(3) = deltad(j)
     rdata(4) = epslon(i,j)
     rdata(5) = BLOCK(i,j)
     rdata(6) = sol(j)
     CALL WRITE (log5,rdata,6,1)
     rdata(1) = dev(i,j)
     rdata(2) = 0.0
     rdata(3) = 0.0
     rdata(4) = 0.0
     CALL WRITE (log5,rdata,4,1)
   END DO
   GO TO 365
   360 WRITE (nx,380) nlines,nterp(i),nmach(i),nloss(i),nl1(i),nl2(i),  &
       neval(i),ncurve(i),nliter(i),ndel(i),nout1(i),nout2(i),  &
       nout3(i),nblad(i),speed,(r(i,j),alpb(i,j),deltad(j),  &
       epslon(i,j),BLOCK(i,j),sol(j),dev(i,j),j=1,nlines)
   365 IF (nx == log3) CYCLE
   nx = log3
   IF (naero /= 0 .AND. ipunch /= 0) GO TO 360
 END DO
 380 FORMAT (2I3,3H  0,i3,3H  6,11I3, /,f12.3,/,(6F12.7, /,f12.7,36X))
 
 RETURN
END SUBROUTINE alg19
