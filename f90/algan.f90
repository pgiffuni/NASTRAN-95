SUBROUTINE algan
     
 REAL :: ix,iy,ixy,ipx,ipy,ixn,iyn,ixyn,ixd,iyd
 DIMENSION       rles(21),tcs(21),tes(21),zzs(21),perspt(21),  &
     alpb(10,21),BLOCK(10,21),epslon(10,21),ccord(21),  &
     ys(21,70),yp(21,70),xp(21,70),xs(21,70),zs(21,70),  &
     zq(21),tq(21),blafor(10,21),xcamb(21,10),  &
     tharr(21,10),title(18),idata(24),rdata(6)
 COMMON /system/ ksystm(90),lpunch
 COMMON /ud3prt/ iprtc,istrml,ipgeom
 COMMON /udstr2/ nbldes,stag(21),chordd(21)
 COMMON /contrl/ nanal,naero,narbit,log1,log2,log3,log4,log5,log6
 COMMON /ud3anc/ epz(80,4),r(10,21),zout(21),ss(100),x(100),  &
     yprime(100),ysemi(21,31),xsemi(21,31),zp(21,70),  &
     zsemi(21,31),title2(18),xhere(10),xtemp(100),  &
     rad(100),temp1(21),temp2(21),temp3(21),temp4(21),  &
     zr(21),b1(21),b2(21),pp(21),qq(21),zz(21),rle(21),  &
     tc(21),te(21),cord(21),delx(21),dely(21),s(21),  &
     bs(21),xsemj(21,31),ysemj(21,31),zsemj(21,31),  &
     xsta(21,10),rsta(21,10),kpts(21),sigma(100),  &
     tanphi(10,21),zcamb(21,10),ycamb(21,10), ifangs(10),theta(21,10),alpha(21,10)
 EQUIVALENCE     (title(1),title2(1))
 
 pi = 4.0*ATAN(1.0)
 c1 = 180.0/pi
 CALL fread (log1,title2,18,1)
 IF (iprtc == 1) WRITE (log2,120) title2
 120  FORMAT (1H1,31X,'PROGRAM ALG - COMPRESSOR DESIGN - ANALYTIC MEAN',  &
     'LINE BLADE SECTION', /32X,65(1H*), //10X,5HTITLE,25X,1H=, 18A4)
 CALL fread (log1,idata,17,1)
 nlines = idata(1)
 nstns  = idata(2)
 nz     = idata(3)
 nspec  = idata(4)
 npoint = idata(5)
 nblade = idata(6)
 istak  = idata(7)
 ipunch = idata(8)
 isecn  = idata(9)
 ifcord = idata(10)
 ifplot = idata(11)
 iprint = idata(12)
 isplit = idata(13)
 inast  = idata(14)
 irle   = idata(15)
 irte   = idata(16)
 nsign  = idata(17)
 nbldes = nblade
 IF (iprtc == 0) iprint = 3
 IF (inast == 0 .AND. ipgeom /= -1) inast = -4
 IF (iprtc == 1) WRITE (log2,140) nlines,nstns,nz,nspec,npoint,  &
     nblade,istak,ipunch,isecn,ifcord,ifplot,iprint,isplit,inast, irle,irte,nsign
 140  FORMAT (10X,24HNUMBER of streamsurfaces,6X,1H=,i3, /10X,18HNUMBER &
     &of stations,12X,1H=,i3, /10X,27HNUMBER of constant-z planes,3X,1H= &
     ,i3, /10X,27HNUMBER of blade DATA points,3X,1H=,i3, /10X,31HNUMBER&
     &of points on surfaces  =,i3, /10X,29HNUMBER of blades in blade row,&
     1X,1H=,i3, /10X,5HISTAK,25X,1H=,i3, /10X,6HIPUNCH,24X,1H=,i3,   &
     /10X,5HISECN,25X,1H=,i3,/,10X,6HIFCORD,24X,1H=,i3, /10X,6HIFPLOT,24X  &
     ,1H=,i3, /10X,6HIPRINT,24X,1H=,i3, /10X,6HISPLIT,24X,1H=,i3, /10X,  &
     5HINAST,25X,1H=,i3, /10X,4HIRLE,26X,1H=,i3, /10X,4HIRTE,26X,1H=,i3  &
     , /10X,5HNSIGN,25X,1H=,i3)
 CALL fread (log1,rdata,5,1)
 zinner = rdata(1)
 zouter = rdata(2)
 scale  = rdata(3)
 stackx = rdata(4)
 pltsze = rdata(5)
 IF (iprtc == 1) WRITE (log2,160) zinner,zouter,scale,stackx, pltsze
 160  FORMAT (/10X,6HZINNER,24X,1H=,f8.4, /10X,6HZOUTER,24X,1H=,f8.4, /  &
     10X,5HSCALE,25X,1H=,f8.4, /10X,6HSTACKX,24X,1H=,f8.4, /10X,6HPLTSZE, &
     24X,1H=,f8.4, //20X,36HSTREAMSURFACE geometry specification)
 lnct = 30
 DO  i = 1,nstns
   CALL fread (log1,idata,2,1)
   kpts(i)   = idata(1)
   ifangs(i) = idata(2)
   kpt = kpts(i)
   DO  k = 1,kpt
     CALL fread (log1,rdata,2,1)
     xsta(k,i) = rdata(1)
     rsta(k,i) = rdata(2)
   END DO
   IF (kpts(i) >= 2) GO TO 170
   kpts(i) = 2
   xsta(2,i) = xsta(1,i)
   rsta(2,i) = rsta(1,i) + 1.0
   170  DO  j = 1,nlines
     CALL fread (log1,rdata,2,1)
     r(i,j)      = rdata(1)
     blafor(i,j) = rdata(2)
   END DO
   idum = kpts(i)
   IF (nlines > idum) idum = nlines
   IF (lnct <= 54-idum) GO TO 210
   IF (iprtc /= 0) WRITE (log2,200)
   200  FORMAT (1H1)
   lnct = 2
   210  lnct = lnct + idum + 7
   IF (inast /= 0) GO TO 240
   IF (iprtc == 1)  WRITE (log2,220) i,kpts(i),i,ifangs(i)
   220  FORMAT (/10X,'COMPUTING STATION',i3,5X,'NUMBER OF DESCRIBING ',  &
       'POINTS=',i3,6X,7HIFANGS(,i2,2H)=,i3, //6X,'DESCRIPTION',  &
       9X,'STREAMLINE',5X,5HRADII,/6X,1HX,9X,1HR,11X,6HNUMBER,//)
   DO  k = 1,idum
     IF (iprtc == 1 .AND. k <= kpts(i).AND.k <= nlines)  &
         WRITE (log2,260) xsta(k,i),rsta(k,i),k,r(i,k)
     IF (iprtc == 1 .AND. k <= kpts(i).AND.k > nlines)  &
         WRITE (log2,270) xsta(k,i),rsta(k,i)
     IF (iprtc == 1 .AND. k > kpts(i).AND.k <= nlines)  &
         WRITE (log2,280) k,r(i,k)
   END DO
   IF (inast == 0) CYCLE
   240  IF (iprtc == 1) WRITE (log2,290) i,kpts(i),i,ifangs(i)
   DO  k = 1,idum
     IF (iprtc == 1 .AND. k <= kpts(i).AND.k <= nlines)  &
         WRITE (log2,260) xsta(k,i),rsta(k,i),k,r(i,k),blafor(i,k)
     IF (iprtc == 1 .AND. k <= kpts(i).AND.k > nlines)  &
         WRITE (log2,270) xsta(k,i),rsta(k,i)
     IF (iprtc == 1 .AND. k > kpts(i).AND.k <= nlines)  &
         WRITE (log2,280) k,r(i,k),blafor(i,k)
   END DO
   260  FORMAT (3X,f8.4,2X,f8.4,8X,i2,9X,f8.4,9X,f8.4)
   270  FORMAT (3X,f8.4,2X,f8.4)
   280  FORMAT (29X,i2,9X,f8.4,9X,f8.4)
   290  FORMAT (/10X,'COMPUTING STATION',i3,5X,'NUMBER OF DESCRIBING ',  &
       'POINTS=',i3,6X,7HIFANGS(,i2,2H)=,i3, //6X,'DESCRIPTION',  &
       9X,'STREAMLINE',5X,5HRADII,9X,'DELTA PRESSURE', /6X,1HX,  &
       9X,1HR,11X,6HNUMBER, //)
 END DO
 sq = 0.0
 sb = 0.0
 IF (isecn == 1 .OR. isecn == 3) GO TO 340
 DO  isbs = 1,nspec
   s(isbs)  = 0.0
   bs(isbs) = 0.0
 END DO
 IF (lnct <= 54-nspec) GO TO 310
 IF (iprtc /= 0) WRITE (log2,200)
 lnct = 1
 310  lnct = lnct + nspec + 6
 DO  j = 1,nspec
   CALL fread (log1,rdata,6,1)
   zr(j)  = rdata(1)
   b1(j)  = rdata(2)
   b2(j)  = rdata(3)
   pp(j)  = rdata(4)
   qq(j)  = rdata(5)
   rle(j) = rdata(6)
   CALL fread (log1,rdata,6,1)
   tc(j)   = rdata(1)
   te(j)   = rdata(2)
   zz(j)   = rdata(3)
   cord(j) = rdata(4)
   delx(j) = rdata(5)
   dely(j) = rdata(6)
 END DO
 IF (iprtc == 1) WRITE (log2,330) (zr(j),b1(j),b2(j),pp(j),qq(j),  &
     rle(j),tc(j),te(j),zz(j),cord(j),delx(j),dely(j),j=1,nspec)
 330  FORMAT (/20X,'SECTION GEOMETRY SPECIFICATION', //10X,'STREAMLINE',  &
     '  INLET',5X,6HOUTLET,4X,6HY2 LE/,4X,6HY2 te/,3X,48HLE radi&
     &us max thick te thick  point of  chord OR,3X,7HX stack,3X,7HY stac&
     &k, /11X,6HNUMBER,5X,5HANGLE,5X,5HANGLE,3X,19HMAX value MAX value,3X  &
     ,6H/chord,4X,6H/chord,3X,8H/2*chord,2X,18HMAX thick axial cd,4X,&
     6Hoffset,4X,6HOFFSET, //,(10X,f7.2,3X,f8.3,f10.3,2F10.4,3F10.5,  &
     2F10.4,f11.6,f10.6))
 GO TO 390
 340  IF (lnct <= 50-2*nspec) GO TO 350
 IF (iprtc /= 0) WRITE (log2,200)
 lnct = 1
 350  lnct = lnct + 10 + 2*nspec
 DO  j = 1,nspec
   CALL fread (log1,rdata,6,1)
   zr(j)  = rdata(1)
   b1(j)  = rdata(2)
   b2(j)  = rdata(3)
   pp(j)  = rdata(4)
   qq(j)  = rdata(5)
   rle(j) = rdata(6)
   CALL fread (log1,rdata,6,1)
   tc(j)   = rdata(1)
   te(j)   = rdata(2)
   zz(j)   = rdata(3)
   cord(j) = rdata(4)
   delx(j) = rdata(5)
   dely(j) = rdata(6)
   CALL fread (log1,rdata,2,1)
   s(j)  = rdata(1)
   bs(j) = rdata(2)
 END DO
 IF (iprtc == 1) WRITE (log2,330) (zr(j),b1(j),b2(j),pp(j),qq(j),  &
     rle(j),tc(j),te(j),zz(j),cord(j),delx(j),dely(j),j=1,nspec)
 IF (iprtc == 1 .AND. isecn == 1) WRITE (log2,370) (zr(j),s(j),  &
     bs(j),j=1,nspec)
 370  FORMAT (/10X,'STREAMLINE  INFLECTION  INFLECTION', /11X,'NUMBER',  &
     8X,5HPOINT,7X,5HANGLE, //,(10X,f7.2,f14.5,f11.3))
 IF (iprtc == 1 .AND. isecn == 3) WRITE (log2,380) (zr(j),s(j),  &
     bs(j),j=1,nspec)
 380  FORMAT (/10X,'STREAMLINE  TRANSITION  DEL ANGLE', /11X,'NUMBER',  &
     8X,5HPOINT,6X,7HFROM LE, //,(10X,f7.2,f14.5,f11.3))
 390  IF (isplit == 0) GO TO 430
 DO  j = 1,nspec
   CALL fread (log1,rdata,5,1)
   rles(j)   = rdata(1)
   tcs(j)    = rdata(2)
   tes(j)    = rdata(3)
   zzs(j)    = rdata(4)
   perspt(j) = rdata(5)
 END DO
 IF (iprtc == 1) WRITE (log2,410)
 410  FORMAT (/20X,13HSPLITTER DATA, //10X,10HSTREAMLINE,2X,47HLE radius &
     &max thick te thick  point of per cent, /11X,6HNUMBER,7X,6H/chord,  &
     4X,6H/chord,3X,8H/2*chord,2X,9HMAX thick,2X,8HSPLITTER, /)
 IF (iprtc == 1) WRITE (log2,420) (zr(j),rles(j),tcs(j),tes(j),  &
     zzs(j),perspt(j),j=1,nspec)
 420  FORMAT (10X,f7.2,3X,f8.3,f10.3,3F10.4)
 430  CONTINUE
 IF (ifplot == 0 .OR. ifplot == 4) GO TO 440
 ikdum = 0
 IF (b1(1) < 0.0) ikdum = 1
 IF (ifplot == 1 .OR. ifplot == 3) CALL alg17 (istak,pltsze,1,  &
     title,ikdum,ifplot)
 440  ndum  = npoint
 iidum = isecn
 DO  j = 1,nlines
   npoint = ndum
   isecn  = iidum
   DO  i = 1,nstns
     kpt = kpts(i)
     CALL alg15 (rsta(1,i),xsta(1,i),kpt,r(i,j),xhere(i),1,0)
   END DO
   x(1)   = xhere(1)
   x(100) = xhere(nstns)
   ax = (x(100)-x(1))/99.0
   DO  i = 2,99
     x(i) = x(i-1) + ax
   END DO
   CALL alg14 (xhere,r(1,j),nstns,x,xdum,yprime,100,1)
   CALL alg14 (xhere,r(1,j),nstns,xhere,xdum,tanphi(1,j),nstns,1)
   ss(1) = 0.0
   DO  i = 2,100
     ss(i) = ss(i-1) + ax*SQRT(1.0+((yprime(i)+yprime(i-1))/2.0)**2)
   END DO
   xj = j
   CALL alg15 (zr,b1,nspec,xj,beta1,1,0)
   CALL alg15 (zr,b2,nspec,xj,beta2,1,0)
   CALL alg15 (zr,pp,nspec,xj,p,1,0)
   CALL alg15 (zr,qq,nspec,xj,q,1,0)
   CALL alg15 (zr,rle,nspec,xj,yzero,1,0)
   CALL alg15 (zr,tc,nspec,xj,t,1,0)
   CALL alg15 (zr,te,nspec,xj,yone,1,0)
   CALL alg15 (zr,delx,nspec,xj,xdel,1,0)
   CALL alg15 (zr,dely,nspec,xj,ydel,1,0)
   CALL alg15 (zr,zz,nspec,xj,z,1,0)
   CALL alg15 (zr,cord,nspec,xj,chd,1,0)
   IF (isecn == 0 .OR. isecn == 2) GO TO 480
   CALL alg15 (zr,s,nspec,xj,sq,1,0)
   CALL alg15 (zr,bs,nspec,xj,sb,1,0)
   480  IF (isplit == 0) GO TO 490
   CALL alg15 (zr,rles,nspec,xj,yzeros,1,1)
   CALL alg15 (zr,tcs,nspec,xj,ts,1,1)
   CALL alg15 (zr,tes,nspec,xj,yones,1,1)
   CALL alg15 (zr,zzs,nspec,xj,zspmxt,1,1)
   CALL alg15 (zr,perspt,nspec,xj,perspj,1,1)
   490  CALL alg15 (x,ss,100,stackx,bx,1,1)
   CALL alg13 (j,ys,yp,xs,xp,ysemi,xsemi,log1,log2,npoint,iprint,  &
       beta1,beta2,p,q,yzero,t,yone,xdel,ydel,z,chd,lnct,ifcord,sq,  &
       sb,isecn,xsemj,ysemj,istak,xhere,x,ss,nstns,r,xtemp,yprime,  &
       rad,epz,bx,sigma,ccord,isplit,yzeros,ts,yones,zspmxt,perspj,  &
       inast,irle,irte,tharr)
   CALL alg15 (x,ss,100,stackx,bx,1,1)
   DO  i = 1,100
     x(i)  = x(i) - stackx
     ss(i) = ss(i)- bx
   END DO
   DO  i = 1,nstns
     xhere(i) = xhere(i) - stackx
   END DO
   IF (ifplot == 0 .OR. ifplot == 2 .OR. ifplot == 4) GO TO 570
   xplot = xs(j,1)*scale
   yplot = ys(j,1)*scale
   DO  i = 2,npoint
     xplot = xs(j,i)*scale
     yplot = ys(j,i)*scale
   END DO
   IF (isecn /= 2) GO TO 540
   DO  i = 2,30
     xplot = xsemj(j,i)*scale
     yplot = ysemj(j,i)*scale
   END DO
   540  DO  ii = 1,npoint
     i = npoint - ii + 1
     xplot = xp(j,i)*scale
     yplot = yp(j,i)*scale
   END DO
   DO  i = 2,30
     xplot = xsemi(j,i)*scale
     yplot = ysemi(j,i)*scale
   END DO
   xplot = xs(j,1)*scale
   yplot = ys(j,1)*scale
   570  ijdum = 0
   DO  i = 1,nstns
     IF (ifangs(i) >= 1) ijdum = 1
   END DO
   IF (ijdum == 0 .AND. inast == 0) GO TO 600
   CALL alg15 (ss,x,100,xtemp,xtemp,100,1)
   DO  i = 1,nstns
     CALL alg15 (xtemp,sigma,100,xhere(i),theta(j,i),1,1)
     CALL alg15 (xtemp,yprime,100,xhere(i),alpha(j,i),1,1)
     zcamb(j,i) = r(i,j)*COS(theta(j,i))
     xcamb(j,i) = xhere(i)
     ycamb(j,i) = r(i,j)*SIN(theta(j,i))
   END DO
   600  DO  i = 1,npoint
     xtemp(i) = xs(j,i)
   END DO
   CALL alg15 (ss,x,100,xtemp,xtemp,npoint,1)
   CALL alg15 (xhere,r(1,j),nstns,xtemp,rad,npoint,0)
   k = 1
   DO  i = 1,npoint
     eps = epz(i,k)
     zs(j,i) = rad(i)*COS(eps)
     ys(j,i) = rad(i)*SIN(eps)
     xs(j,i) = xtemp(i)
   END DO
   DO  i = 1,npoint
     xtemp(i) = xp(j,i)
   END DO
   CALL alg15 (ss,x,100,xtemp,xtemp,npoint,1)
   CALL alg15 (xhere,r(1,j),nstns,xtemp,rad,npoint,0)
   k = 2
   DO  i = 1,npoint
     eps = epz(i,k)
     zp(j,i) = rad(i)*COS(eps)
     yp(j,i) = rad(i)*SIN(eps)
     xp(j,i) = xtemp(i)
   END DO
   DO  i = 1,31
     xtemp(i) = xsemi(j,i)
   END DO
   CALL alg15 (ss,x,100,xtemp,xtemp,31,1)
   CALL alg15 (xhere,r(1,j),nstns,xtemp,rad,31,0)
   k = 3
   DO  i = 1,31
     eps = epz(i,k)
     zsemi(j,i) = rad(i)*COS(eps)
     ysemi(j,i) = rad(i)*SIN(eps)
     xsemi(j,i) = xtemp(i)
   END DO
   IF (isecn /= 2) GO TO 690
   DO  i = 1,31
     xtemp(i) = xsemj(j,i)
   END DO
   CALL alg15 (ss,x,100,xtemp,xtemp,31,1)
   CALL alg15 (xhere,r(1,j),nstns,xtemp,rad,31,0)
   k = 4
   DO  i = 1,31
     eps = epz(i,k)
     zsemj(j,i) = rad(i)*COS(eps)
     ysemj(j,i) = rad(i)*SIN(eps)
     xsemj(j,i) = xtemp(i)
   END DO
   690  IF (iprint >= 2) CYCLE
   IF (lnct  <= 50) GO TO 700
   IF (iprtc /=  1) WRITE (log2,200)
   lnct = 1
   700  lnct = lnct + 5
   IF (iprtc == 1) WRITE (log2,710) j
   710  FORMAT (/10X,38HCARTESIAN coordinates on streamsurface,i3, //10X,  &
       8HPOINT no,5X,2HZ1,12X,2HX1,12X,2HY1,16X,2HZ2,12X,2HX2, 12X,2HY2, /)
   i = 1
   720  IF (iprtc == 1) WRITE (log2,730) i,zs(j,i),xs(j,i),ys(j,i),  &
       zp(j,i),xp(j,i),yp(j,i)
   730  FORMAT (10X,i5,3X,1P,3E14.5,4X,1P,3E14.5)
   i = i + 1
   lnct = lnct + 1
   IF (i > npoint) GO TO 750
   IF (lnct <=  59) GO TO 720
   IF (iprtc ==  1) WRITE (log2,740)
   740  FORMAT (1H1,9X,8HPOINT no,5X,2HZ1,12X,2HX1,12X,2HY1,16X,2HZ2,12X,  &
       2HX2,12X,2HY2, /)
   lnct = 2
   GO TO 720
   750  IF (lnct <= 50) GO TO 760
   IF (iprtc /= 0) WRITE (log2,200)
   lnct = 1
   760  lnct = lnct + 3
   IF (isecn /= 2) GO TO 770
   GO TO 820
   770  IF (iprtc == 1) WRITE (log2,780)
   780  FORMAT (/10X,8HPOINT no,4X,5HZSEMI,9X,5HXSEMI,9X,5HYSEMI, /)
   790  FORMAT (/10X,8HPOINT no,4X,5HZSEMI,9X,5HXSEMI,9X,5HYSEMI,13X,  &
       5HZSEMJ,9X,5HXSEMJ,9X,5HYSEMJ, /)
   i = 1
   800  IF (iprtc == 1) WRITE (log2,810) i,zsemi(j,i),xsemi(j,i), ysemi(j,i)
   810  FORMAT (10X,i5,3X,1P,3E14.5)
   GO TO 850
   820  IF (iprtc == 1) WRITE (log2,790)
   i = 1
   830  IF (iprtc == 1) WRITE (log2,840) i,zsemi(j,i),xsemi(j,i),  &
       ysemi(j,i),zsemj(j,i),xsemj(j,i),ysemj(j,i)
   840  FORMAT (10X,i5,3X,1P,3E14.5,4X,1P,3E14.5)
   850  i = i + 1
   lnct = lnct + 1
   IF (i > 31) CYCLE
   IF (lnct <= 59 .AND. isecn == 2) GO TO 830
   IF (isecn /= 2) GO TO 860
   IF (iprtc /= 0) WRITE (log2,200)
   IF (iprtc == 1) WRITE (log2,790)
   lnct = 4
   GO TO 830
   860  IF (lnct <= 59) GO TO 800
   IF (iprtc /= 0) WRITE (log2,200)
   IF (iprtc == 1) WRITE (log2,780)
   lnct = 4
   GO TO 800
 END DO
 IF (iprint == 1) GO TO 1030
 vol = 0.0
 DO  j = 2,nlines
   vol = vol + (((xs(j,1)-xp(j,1))**2+(ys(j,1)-yp(j,1))**2) +  &
       ((xs(j-1,1)-xp(j-1,1))**2 + (ys(j-1,1)-yp(j-1,1))**2))*  &
       (zs(j,1)+zp(j,1)-zs(j-1,1) - zp(j-1,1))*pi/32.0
   DO  i = 2,npoint
     vol = vol + ((SQRT((xs(j,i)-xp(j,i))**2+(ys(j,i)-yp(j,i))**2) +  &
         SQRT((xs(j,i-1)-xp(j,i-1))**2+(ys(j,i-1)-yp(j,i-1))**2))*  &
         (SQRT((xs(j,i-1)-xs(j,i))**2+(ys(j,i-1)-ys(j,i))**2) +  &
         SQRT((xp(j,i-1)-xp(j,i))**2+(yp(j,i-1)-yp(j,i))**2)) +  &
         (SQRT((xs(j-1,i)-xp(j-1,i))**2+(ys(j-1,i)-yp(j-1,i))**2) +  &
         SQRT((xs(j-1,i-1)-xp(j-1,i-1))**2+(ys(j-1,i-1)-yp(j-1,i-1))&
         **2))*(SQRT((xs(j-1,i-1)-xs(j-1,i))**2+(ys(j-1,i-1)-&
         ys(j-1,i))**2)+SQRT((xp(j-1,i-1)-xp(j-1,i))**2+(yp(j-1,i-1)-&
         yp(j-1,i))**2)))*(zs(j,i)+zs(j,i-1)+zp(j,i)+zp(j,i-1)-&
         zs(j-1,i)-zs(j-1,i-1)-zp(j-1,i)-zp(j-1,i-1))/32.0
   END DO
 END DO
 IF (lnct <= 56) GO TO 890
 lnct = 1
 IF (iprtc /= 0) WRITE (log2,200)
 890  lnct = lnct + 4
 IF (iprtc == 1) WRITE (log2,900) vol
 900  FORMAT (//40X,25HVOLUME of blade section =,1P,e11.4, /40X,36(1H*))
 IF (ijdum  == 0) GO TO 1030
 IF (iprint /= 3) WRITE (log2,200)
 IF (iprint == 3) WRITE (log2,910)
 910  FORMAT (//)
 IF (iprtc == 1) WRITE (log2,920)
 920  FORMAT (1H1,42X,43HBLADE calculations for aerodynamic analysis,  &
     /43X,43(1H*))
 idum = 7
 lnct = lnct + 4
 IF (iprint /= 3) lnct = 3
 DO  i = 1,nstns
   IF (ifangs(i) == 0 .OR. (isplit >= 1 .AND. ifangs(i) == 1)) CYCLE
   DO  j = 1,nlines
     CALL alg15 (rsta(1,i),xsta(1,i),kpts(i),r(i,j),xdum,1,0)
     CALL alg14 (rsta(1,i),xsta(1,i),kpts(i),r(i,j),xdum,zq(j),1,1)
     DO  k = 1,npoint
       ss(k)  = xs(j,k)
       rad(k) = ys(j,k)
       xtemp(k) = xp(j,k)
       x(k) = yp(j,k)
     END DO
     xdum = xdum - stackx
     CALL alg15 (ss,rad,npoint,xdum,yy1,1,1)
     CALL alg15 (xtemp,x,npoint,xdum,yy2,1,1)
     w1 = yy1/r(i,j)
     w2 = yy2/r(i,j)
     tq(j) = ABS(ATAN(w1/SQRT(1.-w1**2))-ATAN(w2/SQRT(1.-w2**2)))/  &
         (2.*pi)*FLOAT(nblade)
   END DO
   CALL alg14 (zcamb(1,i),ycamb(1,i),nlines,zcamb(1,i),xdum,rle, nlines,1)
   IF (lnct+idum+nlines <= 59) GO TO 950
   IF (iprtc /= 0) WRITE (log2,200)
   lnct = 2
   950  lnct = lnct + idum + nlines
   IF (iprtc == 1) WRITE (log2,960) i,nlines
   960  FORMAT (///48X,8HSTATION ,i2,5X,17HNUMBER of radii= ,i2,  //36X,&
       6Hradius,5X,7HSECTION,6X,4HLEAN,9X,5HBLADE,7X,5HTHETA, /48X,5HANGLE,  &
       6X,5HANGLE,7X,8HBLOCKAGE, /)
   DO  j = 1,nlines
     eps = (theta(j,i)-ATAN(rle(j)))*c1
     alphb = alpha(j,i)
     alp = (ATAN((tanphi(i,j)*TAN(eps/c1)+alphb*SQRT(1.+tanphi(i,j)**2))  &
         /(1.-tanphi(i,j)*zq(j))))*c1
     alpb(i,j) = alp
     epslon(i,j) = ATAN(TAN(eps/c1)/SQRT(1.0+zq(j)**2))*c1
     IF (isplit < 1) GO TO 990
     CALL fread (log1,rdata,4,1)
     xb = rdata(4)
     IF (iprtc == 1) WRITE (log2,980) xb,i,j
     980  FORMAT(90X,14HADDIT. BLOCK =,f7.5,3H i=,i2,3H j=,i2)
     tq(j) = tq(j) + xb
     990  IF (iprtc == 1) WRITE (log2,1010) r(i,j),alp,eps,tq(j),theta(j,i)
     BLOCK(i,j) = tq(j)
   END DO
   1010 FORMAT (30X,5F12.4)
 END DO
 1030 IF (ifplot < 2 .OR. ifplot == 4) GO TO 1040
 CALL alg17 (istak,pltsze,2,title,ikdum,ifplot)
 1040 IF (iprint == 1 .OR. iprint == 3) GO TO 1060
 lnct = 2
 IF (iprtc == 1) WRITE (log2,1050)
 1050 FORMAT (1H1,27X,74HBLADE surface geometry in cartesian coordinates &
     &at specified values of  z , /28X,18(4H****),2H**)
 1060 IF (iprint == 1 .AND. ifplot <= 1) GO TO 1470
 xz = nz - 1
 dz = (zouter-zinner)/xz
 zout(1) = zinner
 DO  j = 3,nz
   zout(j-1) = zout(j-2) + dz
 END DO
 zout(nz) = zouter
 DO  i = 1,npoint
   CALL alg15 (zs(1,i),xs(1,i),nlines,zout,temp1,nz,0)
   CALL alg15 (zs(1,i),ys(1,i),nlines,zout,temp2,nz,0)
   CALL alg15 (zp(1,i),xp(1,i),nlines,zout,temp3,nz,0)
   CALL alg15 (zp(1,i),yp(1,i),nlines,zout,temp4,nz,0)
   DO  j = 1,nz
     xs(j,i) = temp1(j)
     ys(j,i) = temp2(j)
     xp(j,i) = temp3(j)
     yp(j,i) = temp4(j)
   END DO
 END DO
 DO  i = 1,31
   CALL alg15 (zsemi(1,i),xsemi(1,i),nlines,zout,temp1,nz,0)
   CALL alg15 (zsemi(1,i),ysemi(1,i),nlines,zout,temp2,nz,0)
   DO  j = 1,nz
     xsemi(j,i) = temp1(j)
     ysemi(j,i) = temp2(j)
   END DO
 END DO
 IF (isecn /= 2) GO TO 1110
 DO  i = 1,31
   CALL alg15 (zsemj(1,i),xsemj(1,i),nlines,zout,temp1,nz,0)
   CALL alg15 (zsemj(1,i),ysemj(1,i),nlines,zout,temp2,nz,0)
   DO  j = 1,nz
     xsemj(j,i) = temp1(j)
     ysemj(j,i) = temp2(j)
   END DO
 END DO
 1110 DO  j = 1,nz
   rd = SQRT((xs(j,1)-xp(j,1))**2+(ys(j,1)-yp(j,1))**2)/2.0
   area = pi*rd**2/2.0
   beta1 = ATAN((ys(j,2)+yp(j,2)-ys(j,1)-yp(j,1))/(xs(j,2)+xp(j,2)-&
       xs(j,1)-xp(j,1)))
   xint = area*((xp(j,1)+xs(j,1))/2.0-COS(beta1)*4.0/(3.0*pi)*rd)
   yint = area*((yp(j,1)+ys(j,1))/2.0-SIN(beta1)*4.0/(3.0*pi)*rd)
   IF (isecn /= 2) GO TO 1120
   n1 = npoint
   n  = n1
   n2 = n1 - 1
   beta2 = ATAN((ys(j,n1)+yp(j,n1)-ys(j,n2)-yp(j,n2))/(xs(j,n1)+&
       xp(j,n1)-xs(j,n2)-xp(j,n2)))
   xint = xint + area*((xp(j,n)+xs(j,n))/2.+COS(beta2)*4./(3.*pi)*rd)
   yint = yint + area*((yp(j,n)+ys(j,n))/2.+SIN(beta2)*4./(3.*pi)*rd)
   area = 2.*area
   1120 DO  i = 2,npoint
     dela = (SQRT((xs(j,i)-xp(j,i))**2+(ys(j,i)-yp(j,i))**2)+&
         SQRT((xs(j,i-1)-xp(j,i-1))**2+(ys(j,i-1)-yp(j,i-1))**2))*&
         (SQRT((xs(j,i-1)-xs(j,i))**2+(ys(j,i-1)-ys(j,i))**2)+&
         SQRT((xp(j,i-1)-xp(j,i))**2+(yp(j,i-1)-yp(j,i))**2))/4.0
     area = area + dela
     xint = xint + dela*(xs(j,i)+xs(j,i-1)+xp(j,i)+xp(j,i-1))/4.0
     yint = yint + dela*(ys(j,i)+ys(j,i-1)+yp(j,i)+yp(j,i-1))/4.0
   END DO
   yint = yint/area
   xint = xint/area
   x1   = (xs(j,1)+xp(j,1))/2.
   y1   = (ys(j,1)+yp(j,1))/2.
   t1   = SQRT((xs(j,1)-xp(j,1))**2+(ys(j,1)-yp(j,1))**2)
   f    = 0.
   u    = 0.
   DO  i = 2,npoint
     t2   = SQRT((xs(j,i)-xp(j,i))**2+(ys(j,i)-yp(j,i))**2)
     x2   = (xs(j,i)+xp(j,i))/2.
     y2   = (ys(j,i)+yp(j,i))/2.
     delu = SQRT((x2-x1)**2+(y2-y1)**2)
     u    = u + delu
     tav3 = (t1**3+t2**3)/2.
     f    = f + tav3*delu
     x1   = x2
     y1   = y2
     t1   = t2
   END DO
   torcon = ((1./3.)*f)/(1.+(4./3.)*f/area/u**2)
   ix   = 0.0
   iy   = 0.0
   ixy  = 0.0
   DO  i = 2,npoint
     xd   = (SQRT((xs(j,i-1)-xp(j,i-1))**2+(ys(j,i-1)-yp(j,i-1))**2)+&
         SQRT((xs(j,i)-xp(j,i))**2+(ys(j,i)-yp(j,i))**2))/2.0
     yd   = (SQRT((xs(j,i)-xs(j,i-1))**2+(ys(j,i)-ys(j,i-1))**2)+&
         SQRT((xp(j,i)-xp(j,i-1))**2+(yp(j,i)-yp(j,i-1))**2))/2.0
     ixd  = yd*yd*yd*xd/12.0
     iyd  = xd*xd*xd*yd/12.0
     ang  = ATAN((ys(j,i)+yp(j,i)-ys(j,i-1)-yp(j,i-1))/(xp(j,i)+xs(j,i)&
         -xp(j,i-1)-xs(j,i-1)))
     cosang = COS(2.0*ang)
     ixn  = (ixd+iyd+(ixd-iyd)*cosang)/2.0
     iyn  = (ixd+iyd-(ixd-iyd)*cosang)/2.0
     ixyn = 0.0
     IF (ang /= 0.0) ixyn = ((ixn-iyn)*cosang-ixd+iyd)/ (2.0*SIN(2.0*ang))
     dela = xd*yd
     ymn  = (ys(j,i)+ys(j,i-1)+yp(j,i)+yp(j,i-1))/4.0-yint
     xmn  = (xs(j,i)+xs(j,i-1)+xp(j,i)+xp(j,i-1))/4.0-xint
     ix   = ix + ixn + dela*ymn*ymn
     iy   = iy + iyn + dela*xmn*xmn
     ixy  = ixy+ ixyn+ dela*ymn*xmn
   END DO
   ang  = ATAN(2.0*ixy/(iy-ix))
   ipx  = (ix+iy)/2.0+(ix-iy)/2.0*COS(ang)-ixy*SIN(ang)
   ipy  = (ix+iy)/2.0-(ix-iy)/2.0*COS(ang)+ixy*SIN(ang)
   ang  = ang/2.0*c1
   IF (iprint == 1 .OR. iprint == 3) GO TO 1320
   IF (lnct <= 45) GO TO 1160
   IF (iprtc /= 0) WRITE (log2,200)
   lnct = 1
   1160 lnct = lnct + 16
   IF (iprtc == 1) WRITE (log2,1170) j,zout(j),area,xint,yint,ix,  &
       iy,ixy,ipx,ang,ipy,ang
   1170 FORMAT (/50X,14HSECTION NUMBER,i3,3X,5H z  =,f9.4, /50X,    &
       34H**********************************, ///20X,18HSECTION properties,7X,&
       12HSECTION area,26X,1H=,1P,e12.4,//45X,20HLOCATION of centroid,11X,  &
       4HXBAR,3X,1H=,e12.4, /45X,22HRELATIVE to stack axis,9X,4HYBAR,3X,&
       1H=,e12.4, //45X,22HSECOND moments of area,9X,2HIX,5X,1H=,e12.4,&
       /45X,14HABOUT centroid,17X,2HIY,5X,1H=,e12.4, /76X,3HIXY,4X,1H=,e12.4 &
       //45X,24HPRINCIPAL second moments,7X,3HIPX,4X,1H=,e12.4,4H (at,&
       70P,f7.2,21H degrees TO  x  axis),/45X,22HOF area about centroid,9X  &
       ,3HIPY,4X,1H=,1P,e12.4,4H (at,0P,f7.2,21H degrees to  y  axis))
   IF (iprtc == 1) WRITE (log2,1180) torcon
   1180 FORMAT (/45X,18HTORSIONAL constant,20X,1H=,1P,e12.4, /)
   lnct = lnct + 3
   IF (lnct <= 50) GO TO 1190
   IF(iprtc /=  0) WRITE (log2,200)
   lnct = 1
   1190 lnct = lnct + 5
   IF (iprtc == 1) WRITE (log2,1200)
   1200 FORMAT (/20X,19HSECTION coordinates, /)
   IF (iprtc == 1) WRITE (log2,1210)
   1210 FORMAT (31X,8HPOINT no,5X,2HX1,12X,2HY1,16X,2HX2,12X,2HY2, /)
   DO  i = 1,npoint
     lnct = lnct + 1
     IF (lnct <= 60) GO TO 1220
     lnct = 4
     IF (iprtc /= 0) WRITE (log2,200)
     IF (iprtc == 1) WRITE (log2,1210)
     1220 IF (iprtc == 1) WRITE (log2,1230) i,xs(j,i),ys(j,i),xp(j,i),  &
         yp(j,i)
   END DO
   1230 FORMAT (31X,i5,3X,1P,2E14.5,4X,1P,2E14.5)
   IF (lnct <= 55) GO TO 1240
   lnct = 1
   IF (iprtc /= 0) WRITE (log2,200)
   1240 lnct = lnct + 3
   IF (iprtc == 1 .AND. isecn == 2) WRITE (log2,1260)
   IF (isecn == 2) GO TO 1270
   IF (iprtc == 1) WRITE (log2,1250)
   1250 FORMAT (/31X,8HPOINT no,5X,5HXSEMI,9X,5HYSEMI, /)
   1260 FORMAT (/31X,8HPOINT no,5X,5HXSEMI,9X,5HYSEMI,12X,5HXSEMJ,9X,  &
       5HYSEMJ, /)
   1270 DO  i = 1,31
     lnct = lnct + 1
     IF (lnct <= 60) GO TO 1290
     IF (iprtc /= 0) WRITE (log2,200)
     IF (iprtc == 1 .AND. isecn == 2) WRITE (log2,1260)
     IF (isecn == 2) GO TO 1280
     IF (iprtc == 1) WRITE (log2,1250)
     1280 lnct = 4
     1290 IF (iprtc == 1 .AND. isecn == 2) WRITE (log2,1230) i,xsemi(j,i),  &
         ysemi(j,i),xsemj(j,i),ysemj(j,i)
     IF (isecn == 2) CYCLE
     IF (iprtc == 1) WRITE (log2,1310) i,xsemi(j,i),ysemi(j,i)
   END DO
   1310 FORMAT (31X,i5,3X,1P,2E14.5)
   1320 IF (ifplot < 2) CYCLE
   IF (ifplot == 4) GO TO 1380
   xplot = xs(j,1)*scale
   yplot = ys(j,1)*scale
   DO  i = 2,npoint
     xplot = xs(j,i)*scale
     yplot = ys(j,i)*scale
   END DO
   IF (isecn /= 2) GO TO 1350
   DO  i = 2,30
     xplot = xsemj(j,i)*scale
     yplot = ysemj(j,i)*scale
   END DO
   1350 DO  ii = 1,npoint
     i = npoint + 1 - ii
     xplot = xp(j,i)*scale
     yplot = yp(j,i)*scale
   END DO
   DO  i = 2,30
     xplot = xsemi(j,i)*scale
     yplot = ysemi(j,i)*scale
   END DO
   xplot = xs(j,1)*scale
   yplot = ys(j,1)*scale
   CYCLE
   1380 CONTINUE
   xj = j
   stager = ATAN((ys(j,npoint)+yp(j,npoint)-ys(j,1)-yp(j,1))/  &
       (xs(j,npoint)+xp(j,npoint)-xs(j,1)-xp(j,1)))*c1
   xsign  = FLOAT(nsign)
   sinstg = SIN(stager/c1)
   cosstg = COS(stager/c1)
   yplot  = 4.75
   xplot  = 4.75*sinstg/cosstg
   IF (ABS(xplot) <= 22.0) GO TO 1390
   xplot  = 22.0
   yplot  =-22.0/sinstg*cosstg
   1390 CONTINUE
   xplot = -xplot
   yplot = -yplot
   xplot = 22.0
   yplot =-22.0*sinstg/cosstg
   IF (ABS(yplot) <= 4.75) GO TO 1400
   yplot =-4.75
   xplot = 4.75/sinstg*cosstg
   1400 CONTINUE
   xplot = -xplot
   yplot = -yplot
   xplot = scale*(xs(j,1)*cosstg+ys(j,1)*sinstg)
   yplot = scale*(ys(j,1)*cosstg-xs(j,1)*sinstg)
   DO  i = 2,npoint
     xplot = scale*(xs(j,i)*cosstg+ys(j,i)*sinstg)
     yplot = scale*(ys(j,i)*cosstg-xs(j,i)*sinstg)
   END DO
   IF (isecn /= 2) GO TO 1430
   DO  i = 2,30
     xplot = scale*(xsemj(j,i)*cosstg+ysemj(j,i)*sinstg)
     yplot = scale*(ysemj(j,i)*cosstg-xsemj(j,i)*sinstg)
   END DO
   1430 DO  ii = 1,npoint
     i = npoint + 1 - ii
     xplot = scale*(xp(j,i)*cosstg+yp(j,i)*sinstg)
     yplot = scale*(yp(j,i)*cosstg-xp(j,i)*sinstg)
   END DO
   DO  i = 2,30
     xplot = scale*(xsemi(j,i)*cosstg+ysemi(j,i)*sinstg)
     yplot = scale*(ysemi(j,i)*cosstg-xsemi(j,i)*sinstg)
   END DO
   xplot = scale*(xs(j,1)*cosstg+ys(j,1)*sinstg)
   yplot = scale*(ys(j,1)*cosstg-xs(j,1)*sinstg)
 END DO
 1470 CONTINUE
 IF (inast == 0) GO TO 1580
 xsign = FLOAT(nsign)
 WRITE  (log2,1471)
 1471 FORMAT (1H0,10X,34HNASTRAN compressor blade bulk DATA , /10X,  &
     36(1H*), //)
 IF (ipgeom == 1) GO TO 1562
 WRITE  (log2,1472)
 1472 FORMAT (11X,30H*** ctria2 AND ptria2 DATA ***, /)
 nstad = irte - irle + 1
 jloop = 0
 nelem = 0
 nstrd = nlines - 1
 irt   = irte - 1
 nt    = 1995
 DO  j = 1,nstrd
   DO  i = irle,irt
     nelem = nelem + 1
     igrd1 = i - 1 + jloop
     igrd3 = igrd1 + nstad
     igrd2 = igrd1 + nstad + 1
     nt    = nt + 5
     WRITE (lpunch,1530) nelem,nt,igrd1,igrd2,igrd3
     WRITE (log2,1531) nelem,nt,igrd1,igrd2,igrd3
     IF (ABS(FLOAT(inast)) > 3.5) GO TO 1480
     thck = (tharr(j,i)+tharr(j+1,i)+tharr(j+1,i+1))/3.
     pres =-xsign*(blafor(i,j)+blafor(i,j+1)+blafor(i+1,j+1))/3.
     GO TO 1490
     1480 thck = (tharr(j,i)+tharr(j+1,i)+tharr(j+1,i+1)+tharr(j,i+1))/4.
     pres =-xsign*(blafor(i,j)+blafor(i,j+1)+blafor(i+1,j+1)+  &
         blafor(i+1,j))/4.
     1490 WRITE (lpunch,1540) nt,thck
     WRITE (log2,1541)   nt,thck
     IF (inast > 0) WRITE (lpunch,1550) pres,nelem
     IF (inast > 0) WRITE (log2,1551)   pres,nelem
     nelem = nelem + 1
     igrd3 = igrd2
     igrd2 = igrd1 + 1
     IF (ABS(FLOAT(inast)) > 3.5) GO TO 1500
     nt   = nt + 5
     thck = (tharr(j,i)+tharr(j,i+1)+tharr(j+1,i+1))/3.
     pres = -xsign*(blafor(i,j)+blafor(i+1,j)+blafor(i+1,j+1))/3.
     WRITE (lpunch,1540) nt,thck
     WRITE (log2,1541)   nt,thck
     1500 WRITE (lpunch,1530) nelem,nt,igrd1,igrd2,igrd3
     WRITE (log2,1531)   nelem,nt,igrd1,igrd2,igrd3
     IF (inast > 0) WRITE (lpunch,1550) pres,nelem
     IF (inast > 0) WRITE (log2,1551)   pres,nelem
   END DO
   jloop = jloop + nstad
 END DO
 1530 FORMAT (6HCTRIA2,7X,i3,4X,i4,3(5X,i3))
 1531 FORMAT (1X,6HCTRIA2,7X,i3,4X,i4,3(5X,i3))
 1540 FORMAT (6HPTRIA2,6X,i4,7X,1H1,f8.4,6X,2H0.)
 1541 FORMAT (1X,6HPTRIA2,6X,i4,7X,1H1,f8.4,6X,2H0.)
 1550 FORMAT (6HPLOAD2,8X,2H60,f8.4,5X,i3)
 1551 FORMAT (1X,6HPLOAD2,8X,2H60,f8.4,5X,i3)
 1560 FORMAT (4HGRID,9X,i3,8X,3F8.4)
 1561 FORMAT (1X,4HGRID,9X,i3,8X,3F8.4)
 1562 CONTINUE
 WRITE  (log2,1563)
 1563 FORMAT (1H0,10X,29H*** blade grid point DATA *** ,/)
 jd = 0
 DO  j = 1,nlines
   DO  i = irle,irte
     jd = jd + 1
     ycamb(j,i) = -xsign*ycamb(j,i)
     WRITE (log2,1561) jd,xcamb(j,i),ycamb(j,i),zcamb(j,i)
     WRITE (lpunch,1560) jd,xcamb(j,i),ycamb(j,i),zcamb(j,i)
   END DO
 END DO
 IF (istrml == -1 .OR. istrml == 2) GO TO 1580
 WRITE  (log2,1571)
 1571 FORMAT (1H0,10X,27H*** blade streaml1 DATA ***,/)
 nstad  = irte - irle + 1
 nstad1 = nstad - 1
 DO  j = 1,nlines
   nd1 = (j-1)*nstad + 1
   nd2 = nd1 + nstad1
   WRITE  (lpunch, 1573) j,nd1,nd2
   WRITE  (log2,1574) j,nd1,nd2
 END DO
 1573 FORMAT (8HSTREAML1,i8,i8,8H thru   ,i8)
 1574 FORMAT (1X,8HSTREAML1,i8,i8,8H thru   ,i8)
 1580 CONTINUE
 IF (naero == 1 .OR. ipunch == 1) CALL alg19 (log1,log2,log3,log5,  &
     nlines,nspec,kpts,rsta,xsta,r,zr,b1,b2,tc,pi,c1,nblade,ccord,  &
     block,alpb,epslon,ifangs,ipunch,naero)
!     IF (IFPLOT .NE. 0) CALL PLOT (0.0,0.0,-3)

 RETURN
END SUBROUTINE algan
