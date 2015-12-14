SUBROUTINE mbreg (ireg,nw1,nwn,nc21,nc2n,nc1,ncn,nd1,ndn,xk,yk,  &
        xk1,yk1,xk2,yk2,xwte,ywte,kte,kte1,kte2,parea)
     
!     SUBROUTINE TO COMPUTE LIMITS OF REGION AND PERCENTAGE OF BOX IN
!     EACH
 
 
 INTEGER, INTENT(OUT)                     :: ireg
 INTEGER, INTENT(OUT)                     :: nw1(1)
 INTEGER, INTENT(OUT)                     :: nwn(1)
 INTEGER, INTENT(OUT)                     :: nc21(1)
 INTEGER, INTENT(OUT)                     :: nc2n(1)
 INTEGER, INTENT(OUT)                     :: nc1(1)
 INTEGER, INTENT(OUT)                     :: ncn(1)
 INTEGER, INTENT(OUT)                     :: nd1(1)
 INTEGER, INTENT(OUT)                     :: ndn(1)
 REAL, INTENT(OUT)                        :: xk(1)
 REAL, INTENT(OUT)                        :: yk(1)
 REAL, INTENT(OUT)                        :: xk1(1)
 REAL, INTENT(OUT)                        :: yk1(1)
 REAL, INTENT(OUT)                        :: xk2(1)
 REAL, INTENT(OUT)                        :: yk2(1)
 REAL, INTENT(OUT)                        :: xwte(1)
 REAL, INTENT(OUT)                        :: ywte(1)
 INTEGER, INTENT(OUT)                     :: kte(1)
 INTEGER, INTENT(OUT)                     :: kte1(1)
 INTEGER, INTENT(OUT)                     :: kte2(1)
 REAL, INTENT(OUT)                        :: parea(50,50,3)
 LOGICAL :: cntrl2,cntrl1,crank1,crank2,asym,debug
 
 
 COMMON /system/ sys,n6
 COMMON /mboxa / x(12),y(12),tang(10),ang(10),cotang(10)
 COMMON /mboxc / njj ,crank1,crank2,cntrl1,cntrl2,nbox,npts0,npts1,  &
     npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,boxl,boxw,  &
     boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 DATA    debug / .false./
 
 iprint = 0
 IF (debug) iprint = 1
 
 ireg = 1
 boxa = boxl*boxw
 kpt  = 0
 kc1t = 0
 kc2t = 0
 DO  i = 1,50
   nw1(i)  = 0
   nwn(i)  = 0
   nc1(i)  = 0
   ncn(i)  = 0
   nc21(i) = 0
   nc2n(i) = 0
   nd1(i)  = 0
   ndn(i)  = 0
   kte(i)  = 0
   kte1(i) = 0
   kte2(i) = 0
   xwte(i) = 0.0
   ywte(i) = 0.0
   DO  j = 1,50
     DO  kp = 1,3
       parea(i,j,kp) = 0.
     END DO
   END DO
 END DO
 DO  i = 1,200
   xk(i) = 0.0
   yk(i) = 0.0
 END DO
 DO  i = 1,125
   xk1(i) = 0.0
   yk1(i) = 0.0
   xk2(i) = 0.0
   yk2(i) = 0.0
 END DO
 
!     LEADING EDGE OF MAIN
 
 xre = 0.0
 ybe = 0.0
 k   = 1
 yr  = -0.5*boxw
 DO  i = 1,nsb
   yl  = yr
   yr  = (FLOAT(i)-0.5)*boxw
   xle = xre
   IF (yr > y(k+1)) GO TO 50
   xre = (yr-y(k))*tang(k) + x(k)
   GO TO 60
   50 xre = (yr-y(k+1))*tang(k+1) + x(k+1)
   kpt = 1
   60 xt  = xle - AMOD(xle,boxl)
   xb  = xt + boxl
   j1  = xb/boxl + 0.01
   
   
   DO  j = j1,ncb
     IF (xre > xb) GO TO 100
     IF (xle > xt) GO TO 90
     
     IF (xt <= x(k+1)) GO TO 70
     yte  = (xt-x(k+1))*cotang(k+1) + y(k+1)
     GO TO 80
     70 yte  = (xt-x(k))*cotang(k) + y(k)
     IF (kpt == 1) kpt =  2
     80 a = 0.5*(yr-yte)*(xre-xt)
     IF (kpt == 2) a = a + (xt*(y(k+1)-yr) - yte*(x(k+1)-xre) +  &
         x(k+1)*yr - xre*y(k+1))/2.0
     pa = 1.0 - a/boxa
     GO TO 180
     
     90 a  = 0.5*(xle+xre-2.0*xt)*(yr-yl)
     IF (kpt > 0) a = a + (xle*(y(k+1)-yr) - yl*(x(k+1)-xre) +  &
         x(k+1)*yr - xre*y(k+1))/2.0
     pa = 1.0 - a/boxa
     GO TO 180
     
     100 IF (xle > xt) GO TO 130
     
     yte = ybe
     IF (xb > x(k+1)) GO TO 110
     ybe = (xb-x(k))*cotang(k) + y(k)
     GO TO 120
     110 ybe = (xb-x(k+1))*cotang(k+1) + y(k+1)
     IF (kpt == 1) kpt  =  2
     120 a = 0.5*(yte+ybe-2.0*yl)*boxl
     IF (kpt == 2) a = a + (xt*(ybe-y(k+1)) - yte*(xb-x(k+1)) +  &
         xb*y(k+1) - ybe*x(k+1))/2.0
     pa = a/boxa
     GO TO 160
     
     130 IF (xb > x(k+1)) GO TO 140
     ybe = (xb-x(k))*cotang(k) + y(k)
     GO TO 150
     140 ybe = (xb-x(k+1))*cotang(k+1) + y(k+1)
     IF (kpt == 1) kpt =  2
     150 a   = 0.5*(xb-xle)*(ybe-yl)
     IF (kpt == 2) a = a + (xle*(ybe-y(k+1)) - yl*(xb-x(k+1)) +  &
         xb*y(k+1) - ybe*x(k+1))/2.0
     pa  = a/boxa
     160 xt  = xb
     xb  = FLOAT(j+1)*boxl
     IF (kpt == 2) kpt = 3
     IF (i == 1) pa = 2.0*pa - 1.0
     parea(j,i,1) = pa
   END DO
   GO TO 190
   
   180 IF (i == 1) pa = 2.0*pa - 1.0
   parea(j,i,1) = pa
   190 yc  = yr - 0.5*boxw
   IF (kpt <= 0) GO TO 200
   IF (yc <= y(k+1)) GO TO 200
   xc  = (yc-y(k+1))*tang(k+1) + x(k+1)
   GO TO 210
   200 xc  =  (yc-y(k))*tang(k) + x(k)
   210 nw1(i) = xle/boxl + 1.0001
   IF (kpt > 0) k = k + 1
   kpt = 0
 END DO
 IF (iprint <= 0) GO TO 250
 WRITE  (n6,240)
 WRITE  (n6,230) (nw1(i),i=1,nsb)
 230 FORMAT (10I12)
 240 FORMAT (4H nw1)
 250 CONTINUE
 
!     TRAILING EDGE OF MAIN
 
 xre = x(4)
 k   = 4
 yr  = 0.0
 DO  i = 1,nsb
   yl  = yr
   yr  = (FLOAT(i)-0.5)*boxw
   xle = xre
   IF (yr > y(k+1)) GO TO 260
   xre = (yr-y(k))*tang(k) + x(k)
   GO TO 270
   260 xre = (yr-y(k+1))*tang(k+1) + x(k+1)
   kpt = 1
   270 xt  = xle - AMOD(xle,boxl)
   xb  = xt + boxl
   j   = xb/boxl + 0.01
   IF (j > 50) GO TO 410
   ipt = 0
   IF (xre > xb .OR. xre < xt) GO TO 280
   a   = 0.5*(xle+xre-2.0*xt)*(yr-yl)
   IF (kpt > 0) a = a + (xle*(y(k+1)-yr) - yl*(x(k+1)-xre) +  &
       x(k+1)*yr - y(k+1)*xre)/2.0
   IF (i == 1) a = 2.0*a
   GO TO 430
   
   280 ipt = 1
   IF (xle < xre) GO TO 340
   ipt = -1
   290 IF (xre < xt) GO TO 300
   
   a  =  0.5*(xb-xre)*(yr-yte)
   IF (kpt > 0 .AND. kpt < 3) a = a - (xre*(yte - y(k+1)) -  &
       yr*(xb-x(k+1)) + xb*y(k+1) - yte*x(k+1))/2.0
   IF (i == 1) a = 2.0*a
   GO TO 420
   
   300 ybe = yte
   IF (xt < x(k+1)) GO TO 310
   yte = (xt-x(k))*cotang(k) + y(k)
   GO TO 320
   310 yte = (xt-x(k+1))*cotang(k+1) + y(k+1)
   IF (kpt == 1) kpt =  2
   320 IF (xle > xb) GO TO 330
   
   a = 0.5*(xle-xt)*(yte-yl)
   IF (kpt == 2) a = a + (xt*(yl-y(k+1)) - yte*(xle-x(k+1)) +  &
       xle*y(k+1) - yl*x(k+1))/2.0
   GO TO 390
   
   330 a = 0.5*boxl*(yte+ybe - 2.0*yl)
   IF (kpt == 2) a = a + (xt*(ybe-y(k+1)) - yte*(xb- x(k+1)) +  &
       xb*y(k+1) - ybe*x(k+1))/2.0
   GO TO 390
   
   340 IF (xre > xb) GO TO 350
   
   a = 0.5*(yr-ybe)*(xre-xt)
   IF (kpt > 0 .AND. kpt < 3) a = a + (xt*(y(k+1)-yr) -  &
       ybe*(x(k+1)-xre ) + x(k+1)*yr - y(k+1)*xre)/2.0
   IF (i == 1) a = 2.0*a
   GO TO 430
   
   350 yte = ybe
   IF (xb > x(k+1)) GO TO 360
   ybe = (xb-x(k))*cotang(k) + y(k)
   GO TO 370
   360 ybe = (xb-x(k+1))*cotang(k+1) + y(k+1)
   IF (kpt ==  1) kpt = 2
   370 IF (xle < xt) GO TO 380
   
   a  = 0.5*(xb-xle)*(ybe-yl)
   IF (kpt == 2) a = a - (xle*(y(k+1)-ybe) - yl* (x(k+1)-xb) +  &
       x(k+1)*ybe - y(k+1)*xb)/2.0
   IF (i == 1) a = 2.0*a
   a  =  boxa - a
   GO TO 400
   
   380 a  =  0.5*boxl*(2.0*yr - yte - ybe)
   IF (kpt == 2) a = a + (xt*(y(k+1)-ybe) - yte*(x(k+1)-xb) +  &
       x(k+1)*ybe - y(k+1)*xb)/2.0
   
   390 IF (i == 1) a = 2.0*a
   400 pa = a/boxa
   a  = 1.0
   IF (parea(j,i,1) > 0.0) a = parea(j,i,1)
   parea(j,i,1) = pa*a
   j  = j + ipt
   IF (j > 50) GO TO 410
   xb = FLOAT(j)*boxl
   xt = xb - boxl
   IF (kpt == 2) kpt = 3
   IF (ipt > 0) GO TO 340
   GO TO 290
   410 ireg = 2
   RETURN
   
   420 a  = boxa - a
   430 yc = yr - 0.5*boxw
   IF (kpt <= 0) GO TO 440
   IF (yc <= y(k+1)) GO TO 440
   xc = (yc-y(k+1))*tang(k+1) + x(k+1)
   GO TO 450
   440 xc = (yc-y(k))*tang(k) + x(k)
   450 nwn(i) = AMAX1(xle,xre)/boxl + 0.9999
   pa = a/boxa
   a  = 1.0
   IF (parea(j,i,1) > 0.0) a = parea(j,i,1)
   parea(j,i,1) = pa*a
   xwte(i) = xc
   ywte(i) = yc
   IF (kpt > 0) k = k + 1
   kpt =  0
 END DO
 IF (iprint <= 0) GO TO 480
 WRITE  (n6,470)
 WRITE  (n6,230) (nwn(i),i=1,nsb)
 470 FORMAT (4H nwn)
 480 CONTINUE
 ntote = nsb
 
!     FILL IN MAIN PERCENTAGES
 
 loop520:  DO  i = 1,nsb
   n1 = nw1(i)
   nn = nwn(i)
   DO  j = n1,nn
     IF (parea(j,i,1) <= 0.0) parea(j,i,1) = 1.0
   END DO
   
!     DIAPHRAGM INDEX
   
   IF (i /= 1) GO TO 500
   nd1(1) = nw1(1)
   ndn(1) = nwn(1)
   CYCLE loop520
   500 nd1(i) = MIN0(nw1(i),nd1(i-1)+1)
   ndn(i) = MAX0(nwn(i),ndn(i-1)-1)
   IF (ndn(i) <= ndn(i-1)+1) CYCLE loop520
   DO  k = 2,i
     kk = i - k + 1
     IF (ndn(kk) >= ndn(kk+1)-1) CYCLE loop520
     ndn(kk) = MAX0(ndn(kk),ndn(kk+1)-1)
   END DO
 END DO loop520
 j  =  nsb + 1
 530 IF (nd1(j-1) >= ndn(j-1)-1) GO TO 540
 nd1(j) = nd1(j-1) + 1
 ndn(j) = ndn(j-1) - 1
 j  = j + 1
 IF (j <= 50) GO TO 530
 ireg =  2
 RETURN
 
 540 nsbd =  j - 1
 IF (iprint <= 0) GO TO 580
 WRITE  (n6,550)
 WRITE  (n6,230) (nd1(i),i=1,nsbd)
 550 FORMAT (4H nd1)
 WRITE  (n6,560)
 WRITE  (n6,230) (ndn(i),i=1,nsbd)
 560 FORMAT (4H ndn)
 
 WRITE  (n6,610)
 DO  i = 1,ncb
   WRITE (n6,640) i
   WRITE (n6,600) (parea(i,j,1),j=1,nsb)
 END DO
 580 IF (cntrl1) CALL mbctr (1,il1,ir1,ncn,nc1,nwn,nw1,parea)
 IF (cntrl2) CALL mbctr (2,il2,ir2,nc2n,nc21,nwn,nw1,parea)
 IF (iprint == 0) GO TO 650
 DO  kxyz = 1,3
   IF (kxyz == 1) WRITE (n6,610)
   IF (kxyz == 2) WRITE (n6,620)
   IF (kxyz == 3 ) WRITE (n6,630)
   DO  i = 1,ncb
     WRITE (n6,640) i
     WRITE (n6,600) (parea(i,j,kxyz),j=1,nsb)
   END DO
 END DO
 600 FORMAT (5X,10F9.5)
 610 FORMAT (12H parea, main)
 620 FORMAT (14H parea, cntrl1)
 630 FORMAT (14H parea, cntrl2)
 640 FORMAT (4H row,i4)
 
!     MAIN BOX CTR. COORDINATES
 
 650 kc = 0
 DO  i = 1,ncb
   ixr = i - 1
   DO  j = 1,nsb
     IF (.NOT.(i >= (nd1(j)) .AND. i <= (ndn(j)))) CYCLE
     IF (parea(i,j,1) < 0.005) CYCLE
     jxr = j  - 1
     kc  = kc + 1
     IF (kc >= 200) GO TO 830
     xk(kc) = boxl*(FLOAT(ixr)+0.5)
     yk(kc) = boxw*FLOAT(jxr)
   END DO
 END DO
 DO  j = 1,nsb
   kc = kc + 1
   IF (kc >= 200) GO TO 830
   kte(j) = kc
   xk(kc) = xwte(j)
   yk(kc) = ywte(j)
 END DO
 kct = kc
 IF (iprint <= 0) GO TO 700
 WRITE  (n6,680) (i,xk(i),i=1,kc)
 680 FORMAT (1H1,23H main box ctr. x coord.,/(10(1X,i4,f8.2)))
 WRITE  (n6,690) (i,yk(i),i=1,kc)
 690 FORMAT (1H1,23H main box ctr. y coord.,/(10(1X,i4,f8.2)))
 700 CONTINUE
 
!     CNTRL1 BOX CTR. COORDINATES
 
 IF (.NOT.cntrl1) GO TO 750
 kc1 = 0
 DO  i = 1,ncb
   ixr = i - 1
   DO  j = il1,ir1
     IF (parea(i,j,2) < 0.005) CYCLE
     jxr = j - 1
     kc1 = kc1 + 1
     IF (kc1 >= 125) GO TO 830
     xk1(kc1) = boxl*(FLOAT(ixr)+0.5)
     yk1(kc1) = boxw*FLOAT(jxr)
   END DO
 END DO
 DO  j = il1,ir1
   kc1 = kc1 + 1
   IF (kc1 >= 125) GO TO 830
   kte1(j)  = kc1
   xk1(kc1) = xwte(j)
   yk1(kc1) = ywte(j)
 END DO
 kc1t = kc1
 IF (iprint <= 0) GO TO 750
 WRITE  (n6,730) (i,xk1(i),i=1,kc1)
 730 FORMAT (1H1,25H cntrl1 box ctr. x coord.,/(10(1X,i4,f8.2)))
 WRITE  (n6,740) (i,yk1(i),i=1,kc1)
 740 FORMAT (1H1,25H cntrl1 box ctr. y coord.,/(10(1X,i4,f8.2)))
 
!     CNTRL2 BOX CTR. COORDINATES
 
 750 IF (.NOT.cntrl2) GO TO 800
 kc2 = 0
 DO  i = 1,ncb
   ixr = i - 1
   DO  j = il2,ir2
     IF (parea(i,j,3) < 0.005) CYCLE
     jxr = j - 1
     kc2 = kc2 + 1
     IF (kc2 >= 125) GO TO 830
     xk2(kc2) = boxl*(FLOAT(ixr)+0.5)
     yk2(kc2) = boxw*FLOAT(jxr)
   END DO
 END DO
 DO  j = il2,ir2
   kc2 = kc2 + 1
   IF (kc2 >= 125) GO TO 830
   kte2(j)  = kc2
   xk2(kc2) = xwte(j)
   yk2(kc2) = ywte(j)
 END DO
 kc2t = kc2
 IF (iprint <= 0) GO TO 800
 WRITE  (n6,780) (i,xk2(i),i=1,kc2)
 780 FORMAT (1H1,25H cntrl2 box ctr. x coord.,/(10(1X,i4,f8.2)))
 WRITE  (n6,790) (i,yk2(i),i=1,kc2)
 790 FORMAT (1H1,25H cntrl2 box ctr. y coord.,/(10(1X,i4,f8.2)))
 800 CONTINUE
 boxl = boxl/cr
 boxw = boxw/cr
 boxa = boxa/cr**2
 DO  i = 1,12
   x(i) = x(i)/cr
   y(i) = y(i)/cr
 END DO
 DO  i = 1,50
   xwte(i) = xwte(i)/cr
   ywte(i) = ywte(i)/cr
 END DO
 GO TO 840
 830 ireg = 2
 840 RETURN
END SUBROUTINE mbreg
