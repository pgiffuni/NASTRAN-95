SUBROUTINE mbctr (ictr,il1,ir1,ncn,nc1,nwn,nw1,parea)
     
!     CONTROL1 SURFACE
 
!     CONVEX ONLY:
!     MUST COMPILE WITH O1 OR LOWER OPTIMIZATION OPTION. IF O2 IS USED,
!     THE COMPILER WOULD GO INTO INFINITE LOOP.
 
 
 INTEGER, INTENT(IN OUT)                  :: ictr
 INTEGER, INTENT(OUT)                     :: il1
 INTEGER, INTENT(OUT)                     :: ir1
 INTEGER, INTENT(OUT)                     :: ncn(1)
 INTEGER, INTENT(OUT)                     :: nc1(1)
 INTEGER, INTENT(IN)                      :: nwn(1)
 INTEGER, INTENT(IN)                      :: nw1(1)
 REAL, INTENT(OUT)                        :: parea(50,50,3)
 LOGICAL :: cntrl2,cntrl1,crank1,crank2,asym
 DIMENSION  x(5),y(5),tang(5),cotang(5)
 COMMON /mboxa/ xx(12),yy(12),tg(10),ang(10),cotg(10)
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2,nbox,npts0,npts1,  &
     npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,boxl,boxw,  &
     boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 
 IF (ictr == 2) GO TO 1000
 x(1) = xx( 7)
 x(2) = xx( 8)
 x(3) = xx( 9)
 x(4) = xx(11)
 y(1) = yy( 7)
 y(2) = yy( 8)
 y(3) = yy( 9)
 y(4) = yy(11)
 tang(1) = tg(6)
 tang(2) = tg(7)
 tang(3) = tg(8)
 cotang(1) = cotg(6)
 cotang(2) = cotg(7)
 cotang(3) = cotg(8)
 GO TO 2000
 
!     CONTROL2 SURFACE
 
 1000 x(1) = xx(11)
 x(2) = xx( 9)
 x(3) = xx(10)
 x(4) = xx(12)
 y(1) = yy(11)
 y(2) = yy( 9)
 y(3) = yy(10)
 y(4) = yy(12)
 tang(1) = tg( 8)
 tang(2) = tg(10)
 tang(3) = tg( 9)
 cotang(1) = cotg( 8)
 cotang(2) = cotg(10)
 cotang(3) = cotg( 9)
 
 2000 x(5) = xx(5)
 y(5) = yy(5)
 tang(4) = tg(4)
 tang(5) = tg(5)
 cotang(4) = cotg(4)
 cotang(5) = cotg(5)
 
 il1 = AMIN1(y(2),y(1))/boxw + 1.5
 ir1 = AMAX1(y(3),y(4))/boxw + 1.4999
 DO  i = il1,ir1
   yr  = (FLOAT(i)-0.5)*boxw
   yl  = yr-boxw
   
   xll = (yl-y(2))*tang(1) + x(2)
   xrl = (yr-y(2))*tang(1) + x(2)
   xlh = (yl-y(2))*tang(2) + x(2)
   xrh = (yr-y(2))*tang(2) + x(2)
   xlr = (yl-y(3))*tang(3) + x(3)
   xrr = (yr-y(3))*tang(3) + x(3)
   
   IF (crank2 .AND. y(5) <= y(1)) GO TO 4515
   
   xlt = (yl-y(1))*tang(4) + x(1)
   xrt = (yr-y(1))*tang(4) + x(1)
   GO TO 4520
   
   4515 xlt = (yl-y(1))*tang(5) + x(1)
   xrt = (yr-y(1))*tang(5) + x(1)
   
   4520 IF (yl <= y(2) .AND. yr >= y(2)) GO TO 4525
   
   IF (yr < y(2)) GO TO 4550
   jt = (xlh-AMOD(xlh,boxl)+boxl)/boxl + 0.01
   GO TO 4600
   
   4525 jt = (x(2)-AMOD(x(2),boxl)+boxl)/boxl + 0.01
   GO TO 4600
   
   4550 jt = (xrl-AMOD(xrl,boxl)+boxl)/boxl + 0.01
   
   4600 IF (yl < y(4) .AND. yr >= y(4).AND.xrt >= xlt) GO TO 4625
   IF (yl >= y(4)) GO TO 4650
   
   jb = (AMAX1(xlt,xrt)-AMOD(AMAX1(xlt,xrt),boxl)+boxl)/boxl + 0.01
   GO TO 4700
   4625 jb = (x(4)-AMOD(x(4),boxl)+boxl)/boxl + 0.01
   GO TO 4700
   
   4650 jb = (xlr-AMOD(xlr,boxl)+boxl)/boxl + 0.01
   
   4700 DO  j = jt,jb
     
     xb = FLOAT(j)*boxl
     xt = xb - boxl
     
     ytl = (xt-x(2))*cotang(1) + y(2)
     ybl = (xb-x(2))*cotang(1) + y(2)
     yth = (xt-x(2))*cotang(2) + y(2)
     ybh = (xb-x(2))*cotang(2) + y(2)
     ytr = (xt-x(3))*cotang(3) + y(3)
     ybr = (xb-x(3))*cotang(3) + y(3)
     
     IF (crank2 .AND. y(5) <= y(1)) GO TO 4706
     
     ytt = (xt-x(1))*cotang(4) + y(1)
     ybt = (xb-x(1))*cotang(4) + y(1)
     GO TO 4708
     
     4706 ytt = (xt-x(1))*cotang(5) + y(1)
     ybt = (xb-x(1))*cotang(5) + y(1)
     
!     FULL BOXES
     
     4708 IF (yl >= ytl .AND. xt >= xrh .AND. yr < ybr .AND.  &
         xb < xrt .AND. xb < xlt) GO TO 4900
     
!     DOUBLE CORNER BOXES
     
     IF (yl <= y(2) .AND. yr >= y(2) .AND. xt < x(2) .AND.  &
         xb >= x(2) .AND. yl <= y(1) .AND. yr >= y(1) .AND.  &
         xt < x(1) .AND. xb >= x(1)) GO TO 4820
     
     IF (yl < y(3) .AND. yr >= y(3) .AND. xt < x(3) .AND.  &
         xb >= x(3) .AND. yl < y(4) .AND. yr >= y(4) .AND.  &
         xt < x(4) .AND. xb >= x(4)) GO TO 4840
     
!     SINGLE CORNER BOXES
     
     IF (yl <= y(2) .AND. yr >= y(2) .AND. xt < x(2) .AND.  &
         xb >= x(2)) GO TO 4710
     IF (yl <= y(1) .AND. yr >= y(1) .AND. xt < x(1) .AND.  &
         xb >= x(1)) GO TO 4730
     IF (yl < y(3) .AND. yr >= y(3) .AND. xt < x(3) .AND.  &
         xb >= x(3)) GO TO 4750
     IF (yl < y(4) .AND. yr >= y(4) .AND. xt < x(4) .AND.  &
         xb >= x(4)) GO TO 4770
     
!     HINGE + T. E. BOXES
     
     IF (xt < xrh .AND. (xb >= xlt .OR. xb >= xrt)) GO TO 4788
     
!     SIDE BOXES
     
     IF (xb >= xlh .AND. xt < xrh .AND. yl >= ytl .AND. (xb < x(3)  &
         .OR. yr < y(3))) GO TO 4745
     IF (yl <= ytl .AND. yr >= ybl .AND. xb >= x(2) .AND. xt < x(1))  &
         GO TO 4740
     IF (yl < ytr .AND. yr >= ybr .AND. xb >= x(3) .AND. xt < x(4)  &
         .AND. yr >= y(4)) GO TO 4765
     GO TO 4747
     
!     FWD LH CORNER
     
     4710 IF (yl <= ybl .AND. yr >= ybl .AND. xt < xrh .AND. xb >= xrh)  &
         pa = .5*((y(2)-ybl)*(xb-x(2))+(2.*xb-x(2)-xrh)*(yr-y(2)))/boxa
     IF (xt < xll .AND. xb >= xll .AND. xt < xrh .AND. xb >= xrh)  &
         pa = .5*((2.*xb-xll-x(2))*(y(2)-yl)+(2.*xb-xrh-x(2))*(yr-y(2))) /boxa
     IF (xt < xll .AND. xb >= xll .AND. yl < ybh .AND. yr >= ybh)  &
         pa = .5*((2.*xb-x(2)-xll)*(y(2)-yl)+(ybh-y(2))*(xb-x(2)))/boxa
     IF (yl <= ybl .AND. yr >= ybl .AND. yl < ybh .AND. yr >= ybh)  &
         pa = 0.5*(xb-x(2))*(ybh-ybl)/boxa
     IF (i-1 == 0) THEN
       GO TO  5000
     ELSE
       GO TO  4799
     END IF
     
     4720 IF (yl <= ytl .AND. yr >= ytl .AND. yl < yth .AND. yr >= yth .AND. &
             xt < xll .AND. xb >= xll .AND. xt < xrh .AND. xb >= xrh)  &
         pa = 1.0 - 0.5*((xll-xt)*(ytl-yl)+(yr-yth)*(xrh-xt))/boxa
     IF (yl <= ybl .AND. yr >= ybl .AND. yl <= ytl .AND. yr >= ytl .AND. &
             yl < yth .AND. yr >= yth .AND. xt < xrh .AND. xb >= xrh)  &
         pa = 1.0 - 0.5*((ytl+ybl-2.0*yl)*boxl+(yr-yth)*(xrh-xt))/boxa
     IF (yl < ybh .AND. yr >= ybh .AND. xt < xll .AND. xb >= xll .AND. &
             yl <= ytl .AND. yr >= ytl .AND. yl < yth .AND. yr >= yth)  &
         pa = 1.0 - 0.5*((ytl-yl)*(xll-xt)+(2.0*yr-yth-ybh)*boxl)/boxa
     IF (yl <= ytl .AND. yr >= ytl .AND. yl < yth .AND. yr >= yth .AND. &
             yl <= ybl .AND. yr >= ybl .AND. yl < ybh .AND. yr >= ybh)  &
         pa = 0.5*(yth+ybh-ytl-ybl)/boxw
     IF (yl <= ytl .AND. yr >= ytl .AND. yl <= ybl .AND. yr >= ybl .AND. &
             yl < ybt .AND. yr >= ybt .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*((2.0*yr-ytl-ybl)*boxl-(yr-ybt)*(xb-xrt))/boxa
     IF (yl < ytl .AND. yr >= ytl .AND. xt < xll .AND. xb >= xll .AND. &
             xt < xlt .AND. xb >= xlt .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 1.0 - 0.5*((ytl-yl)*(xll-xt)+(xb-xlt)*(ybt-yl))/boxa
     IF (xt < xll .AND. xb >= xll .AND. xt < xrl .AND. xb >= xrl .AND. &
             yl < ybt .AND. yr >= ybt .AND. xt < xlt .AND. xb >= xlt)  &
         pa = 0.5*((2.0*xb-xll-xrl)*boxw-(xb-xlt)*(ybt-yl))/boxa
     IF (yl < ytl .AND. yr >= ytl .AND. xt < xll .AND. xb >= xll .AND. &
             xt < xlt .AND. xb >= xlt .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*((xlt+xrt-2.0*xb)*boxw-(ytl-yl)*(xll-xt))/boxa
     IF (xt < xll .AND. xb >= xll .AND. xt < xlt .AND. xb >= xlt .AND. &
             xt < xrl .AND. xb >= xrl .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*(xlt+xrt-xll-xrl)/boxl
     IF (i-1 == 0) THEN
       GO TO  5000
     ELSE
       GO TO  4799
     END IF
     
!     AFT LH CORNER
     
     4730 IF (x(1) >= x(4)) GO TO 4735
     IF (yl < ybt .AND. yr >= ybt .AND. yl <= ytl .AND. yr >= ytl)  &
         pa = .5*((2.*yr-ytl-y(1))*(x(1)-xt)+(2.*yr-y(1)-ybt)*(xb-x(1))) /boxa
     4735 IF (yl <= ytl .AND. yr >= ytl .AND. xt < xrt .AND. xb >= xrt)  &
         pa = .5*((y(1)-ytl)*(x(1)-xt)+(x(1)+xrt-2.*xt)*(yr-y(1)))/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. yl <= ytl .AND. yr >= ytl .AND. &
             ytt >= ytl) pa = 0.5*(ytt-ytl)*(x(1)-xt)/boxa
     IF (xt < xrl .AND. xb >= xrl .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*(xrt-xrl)*(yr-y(1))/boxa
     IF (yl < ybt .AND. yr >= ybt .AND. xt < xrl .AND. xb >= xrl)  &
         pa = .5*((x(1)-xrl)*(yr-y(1))+(2.*yr-y(1)-ybt)*(xb-x(1)))/boxa
     IF (i-1 == 0) THEN
       GO TO  5000
     ELSE
       GO TO  4799
     END IF
     
!     LH EDGE
     
     4740 IF (i == 1) GO TO 4800
     
     IF (xt < xll .AND. xb >= xll .AND. xt < xrl .AND. xb >= xrl)  &
         pa = 0.5*(2.0*xb-xll-xrl)/boxl
     IF (yl <= ybl .AND. yr >= ybl .AND. xt < xrl .AND. xb >= xrl)  &
         pa = 0.5*(yr-ybl)*(xb-xrl)/boxa
     IF (yl <= ytl .AND. yr >= ytl .AND. yl <= ybl .AND. yr >= ybl)  &
         pa = 0.5*(2.0*yr-ytl-ybl)/boxw
     IF (yl <= ytl .AND. yr >= ytl .AND. xt < xll .AND. xb >= xll)  &
         pa = 1.0 - 0.5*(xll-xt)*(ytl-yl)/boxa
     GO TO 4720
     
!     HINGE LINE
     
     4745 IF (yl < yth .AND. yr >= yth .AND. xt < xrh .AND. xb >= xrh)  &
         pa = 1.0 - 0.5*(yr-yth)*(xrh-xt)/boxa
     IF (xt < xlh .AND. xb >= xlh .AND. xt < xrh .AND. xb >= xrh)  &
         pa = 0.5*(2.0*xb-xlh-xrh)/boxl
     IF (yl < ybh .AND. yr >= ybh .AND. xt < xlh .AND. xb >= xlh)  &
         pa = 0.5*(xb-xlh)*(ybh-yl)/boxa
     IF (yl < yth .AND. yr >= yth .AND. yl < ybh .AND. yr >= ybh)  &
         pa = 0.5*(yth+ybh-2.0*yl)/boxw
     GO TO 4760
     
!     TRAILING EDGE
     
     4747 IF (yl < ytt .AND. yr >= ytt .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*(yr-ytt)*(xrt-xt)/boxa
     IF (xt < xlt .AND. xb >= xlt .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*(xlt+xrt-2.0*xt)/boxl
     IF (xt < xlt .AND. xb >= xlt .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 1.0 - 0.5*(xb-xlt)*(ybt-yl)/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 0.5*(2.0*yr-ytt-ybt)/boxw
     IF (yl < ybt .AND. yr >= ybt .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 1.0 - 0.5*(yr-ybt)*(xb-xrt)/boxa
     IF (xt < xlt .AND. xb >= xlt .AND. yl < ytt .AND. yr >= ytt)  &
         pa = 0.5*(xlt-xt)*(ytt-yl)/boxa
     GO TO 4799
     
!     FWD RH CORNER
     
     4750 IF (yr >= ybr .AND. yl < ybr .AND. xt < xlh .AND. xb >= xlh)  &
         pa = 0.5*((x(3)-xlh)*(y(3)-yl)+(y(3)+ybr-2.0*yl)*(xb-x(3))) /boxa
     IF (xt < xlh .AND. xb >= xlh .AND. xt < xlr .AND. xb >= xlr)  &
         pa = 0.5*(xlr-xlh)*(y(3)-yl)/boxa
     IF (yr >= yth .AND. yl < yth .AND. xt < xlr .AND. xb >= xlr)  &
         pa = .5*((yth+y(3)-2.*yl)*(x(3)-xt)+(y(3)-yl)*(xlr-x(3)))/boxa
     IF (yr >= yth .AND. yl < yth .AND. yr >= ybr .AND. yl < ybr)  &
         pa = 0.5*((yth+y(3)-2.0*yl)*(x(3)-xt)+(y(3)+ybr-2.0*yl)* (xb-x(3)))/boxa
     GO TO 4799
     
     4760 IF (yr >= ybr .AND. yl < ybr .AND. xt < xrr .AND. xb >= xrr .AND. &
             xt < xrh .AND. xb >= xrh .AND. yl < yth .AND. yr >= yth)  &
         pa = 1.0 - 0.5*((yr-yth)*(xrh-xt)+(xb-xrr)*(yr-ybr))/boxa
     IF (yr >= ybr .AND. yl < ybr .AND. xt < xrr .AND. xb >= xrr .AND. &
             xt < xrh .AND. xb >= xrh .AND. xt < xlh .AND. xb >= xlh)  &
         pa = 0.5*((2.0*xb-xlh-xrh)*boxw-(yr-ybr)*(xb-xrr))/boxa
     IF (yl < yth .AND. yr >= yth .AND. xt < xrh .AND. xb >= xrh .AND. &
             xt < xrr .AND. xb >= xrr .AND. xt < xlr .AND. xb >= xlr)  &
         pa = 0.5*((xlr+xrr-2.0*xt)*boxw-(yr-yth)*(xrh-xt))/boxa
     IF (xt < xlh .AND. xb >= xlh .AND. xt < xlr .AND. xb >= xlr .AND. &
             xt < xrh .AND. xb >= xrh .AND. xt < xrr .AND. xb >= xrr)  &
         pa = 0.5*(xlr+xrr-xlh-xrh)/boxl
     GO TO 4799
     
!     RH EDGE
     
     4765 IF (xt < xrr .AND. xb >= xrr .AND. xt < xlr .AND. xb >= xlr)  &
         pa = 0.5*(xlr+xrr-2.0*xt)/boxl
     IF (yr >= ytr .AND. yl < ytr .AND. xt < xlr .AND. xb >= xlr)  &
         pa = 0.5*(ytr-yl)*(xlr-xt)/boxa
     IF (yr >= ytr .AND. yl < ytr .AND. yr >= ybr .AND. yl < ybr)  &
         pa = 0.5*(ytr+ybr-2.0*yl)/boxw
     IF (yr >= ybr .AND. yl < ybr .AND. xt < xrr .AND. xb >= xrr)  &
         pa = 1.0 - 0.5*(xb-xrr)*(yr-ybr)/boxa
     GO TO 4780
     
!     AFT RH CORNER
     
     4770 IF (x(4) >= x(1)) GO TO 4775
     IF (yr >= ytr .AND. yl < ytr .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 0.5*((ytr+y(4)-2.0*yl)*(x(4)-xt)+(y(4)+ybt-2.0*yl)* (xb-x(4)))/boxa
     4775 IF (yr >= ytr .AND. yl < ytr .AND. xt < xlt .AND. xb >= xlt)  &
         pa = .5*((xlt-x(4))*(y(4)-yl)+(ytr+y(4)-2.*yl)*(x(4)-xt))/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. yr >= ytr .AND. yl < ytr)  &
         pa = 0.5*(ytr-ytt)*(x(4)-xt)/boxa
     IF (xt < xlt .AND. xb >= xlt .AND. xt < xrr .AND. xb >= xrr)  &
         pa = 0.5*((xlt+x(4)-2.0*xt)*(y(4)-yl)+(x(4)+xrr-2.0*xt)* (yr-y(4)))/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. xt < xrr .AND. xb >= xrr)  &
         pa = .5*((y(4)-ytt)*(x(4)-xt)+(x(4)+xrr-2.*xt)*(yr-y(4)))/boxa
     GO TO 4799
     
     4780 IF (xt < xlt .AND. xb >= xlt .AND. yr >= ybr .AND. yl < ybr .AND. &
             yl < ybt .AND. yr >= ybt .AND. xt < xrr .AND. xb >= xrr)  &
         pa = 1.0 - 0.5*((xb-xlt)*(ybt-yl)+(yr-ybr)*(xb-xrr))/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. yl < ybt .AND. yr >= ybt .AND. &
             yr >= ybr .AND. yl < ybr .AND. xt < xrr .AND. xb >= xrr)  &
         pa = 0.5*((2.0*yr-ytt-ybt)*boxl-(yr-ybr)*(xb-xrr))/boxa
     IF (yr >= ytr .AND. yl < ytr .AND. yr >= ybr .AND. yl < ybr .AND. &
             xt < xlt .AND. xb >= xlt .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 0.5*((ytr+ybr-2.0*yl)*boxl-(xb-xlt)*(ybt-yl))/boxa
     IF (yl < ytt .AND. yr >= ytt .AND. yr >= ytr .AND. yl < ytr .AND. &
             yl < ybt .AND. yr >= ybt .AND. yr >= ybr .AND. yl < ybr)  &
         pa = 0.5*(ytr-ytt+ybr-ybt)/boxw
     GO TO 4799
     
     4788 IF (xb >= xlt .AND. yr >= ybt .AND. yl < yth)  &
         pa = 1.0 - 0.5*((yr-yth)*(xrh-xt)+(xb-xlt)*(ybt-yl))/boxa
     IF (xb >= xrt .AND. yl < ybt .AND. yl < yth)  &
         pa = 1.0 - 0.5*((yr-yth)*(xrh-xt)+(xb-xrt)*(yr-ybt))/boxa
     IF (xt < xlh .AND. xb >= xlt .AND. xb >= xrt)  &
         pa = 0.5*(xrt-xrh+xlt-xlh)/boxl
     IF (xt < xlh .AND. yl < ybt .AND. xb >= xrt)  &
         pa = 1.0 - 0.5*((xlh+xrh-2.0*xt)*boxw+(xb-xrt)*(yr-ybt))/boxa
     IF (yl < yth .AND. xb >= xlt .AND. xb >= xrt)  &
         pa = 1.0 - 0.5*((2.0*xb-xlt-xrt)*boxw+(xrh-xt)*(yr-yth))/boxa
     IF (xt < xlh .AND. yr >= ybt .AND. xb >= xlt)  &
         pa = 1.0 - 0.5*((xlh+xrh-2.0*xt)*boxw+(xb-xlt)*(ybt-yl))/boxa
     
     4799 parea(j,i,2) = pa
     parea(j,i,1) = parea(j,i,1) - pa
     CYCLE
     
     4800 yl1  = 0.0
     xll1 = (yl1-y(2))*tang(1) + x(2)
     
     IF (xt < xll1 .AND. xb >= xll1 .AND. xt < xrl .AND. xb >= xrl)  &
         pa = 0.5*(2.0*xb-xll1-xrl)*(yr-yl1)/boxa
     IF (yl1 <= ybl .AND. yr >= ybl .AND. xt < xrl .AND. xb >= xrl)  &
         pa = 0.5*(xb-xrl)*(yr-ybl)/boxa
     IF (yl1 <= ytl .AND. yr >= ytl .AND. xt < xll .AND. xb >= xll)  &
         pa = 1.0 - 0.5*(ytl-yl1)*(xll1-xt)/boxa
     IF (yl1 <= ytl .AND. yr >= ytl .AND. yl1 <= ybl .AND. yr >= ybl)  &
         pa = 0.5*(2.0*yr-ytl-ybl)/boxw
     GO TO 4720
     
!     LH CORNERS
     
     4820 IF (yl <= yth .AND. yr >= yth .AND. yl <= ybt .AND. yr >= ybt)  &
         pa = 0.5*((2.0*yr-yth-y(2))*(x(2)-xt)+(2.0*yr-y(2)-y(1))*  &
         (x(1)-x(2))+(2.0*yr-y(1)-ybt)*(xb-x(1)))/boxa
     IF (xt < xrh .AND. xb >= xrh .AND. yl <= ybt .AND. yr >= ybt)  &
         pa = 0.5*((x(2)-xrh)*(yr-y(2))+(2.0*yr-y(2)-y(1))*(x(1)-x(2))  &
         +(2.0*yr-y(1)-ybt)*(xb-x(1)))/boxa
     IF (yl <= yth .AND. yr >= yth .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*((2.0*yr-yth-y(2))*(x(2)-xt)+(2.0*yr-y(2)-y(1))*  &
         (x(1)-x(2))+(xrt-x(1))*(yr-y(1)))/boxa
     IF (xt < xrh .AND. xb >= xrh .AND. xt < xrt .AND. xb >= xrt)  &
         pa = 0.5*((x(2)-xrh)*(yr-y(2))+(2.0*yr-y(2)-y(1))*(x(1)-x(2))  &
         +(xrt-x(1))*(yr-y(1)))/boxa
     IF (i-1 == 0) THEN
       GO TO  5000
     ELSE
       GO TO  4799
     END IF
     
!     RH CORNERS
     
     4840 IF (yl < yth .AND. yr >= yth .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 0.5*((yth+y(3)-2.0*yl)*(x(3)-xt)+(y(3)+y(4)-2.0*yl)*(x(4)  &
         -x(3))+(y(4)+ybt-2.0*yl)*(xb-x(4)))/boxa
     IF (xt < xlh .AND. xb >= xlh .AND. yl < ybt .AND. yr >= ybt)  &
         pa = 0.5*((x(3)-xlh)*(y(3)-yl)+(y(3)+y(4)-2.0*yl)*(x(4)-x(3))  &
         +(y(4)+ybt-2.0*yl)*(xb-x(4)))/boxa
     IF (yl < yth .AND. yr >= yth .AND. xt < xlt .AND. xb >= xlt)  &
         pa = 0.5*((yth+y(3)-2.0*yl)*(x(3)-xt)+(y(3)+y(4)-2.0*yl)*(x(4)  &
         -x(3))+(xlt-x(4))*(y(4)-yl))/boxa
     IF (xt < xlh .AND. xb >= xlh .AND. xt < xlt .AND. xb >= xlt)  &
         pa = 0.5*((x(3)-xlh)*(y(3)-yl)+(y(3)+y(4)-2.0*yl)*(x(4)-x(3))  &
         +(xlt-x(4))*(y(4)-yl))/boxa
     GO TO 4799
     
     4900 parea(j,i,2) = 1.0
     parea(j,i,1) = 0.0
     CYCLE
     
     5000 parea(j,i,2) = 2.0*pa
     parea(j,i,1) = parea(j,i,1) - parea(j,i,2)
     
   END DO
   
   yc = yr - boxw/2.0
   xf = (yl-y(2))*tang(2) + x(2)
   
   IF (yc < y(2)) GO TO 5600
   IF (yc < y(3)) GO TO 5800
   IF (yc >= y(4)) CYCLE
   xf2 = (yr-y(3))*tang(3) + x(3)
   nc1(i) = xf2/boxl + 1.0
   IF (yc < y(1)) GO TO 5900
   5500 ncn(i) =  nwn(i)
   CYCLE
   5600 IF (yc < y(1)) CYCLE
   xf1 = (yr-y(2))*tang(1) + x(2)
   nc1(i) = xf1/boxl + 1.0
   5700 nc1(i) = MAX0(nc1(i),nw1(i))
   IF (yc < y(4)) GO TO 5500
   xf2 = (yr-y(3))*tang(3) + x(3)
   ncn(i) = xf2/boxl + 1.0
   CYCLE
   5800 nc1(i) = xf/boxl + 1.0
   IF (yc >= y(1)) GO TO 5700
   5900 xf1 = (yr-y(2))*tang(1) + x(2)
   ncn(i) = xf1/boxl + 1.0
 END DO
 
 RETURN
END SUBROUTINE mbctr
