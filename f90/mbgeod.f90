SUBROUTINE mbgeod
     
!     SUBROUTINE TO COMPUTE GEOMETRY AND INDEXES OF REGIONS
 
 LOGICAL :: cntrl2,cntrl1,crank1,crank2,asym
 COMMON /mboxc/ njj,crank1,crank2,cntrl1,cntrl2,nbox,npts0,npts1,  &
     npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,boxl,boxw,  &
     boxa,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 COMMON /mboxa/ x(12),y(12),tang(10),ang(10),cotang(10)
 
!     MAIN GEOMETRY
 
 big = -1.0E35
 DO  i = 1,10
   tang(i) = 0.0
   ang(i)  = 0.0
 END DO
 y(4) = y(1)
 y(6) = y(3)
 
 IF (crank1) GO TO 400
 x(2) = x(3)
 y(2) = y(3)
 tang(2) = 0.0
 
 400  IF (crank2) GO TO 500
 x(5) = x(6)
 y(5) = y(6)
 tang(5) = 0.0
 
 500  tang(1) = (x(2)-x(1))/(y(2)-y(1))
 ang(1)  = 57.2958*ATAN(tang(1))
 IF (crank1) tang(2) = (x(3)-x(2))/(y(3)-y(2))
 ang(2)  = 57.2958*ATAN(tang(2))
 tang(4) = (x(5)-x(4))/(y(5)-y(4))
 ang(4)  = 57.2958*ATAN(tang(4))
 IF (crank2) tang(5) = (x(6)-x(5))/(y(6)-y(5))
 ang(5)  = 57.2958*ATAN(tang(5))
 
 areaw   = 0.5*(x(1)*(y(1)-y(2)) + x(2)*(y(1)-y(3)) +  &
     x(3)*(y(2)-y(3)) + x(4)*(y(5)-y(1)) + x(5)*(y(3)-y(1)) + x(6)*(y(3)-y(5)))
 
!     CONTROL1 SURFACE GEOMETRY
 
 area1 = 0.0
 IF (.NOT.cntrl1) GO TO 1620
 tang(7) = (x(9)-x(8))/(y(9)-y(8))
 ang(7)  = 57.2958*ATAN(tang(7))
 
 IF (ABS(y(7)-y(8)) > 0.01) GO TO 1000
 y(7) = y(8)
 tm   = big
 IF (y(7) > y(5)) GO TO 900
 800  x(7) = tang(4)*(y(7)-y(4)) + x(4)
 GO TO 1100
 900  x(7) = tang(5)*(y(7)-y(5)) + x(5)
 GO TO 1100
 
 1000 tm   = (x(7)-x(8))/(y(7)-y(8))
 IF (y(5) == y(7) .AND. x(5) == x(7)) GO TO 1100
 y(7) = (tm*y(8)-tang(4)*y(4)+x(4)-x(8))/(tm-tang(4))
 IF (y(7) <= y(5)) GO TO 800
 y(7) = (tm*y(8)-tang(5)*y(5)+x(5)-x(8))/(tm-tang(5))
 GO TO 900
 1100 tang(6) = tm
 
 IF (ABS(y(11)-y(9)) > 0.01) GO TO 1400
 y(11) =  y(9)
 tm    =  big
 IF (y(11) > y(5)) GO TO 1300
 1200 x(11) = tang(4)*(y(11)-y(4)) + x(4)
 GO TO 1500
 1300 x(11) = tang(5)*(y(11)-y(5)) + x(5)
 GO TO 1500
 
 1400 tm    = (x(11)-x(9))/(y(11)-y(9))
 IF (y(5) == y(11) .AND. x(5) == x(11)) GO TO 1500
 y(11) = (tm*y(9)-tang(4)*y(4)+x(4)-x(9))/(tm-tang(4))
 IF (y(11) <= y(5)) GO TO 1200
 y(11) = (tm*y(9)-tang(5)*y(5)+x(5)-x(9))/(tm-tang(5))
 GO TO 1300
 1500 tang(8) = tm
 
 IF (y(7) <= y(5) .AND. y(11) >= y(5)) GO TO 1600
 area1 = 0.5*((x(8)-x(11))*(y(7)-y(9)) + (x(9)-x(7))*(y(8)-y(11)))
 GO TO 1620
 
 1600 area1 = 0.5*(x(5)*(y(11)-y(7)) + x(8)*(y(7)-y(9)) +  &
     x(9)*(y(8)-y(11)) + x(7)*(y(5)-y(8)) + x(11)*(y(9)-y(5)))
 
!     CONTROL2 SURFACE GEOMETRY
 
 1620 area2 = 0.0
 IF (.NOT.cntrl2) GO TO 1700
 tang(10) = (x(10)-x(9))/(y(10)-y(9))
 ang(10)  = 57.2958*ATAN(tang(10))
 IF (ABS(y(12)-y(10)) > 0.01) GO TO 1660
 y(12) = y(10)
 tm    = big
 IF (y(12) > y(5)) GO TO 1650
 1640 x(12) = tang(4)*(y(12)-y(4)) + x(4)
 GO TO 1670
 1650 x(12) = tang(5)*(y(12)-y(5)) + x(5)
 GO TO 1670
 1660 tm    = (x(12)-x(10))/(y(12)-y(10))
 IF (y(5) == y(12) .AND. x(5) == x(12) ) GO TO 1670
 y(12) = (tm*y(10)-tang(4)*y(4)+x(4)-x(10))/(tm-tang(4))
 IF (y(12) <= y(5)) GO TO 1640
 y(12) = (tm*y(10)-tang(5)*y(5)+x(5)-x(10))/(tm-tang(5))
 GO TO 1650
 1670 tang(9) = tm
 
 IF (y(11) <= y(5) .AND. y(12) >= y(5)) GO TO 1680
 area2 = 0.5*((x(9)-x(12))*(y(11)-y(10)) +(x(10)-x(11))*(y(9)-y(12)))
 GO TO 1700
 
 1680 area2 = 0.5*(x(5)*(y(12)-y(11))+x(9)*(y(11)  &
     - y(10))+x(10)*(y(9)-y(12))+x(11)*(y(5) - y(9))+x(12)*(y(10)-y(5)))
 
!     PRINT GEOMETRY DATA
 
 1700 cr    = x(4) - x(1)
 CALL mbprit (areaw,area1,area2)
 gc    = 2.0*cr**2
 xcent = (x(3)+x(4)+x(6))/4.0
 ycent = y(3)*(0.333 + 0.167*(x(6)-x(3))/x(4))
 
 DO  i = 1,10
   IF (tang(i) /= 0) GO TO 1900
   cotang(i) = big
   CYCLE
   1900 IF (tang(i) /= big) GO TO 2000
   cotang(i) = 0.
   CYCLE
   2000 cotang(i) = 1./tang(i)
 END DO
 RETURN
END SUBROUTINE mbgeod
