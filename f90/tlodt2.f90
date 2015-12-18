SUBROUTINE tlodt2 (ts1,ts2)
     
!    CALCULATION OF PTGEN2 - GEN THERMAL LOAD VECTOR DUE TO TRANSVERSE S
 
 
 REAL, INTENT(OUT)                        :: ts1(60)
 REAL, INTENT(OUT)                        :: ts2(20)
 DIMENSION  ts3(20)
 REAL :: j11,j12,j22
 COMMON /emgest/ est(100)
 COMMON /ssgwrk/ x,y,z,dista,distb,distc,a1,a2,a3,b1,b2,b3,g1(3), d(3)
 COMMON /matout/ em(6),dum6(9),rj11,rj12,rj22
 DIMENSION       be(7),ga(7),wt(7),cons(2)
 DATA       be / 0.33333333333333, 0.47014206   , 0.05971588  ,  &
     0.47014206      , 0.101286505  , 0.79742699  , 0.101286505     /
 DATA       ga / 0.33333333333333, 2*.47014206  , 0.05971588  ,  &
     2*0.101286505   , 0.79742699   /
 DATA       wt / 0.1125          , 3*0.066197075, 3*0.06296959/
 
 cons(1)=dista*distc
 cons(2)=distb*distc
 DO  i=1,60
   ts1(i)=0.0
 END DO
 DO  i=1,20
   ts3(i)=0.0
 END DO
 DO  k=1,7
   DO  kase=1,2
     IF (kase == 1) x= be(k)*dista
     IF (kase == 2) x=-be(k)*distb
     y=ga(k)*distc
     cons1=wt(k)*cons(kase)
     thk=a1+a2*x+a3*y
     temp=d(1)+d(2)*x+d(3)*y
     thk1=thk**3/12.0
     d11=em(1)*thk1
     d12=em(2)*thk1
     d13=em(3)*thk1
     d22=em(4)*thk1
     d23=em(5)*thk1
     d33=em(6)*thk1
     d21=d12
     d31=d13
     d32=d23
     j11=1.0/(em(6)*thk)
     j22=j11
     j12=0.0
     a11=-(j11*d11+j12*d13)
     a12=-(j11*d12+j12*d23)
     a13=-(j11*d13+j12*d33)
     a14=-(j11*d31+j12*d21)
     a15=-(j11*d32+j12*d22)
     a16=-(j11*d33+j12*d23)
     a21=-(j12*d11+j22*d13)
     a22=-(j12*d12+j22*d23)
     a23=-(j12*d13+j22*d33)
     a24=-(j12*d13+j22*d12)
     a25=-(j12*d23+j22*d22)
     a26=-(j12*d33+j22*d32)
     a31= a14+2.0*a13
     a32= a12+2.0*a16
     a33= a24+2.0*a23
     a34= a22+2.0*a26
     a35= a33+a11
     a36= a34+a31
     a37= a25+a32
     ts1(31)  =-24.0*a11
     ts1(33)  =-24.0*a21
     ts1(34)  =-6.0*a31
     ts1(35)  =-6.0*a21
     ts1(36)  =-6.0*a35
     ts1(37)  =-4.0*a32
     ts1(38)  =-4.0*a33
     ts1(39)  =-4.0*a36
     ts1(40)  =-6.0*a15
     ts1(41)  =-6.0*a34
     ts1(42)  =-6.0*a37
     ts1(44)  =-24.0*a25
     ts1(45)  =-24.0*a15
     ts1(46)  =-120.0*a11*x
     ts1(48)  =-120.0*a21*x
     ts1(49)  =-12.0*(a32*x+a31*y)
     ts1(50)  =-12.0*(a33*x+a21*y)
     ts1(51)  =-12.0*(a36*x+a35*y)
     ts1(52)  =-12.0*(a15*x+a32*y)
     ts1(53)  =-12.0*(a34*x+a33*y)
     ts1(54)  =-12.0*(a37*x+a36*y)
     ts1(55)  =-24.0*a15*y
     ts1(56)  =-24.0*(a25*x+a34*y)
     ts1(57)  =-24.0*(a15*x+a37*y)
     ts1(59)  =-120.0*a25*y
     ts1(60)  =-120.0*a15*y
     
     
     CALL gmmats (ts1,20,3,0,g_1,3,1,0,ts2)
     DO  i=1,20
       ts2(i)=ts2(i)*temp*thk1*cons1
       ts3(i)=ts3(i)+ts2(i)
     END DO
   END DO
 END DO
 DO  i=1,20
   ts2(i)=ts3(i)
 END DO
 RETURN
END SUBROUTINE tlodt2
