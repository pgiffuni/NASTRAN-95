SUBROUTINE tspl1s (ts1,ts2,ts6,ts6s,ts7,ktr3,ktr31)
     
!    TRANSVERSE SHEAR ROUTINE1 FOR CTRPLT1 - SINGLE PRECISION VERSION
 
 
 REAL, INTENT(OUT)                        :: ts1(60)
 REAL, INTENT(IN OUT)                     :: ts2(60)
 REAL, INTENT(IN OUT)                     :: ts6(40)
 REAL, INTENT(IN OUT)                     :: ts6s(40)
 REAL, INTENT(IN OUT)                     :: ts7(60)
 REAL, INTENT(OUT)                        :: ktr3(400)
 REAL, INTENT(IN)                         :: ktr31(400)
 
 REAL :: j11,j12,j22
 DIMENSION  , gs1(4),ge1(9),be(7),ga(7),wt(7),cons(2)
 COMMON /sma1io/ x,y,z,dista,distb,distc,a1,a2,a3
 COMMON /matout/ em(6),dum6(9),rj11,rj12,rj22
 
 DATA  be          /0.33333333333333E0,0.47014206E0,  &
     0.05971588E0,0.47014206E0,0.101286505E0,0.79742699E0,  &
     0.101286505E0/, ga          /0.33333333333333E0,  &
     2*0.47014206E0,0.05971588E0,2*0.101286505E0,0.79742699E0/,  &
     wt          /0.1125E0,3*0.066197075E0,3*0.06296959E0/
 
 cons(1)=dista*distc
 cons(2)=distb*distc
 DO  i=1,60
   ts1(i)=0.0E0
 END DO
 DO  k=1,7
   DO  kase=1,2
     IF (kase == 1)  x=be(k)*dista
     IF (kase == 2) x=-be(k)*distb
     y=ga(k)*distc
     CALL tspl3s (ts6)
     cons1=wt(k)*cons(kase)
     thk=a1+a2*x+a3*y
     cons14=cons1*thk
     gs1(1)=rj11*cons14
     gs1(2)=rj12*cons14
     gs1(3)=gs1(2)
     gs1(4)=rj22*cons14
     thk1=thk**3/12.0E0
     cons11=cons1*thk1
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
     a31=a14+2.0*a13
     a32=a12+2.0*a16
     a33=a24+2.0*a23
     a34=a22+2.0*a26
     a35=a33+a11
     a36=a34+a31
     a37=a25+a32
     ge1(1)=em(1)*cons11
     ge1(2)=em(2)*cons11
     ge1(3)=em(3)*cons11
     ge1(4)=ge1(2)
     ge1(5)=em(4)*cons11
     ge1(6)=em(5)*cons11
     ge1(7)=ge1(3)
     ge1(8)=ge1(6)
     ge1(9)=em(6)*cons11
     
!        (B1) REFERS TO BENDING STRAIN DUE TO SECOND DERIVATIVES OF W
!        (B2) REFERS TO BENDING STRAINS DUE TO TRANSVERSE SHEAR STRAIN
!        (GAMMA) TRANSPOSE (GS) * (GAMMA) IS CONTRIBUTION OF STIFFNESS
!        MATRIX DUE TO WORK DONE BY SHEARING FORCES UNDERGOING SHEAR DEF
     
     
!  GAMMA TRANSPOSE GS GAMMA
     
     CALL gmmats (ts6,2,20,+1,gs1,2,2,0,ts6s)
     CALL gmmats (ts6s,20,2,-2,ts6,2,20,0,ktr3)
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
     
!  B2 TRANSPOSE D B2
     
     CALL gmmats (ts1,20,3,0,ge1,3,3,0,ts2)
     CALL gmmats (ts2,20,3,-2,ts1,20,3,+1,ktr3)
     
!  B2 TRANSPOSE D B1
     
     CALL tspl2s (ts7)
     CALL gmmats (ts2,20,3, 0,ts7,3,20, 0,ktr31)
     
!  B1 TRANSPOSE D B2
     
     DO  i=1,20
       DO  j=1,20
         ij=(i-1)*20+j
         ji=(j-1)*20+i
         ktr3(ij)=ktr3(ij)+ktr31(ij)+ktr31(ji)
       END DO
     END DO
   END DO
 END DO
 RETURN
END SUBROUTINE tspl1s
