SUBROUTINE strslv (ts6,nots)
     
!    STRESS ROUTINE ,CALLED FROM STRP11, FOR HIGHER ORDER PLATE ELEMENT
 
 
 REAL, INTENT(OUT)                        :: ts6(40)
 LOGICAL, INTENT(IN OUT)                  :: nots
 REAL :: j11,j12,j22
 
 
 COMMON /matout/ em(6)
 COMMON /sdr2x5/ dumsd(1324) ,    x,y,z,dista,distb,distc,a1,a2,a3,b1,b2,b3
 
 DO  i=1,40
   ts6(i)=0.0
 END DO
 thk=a1+a2*x+a3*y
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
 IF (nots) GO TO 146
 thk=b1+b2*x+b3*y
 j11=1.0/(em(6)*thk)
 j12=0.0
 j22=j11
 GO TO 148
 146 CONTINUE
 j11=1.0
 j12=0.0
 j22=1.0
 148 CONTINUE
 
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
 
 x2=x*x
 xy=x*y
 y2=y*y
 a38=a13+a14
 a39=a12+a16
 a40=a23+a24
 a41=a22+a26
 ts6( 7)=6.0*a11
 ts6( 8)=2.0*a31
 ts6( 9)=2.0*a32
 ts6(10)=6.0*a15
 ts6(11)=24.0*a11*x
 ts6(12)=6.0*(a31*x+a11*y)
 ts6(13)=4.0*(a32*x+a31*y)
 ts6(14)=6.0*(a15*x+a32*y)
 ts6(15)=24.0*a15*y
 IF (nots) GO TO 156
 ts6(16)=120.0*(-a11*a11-a13*a21+0.5*a11*x2)
 ts6(17)=12.0*(-a11*a32-a13*a34-a38*a31-a39*a33-a16*a11-a15*a21)  &
     +6.0*(a32*x2+2.0*a31*xy+a11*y2)
 ts6(18)=12.0*(-a11*a15-a13*a25-a38*a32-a39*a34-a16*a31-a15*a33)  &
     +6.0*(a15*x2+2.0*a32*xy+a31*y2)
 ts6(19)=24.0*(-a39*a25-a16*a32-a15*a34+a15*xy+0.5*a32*y2-a38*a15)
 ts6(20)=-120.0*(a16*a15+a15*a25-0.5*a15*y2)
 GO TO 158
 156 CONTINUE
 ts6(16)=60.0*a11*x2
 ts6(17)=6.0*(a32*x2+2.0*a31*xy+a11*y2)
 ts6(18)=6.0*(a15*x2+2.0*a32*xy+a31*y2)
 ts6(19)=12.0*(2.0*a15*xy+a32*y2)
 ts6(20)=60.0*a15*y2
 158 CONTINUE
 ts6(27)=6.0*a21
 ts6(28)=2.0*a33
 ts6(29)=2.0*a34
 ts6(30)=6.0*a25
 ts6(31)=24.0*a21*x
 ts6(32)=6.0*(a33*x+a21*y)
 ts6(33)=4.0*(a34*x+a33*y)
 ts6(34)=6.0*(a25*x+a34*y)
 ts6(35)=24.0*a25*y
 IF (nots) GO TO 166
 ts6(36)=120.0*(-a21*a11-a23*a21+0.5*a21*x2)
 ts6(37)=12.0*(-a21*a32-a23*a34-a40*a31-a41*a33-a26*a11-a25*a21)  &
     +6.0*(a34*x2+2.0*a33*xy+a21*y2)
 ts6(38)=12.0*(-a21*a15-a23*a25-a40*a32-a41*a34-a26*a31-a25*a33)  &
     +6.0*(a25*x2+2.0*a34*xy+a33*y2)
 ts6(39)=24.0*(-a41*a25-a26*a32-a25*a34+a25*xy+0.5*a34*y2-a40*a15)
 ts6(40)=-120.0*(a26*a15+a25*a25-0.5*a25*y2)
 GO TO 168
 166 CONTINUE
 ts6(36)=60.0*a21*x2
 ts6(37)=6.0*(a34*x2+2.0*a33*xy+a21*y2)
 ts6(38)=6.0*(a25*x2+2.0*a34*xy+a33*y2)
 ts6(39)=12.0*(2.0*a25*xy+a34*y2)
 ts6(40)=60.0*a25*y2
 168 CONTINUE
 RETURN
END SUBROUTINE strslv
