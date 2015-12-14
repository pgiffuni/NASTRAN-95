SUBROUTINE pload1 (opt,islt,v,sa,sb,ba,bb,pa,pb,ta,tb,slt,ept)
     
!     PLOAD1 CALCULATES THE END LOADS ON A BAR ELEMENT FOR PLOAD1 LOADS
!     IT IS CALLED ONLY BY PLBAR1
 
!     OPT   = 1, CALLED FROM PLBAR1/EXTERN,  2, CALLED FROM SDRX
!     SLT   = PLOAD1 CARD
!     V     = REFERENCE VECTOR IN BASIC
!     SA    = OFFSET VECTOR IN BASIC POINT A
!     SB    = OFFSET VECTOR IN BASIC POINT B
!     BA    = BASIC  COORD  FOR POINT A
!     BB    = BASIC  COORD  FOR POINT B
!     PA    = LOAD   VECOTR FOR POINT A
!     PB    = LOAD   VECOTR FOR POINT B
!     TA,TB = TRANSFORMATION MATRICES FOR A AND B ONLY USED WITH OPT 1
!     EPT   = POINTER TO EST
 
 
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER, INTENT(IN)                      :: islt(7)
 REAL, INTENT(IN OUT)                     :: v(3)
 REAL, INTENT(IN)                         :: sa(3)
 REAL, INTENT(IN)                         :: sb(3)
 REAL, INTENT(IN)                         :: ba(3)
 REAL, INTENT(IN)                         :: bb(3)
 REAL, INTENT(OUT)                        :: pa(6)
 REAL, INTENT(OUT)                        :: pb(6)
 REAL, INTENT(IN OUT)                     :: ta(9)
 REAL, INTENT(IN OUT)                     :: tb(9)
 REAL, INTENT(IN)                         :: slt(7)
 REAL, INTENT(IN)                         :: ept(32)
 INTEGER :: oldid,TYPE,scale
 REAL :: LEN
 DOUBLE PRECISION :: ax,ay,az,bx,by,bz,dx1,dx2,dl,dt,dfx1,dfy1,dfz1,  &
     dfx2,dfy2,dfz2,s1,s2,s3,s4,s5, i01,i11,i21,i31,i41,i02,i12,i22,i32,i42
 DIMENSION  a(3),b(3),c(3),d(9),e(9), tp(3)
 COMMON /matout/  f,g
 EQUIVALENCE      (a(1),e(1)),(b(1),e(7)),(c(1),e(4))
 DATA    oldid ,  d / 10*0 /
 
 IF (oldid == islt(1)) GO TO 20
 oldid = islt(1)
 
!     CALCULATE AXIS AND LENGTH, AND THE E MATRIX
 
 a(1) = bb(1)-ba(1) + sb(1)-sa(1)
 a(2) = bb(2)-ba(2) + sb(2)-sa(2)
 a(3) = bb(3)-ba(3) + sb(3)-sa(3)
 LEN  = SQRT(sadotb(a,a))
 IF (LEN == 0.0) GO TO 380
 a(1) = a(1)/LEN
 a(2) = a(2)/LEN
 a(3) = a(3)/LEN
 CALL saxb (a,v,b)
 
 temp = SQRT(sadotb(b,b))
 IF (temp == 0.0) GO TO 380
 b(1) = b(1)/temp
 b(2) = b(2)/temp
 b(3) = b(3)/temp
 CALL saxb (b,a,c)
 
!     TRANSVERSE SHEAR
 
 temp = ept(31)*ept(17)*g*LEN**2
 tmp  = 12.0*f*ept(18)
 aly  = 0.0
 IF (ABS(temp+tmp) > 1.0E-14) aly = tmp/(tmp+temp)
 omaly = 1.0 - aly
 
 temp = (temp/ept(31))*ept(32)
 tmp  = 12.0*f*ept(19)
 alz  = 0.0
 IF (ABS(temp+tmp) > 1.0E-14) alz = tmp/(tmp+temp)
 omalz = 1.0 - alz
 
!     START BUILDING THE FORCES AND MOMENTS
 
 20 TYPE  = islt(2)
 scale = islt(3)
 x1  =  slt(4)
 f1  =  slt(5)
 x2  =  slt(6)
 f2  =  slt(7)
 i   = (TYPE-1)/6 + 1
 j   = MOD(TYPE,6)
 IF (j == 0) GO TO 60
 SELECT CASE ( j )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
   CASE (    5)
     GO TO 50
   CASE (    6)
     GO TO 60
 END SELECT
 30 fx1 = a(j)*f1
 fy1 = c(j)*f1
 fz1 = b(j)*f1
 fx2 = a(j)*f2
 fy2 = c(j)*f2
 fz2 = b(j)*f2
 GO TO 100
 40 fx1 = f1
 fy1 = 0.0
 fz1 = 0.0
 fx2 = f2
 fy2 = 0.0
 fz2 = 0.0
 GO TO 100
 50 fx1 = 0.0
 fy1 = f1
 fz1 = 0.0
 fx2 = 0.0
 fy2 = f2
 fz2 = 0.0
 GO TO 70
 60 fx1 = 0.0
 fy1 = 0.0
 fz1 = f1
 fx2 = 0.0
 fy2 = 0.0
 fz2 = f2
 70 j   = 4
 
!     SCALED
 
 100 IF (scale == 2 .OR. scale == 4) GO TO 110
 x1  = x1/LEN
 x2  = x2/LEN
 
!     DISTRIBUTED SCALED LOADS
 
 110 IF (x1 == x2) GO TO 220
 IF (scale <= 2 .OR. j == 4) GO TO 120
 fscale = SQRT(1.0-a(j)**2)
 fx1  = fscale*fx1
 fy1  = fscale*fy1
 fz1  = fscale*fz1
 fx2  = fscale*fx2
 fy2  = fscale*fy2
 fz2  = fscale*fz2
 
!     DISTRIBUTED LOADS
 
 120 dx1  = x1
 dx2  = x2
 dl   = LEN
 dfx1 = fx1
 dfy1 = fy1
 dfz1 = fz1
 dfx2 = fx2
 dfy2 = fy2
 dfz2 = fz2
 s1   = dx2 - dx1
 s2   = .5000000D0*(dx2**2 - dx1**2)
 s3   = .3333333D0*(dx2**3 - dx1**3)
 s4   = .2500000D0*(dx2**4 - dx1**4)
 s5   = .2000000D0*(dx2**5 - dx1**5)
 IF (i == 2) GO TO 140
 
!     FORCES
 
 i01  = dl*(s1-s2)
 i11  = dl*(s1-3.0D0*s3 + 2.0D0*s4)
 i21  = dl*(   3.0D0*s3 - 2.0D0*s4)
 i31  = dl*(s2-2.0D0*s3 + s4)
 i41  = dl*(s4-s3)
 dt   = dl*dl
 i02  = dt*(s2-s3)
 IF (f1 == f2) GO TO 130
 i12  = dt*(s2-3.0D0*s4 + 2.0D0*s5)
 i22  = dt*(   3.0D0*s4 - 2.0D0*s5)
 i32  = dt*(s3-2.0D0*s4 + s5)
 i42  = dt*(s5-s4)
 dt   = dl*(dx2-dx1)
 bx   = (dfx2-dfx1)/dt
 by   = (dfy2-dfy1)/dt
 bz   = (dfz2-dfz1)/dt
 ax   = dfx1 - dx1*bx*dl
 ay   = dfy1 - dx1*by*dl
 az   = dfz1 - dx1*bz*dl
 GO TO 170
 130 ax   = dfx1
 ay   = dfy1
 az   = dfz1
 GO TO 160
 
!     MOMENTS
 
 140 i01  = dl*(s1-s2)
 i11  =-6.0D0*(s2-s3)
 i21  =-i11
 i31  = s1 - 4.0D0*s2 + 3.0D0*s3
 i41  =    - 2.0D0*s2 + 3.0D0*s3
 IF (f1 == f2) GO TO 150
 i02  = (s2-s3)*dl**2
 i12  =-6.0D0*dl*(s3-s4)
 i22  =-i12
 i32  = dl*(s2-4.0D0*s3 + 3.0D0*s4)
 i42  =-dl*(   2.0D0*s3 - 3.0D0*s4)
 dt   = (dx2 -dx1)*dl
 bx   = (dfx2-dfx1)/dt
 by   = (dfz2-dfz1)/dt
 bz   =-(dfy2-dfy1)/dt
 ax   = dfx1 + dx1*bx*dl
 ay   = dfz1 + dx1*by*dl
 az   =-dfy1 + dx1*bz*dl
 GO TO 170
 150 ax   = dfx1
 ay   = dfz1
 az   =-dfy1
 160 bx   = 0.0D0
 by   = 0.0D0
 bz   = 0.0D0
 i12  = 0.0D0
 i22  = 0.0D0
 i32  = 0.0D0
 i42  = 0.0D0
 
!     LOADS
 
 170 pa(1) = i01*ax + i02*bx
 pa(2) = i11*ay + i12*by
 pa(3) = i11*az + i22*bz
 pa(4) = 0.0
 pa(5) =-dl*(i31*az + i32*bz)
 pa(6) = dl*(i31*ay + i32*by)
 dt    = dl*dl
 pb(1) = dl*s2*ax + dt*s3*bx
 pb(2) = i21  *ay + i22  *by
 pb(3) = i21  *az + i22  *bz
 pb(4) = 0.0
 pb(5) =-dl*(i41*az + i42*bz)
 pb(6) = dl*(i41*ay + i42*by)
 IF (i == 2) GO TO 190
 IF (aly == 0.0) GO TO 180
 pa(2) = omaly*pa(2) + aly*(  i01*ay + i02*by   )
 pa(6) = omaly*pa(6) + aly*(  i02*ay - i41*by*dt)*.50
 pb(2) = omaly*pb(2) + aly*(dl*s2*ay + s3 *by*dt)
 pb(6) = omaly*pb(6) - aly*(  i02*ay - i41*by*dt)*.50
 180 IF (alz == 0.0) GO TO 300
 pa(3) = omalz*pa(3) + alz*(  i01*az + i02*bz   )
 pa(5) = omalz*pa(5) - alz*(  i02*az - i41*bz*dt)*.50
 pb(3) = omalz*pb(3) + alz*(dl*s2*az + s3 *bz*dt)
 pb(5) = omalz*pb(5) + alz*(  i02*az - i41*bz*dt)*.50
 GO TO 300
 190 temp  = pa(1)
 pa(1) = pa(4)
 pa(4) = temp
 temp  = pb(1)
 pb(1) = pb(4)
 pb(4) = temp
 IF (aly == 0.0) GO TO 200
 pa(2) = omaly*pa(2)
 pa(6) = omaly*pa(6) + aly*(i01*ay + i02*by)
 pb(2) = omaly*pb(2)
 pb(6) = omaly*pb(6) + aly*(dl*s2*ay + s3*by*dt)
 200 IF (alz == 0.0) GO TO 300
 pa(3) = omalz*pa(3)
 pa(5) = omalz*pa(5) + alz*(i01*az + i02*bz)
 pb(3) = omalz*pb(3)
 pb(5) = omalz*pb(5) + alz*(dl*s2*az + s3*bz*dt)
 GO TO 300
 
!     CONCENTRATED LOADS
 
 220 tmp = 1.0 - x1
 IF (i == 2) GO TO 230
 
!     FORCES
 
 temp  = 1.0 - 3.0*x1**2 + 2.0*x1**3
 pa(1) = tmp*fx1
 pa(2) = temp*fy1*omaly + fy1*tmp*aly
 pa(3) = temp*fz1*omalz + fz1*tmp*alz
 pa(4) = 0.0
 temp  =-LEN*x1*tmp**2
 tmp   = tmp*LEN*x1*.50
 pa(5) = temp*fz1*omalz - fz1*tmp*alz
 pa(6) =-temp*fy1*omaly + fy1*tmp*aly
 temp  = 3.0*x1**2 - 2.0*x1**3
 pb(1) = x1*fx1
 pb(2) = temp*fy1*omaly + fy1*x1*aly
 pb(3) = temp*fz1*omalz + fz1*x1*alz
 pb(4) = 0.0
 temp  = (1.0-x1)*LEN*x1**2
 pb(5) = temp*fz1*omalz + fz1*tmp*alz
 pb(6) =-temp*fy1*omaly - fy1*tmp*aly
 GO TO 300
 
!     MOMENTS
 
 230 temp  =-(6.0/LEN*x1)*tmp
 pa(1) = 0.0
 pa(2) = temp*fz1*omaly
 pa(3) =-temp*fy1*omalz
 pa(4) = tmp*fx1
 temp  = 1.0 - 4.0*x1 + 3.0*x1**2
 pa(5) = temp*fy1*omalz + fy1*tmp*alz
 pa(6) = temp*fz1*omaly + fz1*tmp*aly
 pb(1) = 0.0
 pb(2) =-pa(2)
 pb(3) =-pa(3)
 pb(4) = x1*fx1
 temp  = 3.0*x1**2 - 2.0*x1
 pb(5) = temp*fy1*omalz + fy1*x1*alz
 pb(6) = temp*fz1*omaly + fz1*x1*aly
 GO TO 300
 
!     PIN FLAGS
 
 300 CALL ploapf (ept,ept,LEN,pa,pb)
 
!     LOAD VECTORS DONE FOR SDRX
 
 IF (opt == 2) GO TO 400
 
!     TRANSFORM LOAD VECTOR TO GLOBAL
 
 CALL gmmats (e ,3,3,1,pa(1),3,1,0,tp   )
 CALL gmmats (ta,3,3,1,tp   ,3,1,0,pa(1))
 CALL gmmats (e ,3,3,1,pb(1),3,1,0,tp   )
 CALL gmmats (tb,3,3,1,tp   ,3,1,0,pb(1))
 CALL gmmats (e ,3,3,1,pa(4),3,1,0,tp   )
 CALL gmmats (ta,3,3,1,tp   ,3,1,0,pa(4))
 CALL gmmats (e ,3,3,1,pb(4),3,1,0,tp   )
 CALL gmmats (tb,3,3,1,tp   ,3,1,0,pb(4))
 
 DO  i = 1,3
   IF (sa(i) /= 0.0) GO TO 320
 END DO
 GO TO 330
 320 d(2)  =-sa(3)
 d(3)  = sa(2)
 d(4)  = sa(3)
 d(6)  =-sa(1)
 d(7)  =-sa(2)
 d(8)  = sa(1)
 CALL gmmats (d,3,3,0,pa(1),3,1,0,tp)
 pa(4) = pa(4) + tp(1)
 pa(5) = pa(5) + tp(2)
 pa(6) = pa(6) + tp(3)
 
 330 DO  i = 1,3
   IF (sb(i) /= 0.0) GO TO 350
 END DO
 GO TO 400
 350 d(2)  =-sb(3)
 d(3)  = sb(2)
 d(4)  = sb(3)
 d(6)  =-sb(1)
 d(7)  =-sb(2)
 d(8)  = sb(1)
 CALL gmmats (d,3,3,0,pb(1),3,1,0,tp)
 pb(4) = pb(4) + tp(1)
 pb(5) = pb(5) + tp(2)
 pb(6) = pb(6) + tp(3)
 GO TO 400
 
!     ERROR
 
 380 CALL mesage (-30,31,oldid)
 
 400 RETURN
END SUBROUTINE pload1
