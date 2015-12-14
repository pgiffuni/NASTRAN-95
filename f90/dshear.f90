SUBROUTINE dshear
     
!     THIS COMPUTES THE THE TWO 6X6 DIFFERENTIAL STIFFNESS MATRICES
!     K(NPVT,NPVT) AND K(NPVT,J) WHERE J = 3,4,1,2 IF NPVT = 1,2,3,4
!     RESPECTIVELY.
 
!     ECPT FOR BOTH PANELS
 
!     ECPT( 1)  -  IELID          ELEMENT ID. NO.
!     ECPT( 2)  -  ISILNO(4)      SCALAR INDEX NUMBERS
!     ECPT( 3)  -   ...                   ...
!     ECPT( 4)  -   ...                   ...
!     ECPT( 5)  -   ...                   ...
!     ECPT( 6)  -  MATID          MATERIAL ID.
!     ECPT( 7)  -  T              THICKNESS
!     ECPT( 8)  -  FMU            NON-STRUCTURAL MASS
!     ECPT( 9)  -  ICSID1         COOR. SYS. ID. FOR GRID POINT 1
!     ECPT(10)  -  GP1(3)         BASIC COORDINATES FOR GRID POINT 1
!     ECPT(11)  -   ...                      ...
!     ECPT(12)  -   ...                      ...
!     ECPT(13)  -  ICSID2         COOR. SYS. ID. FOR GRID POINT 2
!     ECPT(14)  -  GP2(3)         BASIC COORDINATES FOR GRID POINT 2
!     ECPT(15)  -   ...                      ...
!     ECPT(16)  -   ...                      ...
!     ECPT(17)  -  ICSID3         COOR. SYS. ID. FOR GRID POINT 3
!     ECPT(18)  -  GP3(3)         BASIC COORDINATES FOR GRID POINT 3
!     ECPT(19)  -   ...                      ...
!     ECPT(20)  -   ...                      ...
!     ECPT(21)  -  ICSID4         COOR. SYS. ID. FOR GRID POINT 4
!     ECPT(22)  -  GP4(3)         BASIC COORDINATES FOR GRID POINT 4
!     ECPT(23)  -   ...                      ...
!     ECPT(24)  -   ...                      ...
!     ECPT(25)  -  TEMPEL         ELEMENT TEMPERATURE
!     ECPT(26)  -  DEFORM         ELEMENT DEFORMATION (NOT USED)
!     ECPT(27)  -  AVGLTP         AVG.ELEM LOADING TEMPERATURE, NOT USED
!     ECPT(28)  -  U1(3)          TRANSLATION DISPLACEMENTS AT PT. 1
!     ECPT(29)  -  ...                         ...
!     ECPT(30)  -  ...                         ...
!     ECPT(31)  -  U2(3)          TRANSLATION DISPLACEMENTS AT PT. 2
!     ECPT(32)  -  ...                         ...
!     ECPT(33)  -  ...                         ...
!     ECPT(34)  -  U3(3)          TRANSLATION DISPLACEMENTS AT PT. 3
!     ECPT(35)  -  ...                         ...
!     ECPT(36)  -  ...                         ...
!     ECPT(37)  -  U4(3)          TRANSLATION DISPLACEMENTS AT PT. 4
!     ECPT(38)  -  ...                         ...
!     ECPT(39)  -  ...                         ...
 
 REAL :: nusp
 DOUBLE PRECISION :: ke(36),ti(9),vleft(6),vd1,vd2,vkn,vk,v12,v41,  &
     vp12,vi,vj,avec,smallu,smallv,p,x1,x2,x3,x4,  &
     y1,y2,y3,y4,vkl,pa,v12dk,cep1,cep2,ep,temp
 DOUBLE PRECISION :: yp,xp,sa,xq,b,xl,a,a2,a3,a4,a5,b2,b3,b4,b5,c,c2,  &
     c3,c4,c5,d,d2,d3,d4,d5,term1,term2,term3,term4, term5,xl13,xl24
 DOUBLE PRECISION :: vp12l,vjl,z,term,f,e,g,nu,t,c23,nuc
 DOUBLE PRECISION :: ui(3),dpterm,sum,f13,f24,fxx,jj(3),j3x3(9), k3x3(9)
 DIMENSION        vd1(3),vd2(3),vkn(3),vk(3),v12(3),v41(3),vp12(3),  &
     vi(3),vj(3),avec(4),smallu(4),smallv(4),p(4), iecpt(100),ecpt(100),iz(1)
 COMMON /zzzzzz/  zz(1)
 COMMON /ds1aaa/  npvt,icstm,ncstm,dumcl(32),nogo
 COMMON /ds1aet/  ielid,isilno(4),matid,tsp,fmu,icsid1,gp1(3),  &
     icsid2,gp2(3),icsid3,gp3(3),icsid4,gp4(3),tempel,  &
     deform,avgltp,u1(3),u2(3),u3(3),u4(4)
 COMMON /ds1adp/  ke,ti,vleft,vd1,vd2,vkn,vk,v12,v41,vp12,vi,vj,  &
     avec,smallu,smallv,p,x1,x2,x3,x4,y1,y2,y3,y4, vkl,pa,v12dk,cep1,cep2,ep,temp
 COMMON /ds1adp/  yp,xp,sa,xq,b,xl,a,a2,a3,a4,a5,b2,b3,b4,b5,c,c2,  &
     c3,c4,c5,d,d2,d3,d4,d5,term1,term2,term3,term4, term5,xl13,xl24
 COMMON /ds1adp/  vp12l,vjl,z,term,f,e,g,nu,t,c23,nuc
 COMMON /ds1adp/  ui,dpterm,sum,f13,f24,fxx,jj,j3x3,k3x3
 COMMON /matin /  matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/  esp,gsp,nusp,rho,alpha,tsubo,gsube,sigt,sigc,sigs
 EQUIVALENCE      (iz(1),zz(1)),(ielid,iecpt(1),ecpt(1))
 
!     CALL MAT TO GET MATERIAL PROPERTIES.
 
 matidc = matid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 
!     STORE ECPT AND MPT VARIABLES IN DOUBLE PRECISION LOCATIONS
 
 e    = esp
 g    = gsp
 nu   = nusp
 t    = tsp
 c23  = 2.0D0/3.0D0
 nuc  = 1.0D0/(1.0D0+nu)
 
!     COMPUTE DIAGONAL VECTORS.
 
 DO  i = 1,3
   vd1(i) = gp3(i) - gp1(i)
   vd2(i) = gp4(i) - gp2(i)
 END DO
 
!     COMPUTE THE NORMAL VECTOR VKN, NORMALIZE, AND COMPUTE THE
!     PROJECTED AREA, PA
 
 vkn(1) = vd1(2)*vd2(3) - vd1(3)*vd2(2)
 vkn(2) = vd1(3)*vd2(1) - vd1(1)*vd2(3)
 vkn(3) = vd1(1)*vd2(2) - vd1(2)*vd2(1)
 vkl    = DSQRT(vkn(1)**2 + vkn(2)**2 + vkn(3)**2)
 IF (vkl == 0.0D0) GO TO 1010
 vk(1)  = vkn(1)/vkl
 vk(2)  = vkn(2)/vkl
 vk(3)  = vkn(3)/vkl
 pa     = .5D0*vkl
 
!     COMPUTE  SIDES -12- AND -41-
 
 DO  i = 1,3
   v12(i) = gp2(i) - gp1(i)
   v41(i) = gp1(i) - gp4(i)
 END DO
 
!     COMPUTE DOT PRODUCT, V12DK, OF V12 AND VK, THE VECTORS VP12,VI,VJ
 
 v12dk   = v12(1)*vk(1) + v12(2)*vk(2) + v12(3)*vk(3)
 vp12(1) = v12(1) - v12dk*vk(1)
 vp12(2) = v12(2) - v12dk*vk(2)
 vp12(3) = v12(3) - v12dk*vk(3)
 vp12l   = DSQRT(vp12(1)**2 + vp12(2)**2 + vp12(3)**2)
 IF (vp12l == 0.0D0) GO TO 1020
 vi(1)   = vp12(1)/vp12l
 vi(2)   = vp12(2)/vp12l
 vi(3)   = vp12(3)/vp12l
 vj(1)   = vk(2)*vi(3) - vk(3)*vi(2)
 vj(2)   = vk(3)*vi(1) - vk(1)*vi(3)
 vj(3)   = vk(1)*vi(2) - vk(2)*vi(1)
 
!     NORMALIZE J FOR GOOD MEASURE
 
 vjl   = DSQRT(vj(1)**2 + vj(2)**2 + vj(3)**2)
 IF (vjl == 0.0D0) GO TO 1030
 vj(1) = vj(1)/vjl
 vj(2) = vj(2)/vjl
 vj(3) = vj(3)/vjl
 x1    = 0.0D0
 y1    = 0.0D0
 x2    = vp12l
 y2    = 0.0D0
 x3    = vi(1)*vd1(1) + vi(2)*vd1(2) + vi(3)*vd1(3)
 y3    = vj(1)*vd1(1) + vj(2)*vd1(2) + vj(3)*vd1(3)
 x4    =-vi(1)*v41(1) - vi(2)*v41(2) - vi(3)*v41(3)
 y4    =-vj(1)*v41(1) - vj(2)*v41(2) - vj(3)*v41(3)
 
!     CHECK TO SEE IF INTERIOR ANGLES ARE LESS THAN 180 DEGREES.
!     IF NOT, CALL FATAL ERROR MESSAGE.
 
 IF (y3 <=    0.0D0) GO TO 1040
 IF (x3 <= y3*x4/y4) GO TO 1050
 IF (y4 <=    0.0D0) GO TO 1060
 IF (x4 >= x2-(x2-x3)*y4/y3) GO TO 1070
 
!     TEST FOR PARALLEL EFFECTS.
 
 cep1 = DABS((y3-y4)/(x3-x4))
 temp = x3 - x2
 cep2 = DABS((y4*temp-y3*x4)/(x4*temp+y4*y3))
 ep   = 1.0D-1
 IF (cep1 < ep) GO TO 15
 IF (cep2 < ep) GO TO 30
 GO TO 50
 15 IF (cep2 < ep) GO TO 40
 
!     AT THIS POINT THE LINE CONNECTING POINTS 3 AND 4 IS -PARALLEL- TO
!     THE LINE CONNECTING POINTS 1 AND 2.
 
 temp = y3*x4 - y4*(x3-x2)
 yp   = x2*y3*y4/temp
 p(1) = yp - y1
 p(2) = yp - y2
 p(3) = yp - y3
 p(4) = yp - y4
 xp   = x2*y3*x4/temp
 sa   = (x2 - xp)/yp
 c    = (x1 - xp)/yp
 z    = ((p(1)*p(2)*pa)/(p(3)*p(4)*2.0D0*g*t))*  &
     (1.0D0 + c23*nuc*(sa**2 + sa*c + c**2))
 GO TO 60
 
!     AT THIS POINT THE LINE CONNECTING POINTS 1 AND 4 IS -PARALLEL- TO
!     THE LINE CONNECTING POINTS 2 AND 3.
 
 30 d    = -.5D0*(x4/y4 + (x3-x2)/y3)
 xq   = x4 - y4*(x3-x4)/(y3-y4)
 temp = 1.0D0/DSQRT(1.0D0 + d**2)
 p(1) = (xq - x1 - d*y1)*temp
 p(2) = (xq - x2 - d*y2)*temp
 p(3) = (xq - x3 - d*y3)*temp
 p(4) = (xq - x4 - d*y4)*temp
 temp =  xq - x4
 b    = (temp*d + y4)/(temp - y4*d)
 z    = ((p(1)*p(2)*pa)/(p(3)*p(4)*2.0D0*g*t))*  &
     (1.0D0 + c23*nuc*(b**2 + b*d + d**2))
 GO TO 60
 
!     IN THIS CASE THE PANEL APPROXIMATES A PARALLELOGRAM.
 
 40 DO  i = 1,4
   p(i) = 1.0D0
 END DO
 d = -.5D0*(x4/y4 + (x3-x2)/y3 + (y3-y4)/(x3-x4))
 z = pa/(2.0D0*g*t)*(1.0D0 + 2.0D0*d**2*nuc)
 GO TO 60
 
!     IN THIS CASE NO PARALLEL EFFECTS EXIST.
 
 50 xq    = x4 - (x3-x4)/(y3-y4)*y4
 temp  = y3*x4 - y4*(x3-x2)
 xp    = x2*y3*x4/temp
 yp    = x2*y3*y4/temp
 xl    = DSQRT((xq-xp)**2 + yp**2)
 d     = (xq-xp)/yp
 temp  = yp/xl
 p(1)  = temp*(xq - x1 - d*y1)
 p(2)  = temp*(xq - x2 - d*y2)
 p(3)  = temp*(xq - x3 - d*y3)
 p(4)  = temp*(xq - x4 - d*y4)
 c     = xl/p(1) - d
 b     = xl/p(4) - c
 a     = xl/p(2) - d
 a2    = a**2
 b2    = b**2
 c2    = c**2
 d2    = d**2
 a3    = a2*a
 b3    = b2*b
 c3    = c2*c
 d3    = d2*d
 a4    = a3*a
 b4    = b3*b
 c4    = c3*c
 d4    = d3*d
 a5    = a4*a
 b5    = b4*b
 c5    = c4*c
 d5    = d4*d
 temp  = .5D0*p(1)*p(2)*p(3)*p(4)/xl**2
 term  = a + b + c23*(a3+b3) + .2D0*(a5+b5)
 term1 = c + d + c23*(c3+d3) + .2D0*(c5+d5)
 term2 = b + c + c23*(b3+c3) + .2D0*(b5+c5)
 term3 = d + a + c23*(d3+a3) + .2D0*(d5+a5)
 term  = term *DLOG(DABS(a+b))
 term1 = term1*DLOG(DABS(c+d))
 term2 = term2*DLOG(DABS(b+c))
 term3 = term3*DLOG(DABS(d+a))
 term4 = .1D0*((a2-c2)*(b3-d3) + (b2-d2)*(a3-c3))
 term5 = .2D0*((a -c )*(b4-d4) + (b -d )*(a4-c4))
 f     = temp*(term + term1 - term2 - term3 + term4 - term5)
 z     = p(1)*p(2)/(p(3)*p(4)*2.0D0*g*t)*(pa+4.0D0*nuc*(f-c23*pa))
 60 xl13  = DSQRT(x3**2 + y3**2)
 xl24  = DSQRT((x4-x2)**2 + y4**2)
 smallu(1) = x3/xl13
 smallu(2) = (x4-x2)/xl24
 smallu(3) = smallu(1)
 smallu(4) = smallu(2)
 smallv(1) = y3/xl13
 smallv(2) = y4/xl24
 smallv(3) = smallv(1)
 smallv(4) = smallv(2)
 temp      = x4*y3 - x3*y4
 avec(1)   =-.5D0*x2*y4*xl13/temp
 avec(2)   = .5D0*x2*y3*xl24/(temp-x2*(y3-y4))
 avec(3)   =-avec(1)
 avec(4)   =-avec(2)
 
!     COMPUTE THE SUM GIVEN ON P. 16 OF FMMS-39
 
 sum = 0.0D0
 DO  i = 1,4
   ivlbeg   = 1
   vleft(1) = smallu(i)*vi(1) + smallv(i)*vj(1)
   vleft(2) = smallu(i)*vi(2) + smallv(i)*vj(2)
   vleft(3) = smallu(i)*vi(3) + smallv(i)*vj(3)
   IF (iecpt(4*i+5) == 0) GO TO 70
   CALL transd (iecpt(4*i+5),ti)
   ivlbeg = 4
   CALL gmmatd (vleft(1),3,1,1, ti,3,3,0, vleft(4))
   70 k = 24 + 3*i
   ui(1) = ecpt(k+1)
   ui(2) = ecpt(k+2)
   ui(3) = ecpt(k+3)
   CALL gmmatd (vleft(ivlbeg),3,1,1, ui,3,1,0, dpterm)
   sum = sum + avec(i)*dpterm
 END DO
 f13 =-avec(1)*sum/(2.0D0*z)
 f24 = avec(2)*f13/avec(1)
 
!     SEARCH LIST OF SIL NOS. IN THE ECPT FOR THE PIVOT POINT.
 
 DO  i = 1,4
   ii = i
   IF (npvt == iecpt(i+1)) GO TO 100
 END DO
 CALL mesage (-30,34,iecpt(1))
 100 IF (ii == 2 .OR. ii == 4) GO TO 110
 fxx = f13/xl13
 i   = 1
 GO TO 120
 110 fxx = f24/xl24
 i   = 2
 120 jj(1) = -vi(1)*smallv(i) + vj(1)*smallu(i)
 jj(2) = -vi(2)*smallv(i) + vj(2)*smallu(i)
 jj(3) = -vi(3)*smallv(i) + vj(3)*smallu(i)
 
!                     T            T
!     COMPUTE  JJ X JJ  AND VK X VK
 
 CALL gmmatd (jj,3,1,0, jj,3,1,1, j3x3)
 CALL gmmatd (vk,3,1,0, vk,3,1,1, k3x3)
 
!     SUM THE TWO IN J3X3
 
 DO  j = 1,9
   j3x3(j) = j3x3(j) + k3x3(j)
 END DO
 SELECT CASE ( ii )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 160
   CASE (    4)
     GO TO 170
 END SELECT
 140 kk = 3
 GO TO 180
 150 kk = 4
 GO TO 180
 160 kk = 1
 GO TO 180
 170 kk = 2
 
!     ZERO OUT KE
 
 180 DO  i = 1,36
   ke(i) = 0.0D0
 END DO
 
!                 D
!     SET UP THE K   MATRIX
!                 II
 
 mpoint = 1
 IF (iecpt(4*ii+5) == 0) GO TO 200
 CALL transd (ecpt(4*ii+5),ti)
 mpoint = 10
 CALL gmmatd (ti,3,3,1, j3x3(1),3,3,0, k3x3(1))
 CALL gmmatd (k3x3(1),3,3,0, ti,3,3,0, j3x3(1))
 200 k = 1
 j = ii
 210 ke( 1) = fxx*j3x3(k  )
 ke( 2) = fxx*j3x3(k+1)
 ke( 3) = fxx*j3x3(k+2)
 ke( 7) = fxx*j3x3(k+3)
 ke( 8) = fxx*j3x3(k+4)
 ke( 9) = fxx*j3x3(k+5)
 ke(13) = fxx*j3x3(k+6)
 ke(14) = fxx*j3x3(k+7)
 ke(15) = fxx*j3x3(k+8)
 CALL ds1b (ke,iecpt(j+1))
 IF (j == kk) RETURN
 
!                 D
!     SET UP THE K   MATRIX
!                 IJ
 
 j = kk
 IF (iecpt(4*j+5) == 0) GO TO 220
 CALL transd (ecpt(4*j+5),ti)
 npoint = 10
 IF (mpoint == 10) npoint = 1
 CALL gmmatd (j3x3(mpoint),3,3,0, ti,3,3,0, j3x3(npoint))
 k = npoint
 GO TO 230
 220 k = mpoint
 230 fxx = -fxx
 GO TO 210
 
!     ERROR RETURNS
 
 1010 CONTINUE
 1020 CONTINUE
 1030 CALL mesage (30,26,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
 1040 iecpt(2) = 2
 GO TO 2000
 1050 iecpt(2) = 4
 GO TO 2000
 1060 iecpt(2) = 1
 GO TO 2000
 1070 iecpt(2) = 3
 2000 CALL mesage (30,27,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
END SUBROUTINE dshear
