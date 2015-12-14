SUBROUTINE spanl1(iarg)
!*****
! THIS ROUTINE COMPUTES PHASE I PARAMETERS FOR STRESS DATA RECOVERY FOR
! THE SHEAR PANEL (IF IARG = 4) AND THE TWIST PANEL (IF IARG = 5).
! MUCH OF THE CODE WAS LIFTED FROM SUBROUTIVE KPANEL
!*****
 
!                 E C P T  F O R  B O T H  P A N E L S
! ECPT( 1)  -  IELID          ELEMENT ID. NO.
! ECPT( 2)  -  ISILNO(4)      SCALAR INDEX NUMBERS
! ECPT( 3)  -   ...                   ...
! ECPT( 4)  -   ...                   ...
! ECPT( 5)  -   ...                   ...
! ECPT( 6)  -  MATID          MATERIAL ID.
! ECPT( 7)  -  T              THICKNESS
! ECPT( 8)  -  FMU            NON-STRUCTURAL MASS
! ECPT( 9)  -  ICSID1         COOR. SYS. ID. FOR GRID POINT 1
! ECPT(10)  -  GP1(3)         BASIC COORDINATES FOR GRID POINT 1
! ECPT(11)  -   ...                      ...
! ECPT(12)  -   ...                      ...
! ECPT(13)  -  ICSID2         COOR. SYS. ID. FOR GRID POINT 2
! ECPT(14)  -  GP2(3)         BASIC COORDINATES FOR GRID POINT 2
! ECPT(15)  -   ...                      ...
! ECPT(16)  -   ...                      ...
! ECPT(17)  -  ICSID3         COOR. SYS. ID. FOR GRID POINT 3
! ECPT(18)  -  GP3(3)         BASIC COORDINATES FOR GRID POINT 3
! ECPT(19)  -   ...                      ...
! ECPT(20)  -   ...                      ...
! ECPT(21)  -  ICSID4         COOR. SYS. ID. FOR GRID POINT 4
! ECPT(22)  -  GP4(3)         BASIC COORDINATES FOR GRID POINT 4
! ECPT(23)  -   ...                      ...
! ECPT(24)  -   ...                      ...
! ECPT(25)  -  TEMPEL         ELEMENT TEMPERATURE
 
 
 
 
 INTEGER, INTENT(IN)                      :: iarg
 REAL :: nu
 
 
 
 
 
 
 DIMENSION vd1(3)             ,vd2(3)  &
     ,                  vkn(3)             ,vk(3)  &
     ,                  v12(3)             ,v41(3)  &
     ,                  vp12(3)            ,vi(3)  &
     ,                  vj(3)              ,avec(4)  &
     ,                  smallu(4)          ,smallv(4)  &
     ,                  p(4)               ,iecpt(100)  &
     ,                  ecpt(100) ,                  vleft(6)           ,ti(9)
 
! SDR2 PHASE I INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x5/ ielid              ,isilno(4)  &
     ,                  matid              ,t  &
     ,                  fmu                ,icsid1  &
     ,                  gp1(3)             ,icsid2  &
     ,                  gp2(3)             ,icsid3  &
     ,                  gp3(3)             ,icsid4  &
     ,                  gp4(3)             ,tempel ,                  xxxxxx(75)
 COMMON   /sdr2x5/ jelid              ,jsilno(4)  &
     ,                  s(3,4)             ,out(15) ,                  yyyyyy(93)
 
! SDR2 SCRATCH BLOCK
 
 COMMON   /sdr2x6/ vleft              ,ti  &
     ,                  spcon ,                  vd1                ,vd2  &
     ,                  vkn                ,vk  &
     ,                  v12                ,v41  &
     ,                  vp12               ,vi  &
     ,                  vj                 ,avec  &
     ,                  smallu             ,smallv  &
     ,                  p                  ,x1  &
     ,                  x2                 ,x3  &
     ,                  x4                 ,y1  &
     ,                  y2                 ,y3  &
     ,                  y4                 ,vkl  &
     ,                  pa                 ,v12dk  &
     ,                  cep1               ,cep2  &
     ,                  ep                 ,temp
 COMMON   /sdr2x6/ yp                 ,xp  &
     ,                  sa                 ,xq  &
     ,                  b                  ,xl  &
     ,                  a                  ,a2  &
     ,                  a3                 ,a4  &
     ,                  a5                 ,b2  &
     ,                  b3                 ,b4  &
     ,                  b5                 ,c  &
     ,                  c2                 ,c3  &
     ,                  c4                 ,c5  &
     ,                  d                  ,d2  &
     ,                  d3                 ,d4  &
     ,                  d5                 ,term1  &
     ,                  term2              ,term3  &
     ,                  term4              ,term5  &
     ,                  xl13               ,xl24
 
! INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON   /matin/ matidc             ,matflg  &
     ,                  eltemp             ,stress  &
     ,                  sinth              ,costh
 
 
 
 COMMON   /matout/ e                  ,g  &
     ,                  nu                 ,rho  &
     ,                  alpha              ,tsubo  &
     ,                  gsube              ,sigt  &
     ,                  sigc               ,sigs
 
 
 
 EQUIVALENCE (ielid,iecpt(1),ecpt(1))
 
! CALL MAT TO GET MATERIAL PROPERTIES.
 
 matidc = matid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 
! COMPUTE DIAGONAL VECTORS.
 
 DO  i=1,3
   vd1(i) = gp3(i) - gp1(i)
   vd2(i) = gp4(i) - gp2(i)
 END DO
 
! COMPUTE THE NORMAL VECTOR VKN, NORMALIZE, AND COMPUTE THE PROJECTED
! AREA, PA
 
 vkn(1) = vd1(2)*vd2(3)-vd1(3)*vd2(2)
 vkn(2) = vd1(3)*vd2(1)-vd1(1)*vd2(3)
 vkn(3) = vd1(1)*vd2(2)-vd1(2)*vd2(1)
 vkl = SQRT(vkn(1)**2+vkn(2)**2+vkn(3)**2)
 IF (vkl == 0.0) GO TO 160
 vk(1) = vkn(1)/vkl
 vk(2) = vkn(2)/vkl
 vk(3) = vkn(3)/vkl
 pa = .5 * vkl
 
! COMPUTE  SIDES -12- AND -41-
 
 DO  i=1,3
   v12(i) = gp2(i) - gp1(i)
   v41(i) = gp1(i) - gp4(i)
 END DO
 
! COMPUTE DOT PRODUCT, V12DK, OF V12 AND VK, THE VECTORS VP12, VI, VJ
 
 v12dk = v12(1)*vk(1)+v12(2)*vk(2)+v12(3)*vk(3)
 vp12(1) = v12(1)-v12dk*vk(1)
 vp12(2) = v12(2)-v12dk*vk(2)
 vp12(3) = v12(3)-v12dk*vk(3)
 vp12l = SQRT(vp12(1)**2+vp12(2)**2+vp12(3)**2)
 IF (vp12l == 0.0) GO TO 170
 vi(1) = vp12(1)/vp12l
 vi(2) = vp12(2)/vp12l
 vi(3) = vp12(3)/vp12l
 vj(1) = vk(2)*vi(3)-vk(3)*vi(2)
 vj(2) = vk(3)*vi(1)-vk(1)*vi(3)
 vj(3) = vk(1)*vi(2)-vk(2)*vi(1)
 
! NORMALIZE J FOR GOOD MEASURE
 
 vjl = SQRT (vj(1)**2  +  vj(2)**2  +  vj(3)**2)
 IF (vjl == 0.0) GO TO 180
 vj(1) = vj(1) / vjl
 vj(2) = vj(2) / vjl
 vj(3) = vj(3) / vjl
 x1 = 0.0
 y1 = 0.0
 x2 = vp12l
 y2 = 0.0
 x3 = vi(1) * vd1(1)  +  vi(2) * vd1(2)  +  vi(3) * vd1(3)
 y3 = vj(1) * vd1(1)  +  vj(2) * vd1(2)  +  vj(3) * vd1(3)
 x4 =-vi(1) * v41(1)  -  vi(2) * v41(2)  -  vi(3) * v41(3)
 y4 =-vj(1) * v41(1)  -  vj(2) * v41(2)  -  vj(3) * v41(3)
 
! CHECK TO SEE IF INTERIOR ANGLES ARE LESS THAN 180 DEGREES.  IF NOT,
! CALL FATAL ERROR MESSAGE.
 
 IF (y3 <= 0.0) GO TO 190
 IF (x3 <= y3*x4/y4) GO TO 200
 IF (y4 <= 0.0) GO TO 210
 IF (x4 >= x2 - (x2-x3)*y4/y3) GO TO 220
 
! TEST FOR PARALLEL EFFECTS.
 
 temp = x3 - x2
 ep = 0.01
 IF (ABS(y3-y4) < ABS(x3-x4)*ep) GO TO 30
 IF (ABS(y4*temp-y3*x4) < ABS(x4*temp+y4*y3)*ep) GO TO 40
 GO TO 70
 30 IF (ABS(y4*temp-y3*x4) < ABS(x4*temp+y4*y3)*ep) GO TO 50
 
! AT THIS POINT THE LINE CONNECTING POINTS 3 AND 4 IS -PARALLEL- TO THE
! LINE CONNECTING POINTS 1 AND 2.
 
 temp = y3*x4  -  y4 * (x3-x2)
 yp   = x2*y3*y4 / temp
 p(1) = yp - y1
 p(2) = yp - y2
 p(3) = yp - y3
 p(4) = yp - y4
 xp   = x2*y3*x4 / temp
 sa   =(x2 - xp) / yp
 c    =(x1 - xp) / yp
 z    =(  (p(1)*p(2)*pa) / (p(3)*p(4)*2.0*g*t)  ) *  &
     (  1.0  +  2.0/(3.0 + 3.0*nu) * (sa**2 + sa*c + c**2)  )
 GO TO 80
 
! AT THIS POINT THE LINE CONNECTING POINTS 1 AND 4 IS -PARALLEL- TO THE
! LINE CONNECTING POINTS 2 AND 3.
 
 40 d    = -.5 * (  x4/y4  +  (x3-x2)/y3  )
 xq   = x4  - y4  *  (x3-x4)/(y3-y4)
 temp = 1.0 / SQRT (1.0 + d**2)
 p(1) = ( xq - x1 - d*y1) * temp
 p(2) = ( xq - x2 - d*y2) * temp
 p(3) = ( xq - x3 - d*y3) * temp
 p(4) = ( xq - x4 - d*y4) * temp
 temp =   xq - x4
 b    =   (temp * d  +  y4)  /  (temp  -  y4*d)
 z    =(  (p(1)*p(2)*pa) / (p(3)*p(4)*2.0*g*t)  ) *  &
     (  1.0  +  2.0/(3.0 + 3.0*nu) * (b**2 + b*d + d**2)  )
 GO TO 80
 
! IN THIS CASE THE PANEL APPROXIMATES A PARALLELOGRAM.
 
 50 DO  i=1,4
   p(i) = 1.0
 END DO
 d    = -.5 * (  x4/y4  +  (x3-x2)/y3  +  (y3-y4)/(x3-x4)  )
 z    = pa / (2.0*g*t) * (1.0 + 2.0*d**2/(1.0+nu))
 GO TO 80
 
! IN THIS CASE NO PARALLEL EFFECTS EXIST.
 
 70 xq   = x4  -  (x3-x4)/(y3-y4) * y4
 temp = y3*x4  -  y4*(x3-x2)
 xp   = x2*y3*x4 / temp
 yp   = x2*y3*y4 / temp
 xl   = SQRT ( (xq-xp)**2 + yp**2 )
 d    = (xq-xp)/yp
 temp = yp/xl
 p(1) = temp * (xq - x1 - d*y1)
 p(2) = temp * (xq - x2 - d*y2)
 p(3) = temp * (xq - x3 - d*y3)
 p(4) = temp * (xq - x4 - d*y4)
 c    = xl/p(1) - d
 b    = xl/p(4) - c
 a    = xl/p(2) - d
 a2   = a**2
 b2   = b**2
 c2   = c**2
 d2   = d**2
 a3   = a2*a
 b3   = b2*b
 c3   = c2*c
 d3   = d2*d
 a4   = a3*a
 b4   = b3*b
 c4   = c3*c
 d4   = d3*d
 a5   = a4*a
 b5   = b4*b
 c5   = c4*c
 d5   = d4*d
 temp = .5 * p(1) * p(2) * p(3) * p(4) / xl**2
 term =    a  +  b  +  2.0*(a3+b3)/3.0  +  .2*(a5+b5)
 term1=    c  +  d  +  2.0*(c3+d3)/3.0  +  .2*(c5+d5)
 term2=    b  +  c  +  2.0*(b3+c3)/3.0  +  .2*(b5+c5)
 term3=    d  +  a  +  2.0*(d3+a3)/3.0  +  .2*(d5+a5)
 term =   term  * ALOG(ABS(a+b))
 term1=   term1 * ALOG(ABS(c+d))
 term2=   term2 * ALOG(ABS(b+c))
 term3=   term3 * ALOG(ABS(d+a))
 term4=  .1*( (a2-c2)*(b3-d3)  +  (b2-d2)*(a3-c3) )
 term5=  .2*( (a -c )*(b4-d4)  +  (b -d )*(a4-c4) )
 f    =  temp * (term + term1 - term2 - term3 + term4 - term5)
 z    =  p(1)*p(2) / (p(3)*p(4)*2.0*g*t) * (pa + 4.0/(1.0+nu) *  &
     (f - 2.0*pa/3.0))
 80 xl13 =  SQRT (x3**2 + y3**2)
 xl24 =  SQRT (  (x4-x2)**2  +  y4**2  )
 smallu(1) = x3/xl13
 smallu(2) = (x4-x2)/xl24
 smallu(3) = smallu(1)
 smallu(4) = smallu(2)
 smallv(1) = y3/xl13
 smallv(2) = y4/xl24
 smallv(3) = smallv(1)
 smallv(4) = smallv(2)
 temp = x4 * y3  -  x3 * y4
 avec(1) = -.5 * x2 * y4 * xl13 / temp
 avec(2) = .5 * x2 * y3 * xl24 / (temp - x2 * (y3-y4) )
 avec(3) = - avec(1)
 avec(4) = - avec(2)
 
! IF IARG = 4, WE HAVE A SHEAR PANEL, AND IF IARG = 5, A TWIST PANEL.
 
 IF (iarg == 4) GO TO 100
 
! SINCE WE ARE DEALING WITH A TWIST PANEL STORE -SMALLV IN SMALLU AND
! SMALLU IN SMALLV.
 
 DO  i=1,4
   temp = smallu(i)
   smallu(i) = -smallv(i)
   smallv(i) = temp
 END DO
 
! COMPUTE THE SINGLE PRECISION CONSTANT SPCON
 
 100 IF (iarg == 5) GO TO 110
 spcon = -1.0/ (2.0 * z * t)
 GO TO 120
 110 spcon = -1.0/ (4.0 * z)
 
! COMPUTE THE FOUR 1 X 3 MATRICES S
 
 120 DO  i=1,4
   ivlbeg = 1
   vleft(1) = smallu(i) * vi(1)  +  smallv(i) * vj(1)
   vleft(2) = smallu(i) * vi(2)  +  smallv(i) * vj(2)
   vleft(3) = smallu(i) * vi(3)  +  smallv(i) * vj(3)
   IF (iecpt(4*i+5) == 0) GO TO 130
   ivlbeg = 4
   CALL transs (iecpt(4*i+5),ti)
   CALL gmmats (vleft(1),3,1,1, ti,3,3,0, vleft(4) )
   130 CONTINUE
   s(1,i) = spcon * vleft(ivlbeg  ) * avec(i)
   s(2,i) = spcon * vleft(ivlbeg+1) * avec(i)
   s(3,i) = spcon * vleft(ivlbeg+2) * avec(i)
 END DO
 out(1) = avec(1)
 out(2) = avec(2)
 out(3) = t
 out(4) = p(2) / p(1)
 out(5) = p(1) * p(2) / p(3)**2
 out(6) = p(1) * p(2) / p(4)**2
 out(7) = sigs
 jelid = ielid
 DO  i=1,4
   jsilno(i) = isilno(i)
 END DO
 IF( iarg /= 4 ) RETURN
!*****
!  ADDITIONAL PHASE-1 OUTPUTS FOR SHEAR PANEL FORCES  IN PHASE 2
!*****
 out(8) = p(1) / p(3) *t
 out(9) = ( p(1)*p(2) ) / ( p(3)*p(4) ) * t
 out(10) = p(2) / p(4) * t
 out(11)= -v12dk / 2.0
 out(12)= x2 / 2.0
 out(13)= SQRT( (x3-x2)**2 + y3**2 ) / 2.0
 out(14)= SQRT( (x4-x3)**2 + (y4-y3)**2 ) / 2.0
 out(15)= SQRT( x4**2 + y4**2 ) / 2.0
 RETURN
 160 CONTINUE
 170 CONTINUE
 180 CALL mesage (-30,26,iecpt(1))
 190 iecpt(2) = 2
 GO TO 230
 200 iecpt(2) = 4
 GO TO 230
 210 iecpt(2) = 1
 GO TO 230
 220 iecpt(2) = 3
 230 CALL mesage (-30,27,iecpt(1))
 RETURN
END SUBROUTINE spanl1
