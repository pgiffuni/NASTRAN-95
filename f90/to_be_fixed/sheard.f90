SUBROUTINE sheard
     
!     THIS SUBROUTINE COMPUTES THE 12 X 12 STIFFNESS MATRIX FOR THE
!     SHEAR PANEL ELEMENT, AS WELL AS ITS DIAGONALIZED MASS MATRIX.
 
!     DOUBLE PRECISION VERSION
 
!     ECPT FOR THE SHEAR PANEL ELEMENT
 
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
 
 LOGICAL :: iheat,nogo
 INTEGER :: ipart(4),dict(11),estid
 REAL :: nu
 DOUBLE PRECISION :: cepx,cepy,ep,me(144),kout(144),mout(144),ke(144),  &
     t,nuc,g,e,c23,vleft(6),vright(6),ti(9),p(4),x1,  &
     y1,x2,y2,x3,y3,x4,y4,cep1,cep2,temp,yp,xp,sa,a,b,  &
     c,d,term,term1,term2,term3,term4,term5,f,xl13,  &
     xl24,con,z,xl,vd1(3),vd2(3),vkn(3),vk(3),v12(3),  &
     v41(3),vp12(3),vi(3),vj(3),avec(4),smallu(4),  &
     smallv(4),a2,a3,a4,a5,b2,b3,b4,b5,c2,c3,c4,c5, d2,d3,d4,d5
 DIMENSION        iecpt(100),ecpt(100)
 COMMON /system/  ksystm(55),iheat
 COMMON /emgprm/  ixr,jcore,ncore,dum(12),ismb(3),iprec,nogo,heat
 COMMON /emgdic/  idm, ldict,ngrids,elid,estid
 COMMON /emgest/  ielid,isilno(4),matid,tsp,fmu,icsid1,gp1(3),  &
     icsid2,gp2(3),icsid3,gp3(3),icsid4,gp4(3),tempel
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matout/  esp,gsp,nu,rho,alpha,tsub0,gsube,sigt,sigc,sigs
 COMMON /matin /  matidc,matflg,eltemp,stress,sinth,costh
 EQUIVALENCE      (me(1),ke(1)),(kout(1),mout(1)),  &
     (iecpt(1),ecpt(1),ielid),(dict(5),dict5)
 DATA    ipart /  1,2,3,4 /
 
 ngrids = 4
 ldict  = 5 + ngrids
 
!     IF STIFFNESS MATRIX NOT NEEDED GO TO PERFORM MASS CALCULATIONS
 
 IF (ismb(1) == 0) GO TO 400
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 12
 dict(4) = 7
 ip      = iprec
 isort   = 0
 
!     CALL MAT TO GET MATERIAL PROPERTIES.
 
 matidc = matid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 dict5  = gsube
 
 t   = tsp
 g   = gsp
 e   = esp
 IF (t*g == 0.0) GO TO 7770
 c23 = 2.d0/3.d0
 nuc = 1.d0/(1.d0+nu)
 
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
 vkl    = DSQRT(vkn(1)**2 + vkn(2)**2  + vkn(3)**2)
 IF (vkl == 0.) GO TO 7770
 vk(1)  = vkn(1)/vkl
 vk(2)  = vkn(2)/vkl
 vk(3)  = vkn(3)/vkl
 pa     = vkl/2.
 
!     COMPUTE  SIDES -12- AND -41-
 
 DO  i = 1,3
   v12(i) = gp2(i) - gp1(i)
   v41(i) = gp1(i) - gp4(i)
 END DO
 
!     COMPUTE DOT PRODUCT, V12DK, OF V12 AND VK, THE VECTORS VP12, VI,
!     VJ
 
 v12dk   = v12(1)*vk(1) + v12(2)*vk(2) + v12(3)*vk(3)
 vp12(1) = v12(1) - v12dk*vk(1)
 vp12(2) = v12(2) - v12dk*vk(2)
 vp12(3) = v12(3) - v12dk*vk(3)
 vp12l   = DSQRT(vp12(1)**2 + vp12(2)**2 + vp12(3)**2)
 IF (vp12l == 0.) GO TO 7770
 vi(1) = vp12(1)/vp12l
 vi(2) = vp12(2)/vp12l
 vi(3) = vp12(3)/vp12l
 vj(1) = vk(2)*vi(3) - vk(3)*vi(2)
 vj(2) = vk(3)*vi(1) - vk(1)*vi(3)
 vj(3) = vk(1)*vi(2) - vk(2)*vi(1)
 
!     NORMALIZE J FOR GOOD MEASURE
 
 vjl = DSQRT(vj(1)**2 + vj(2)**2 + vj(3)**2)
 IF (vjl == 0.) GO TO 7770
 vj(1) = vj(1)/vjl
 vj(2) = vj(2)/vjl
 vj(3) = vj(3)/vjl
 x1 = 0.
 y1 = 0.
 x2 = vp12l
 y2 = 0.
 x3 = vi(1)*vd1(1) + vi(2)*vd1(2) + vi(3)*vd1(3)
 y3 = vj(1)*vd1(1) + vj(2)*vd1(2) + vj(3)*vd1(3)
 x4 =-vi(1)*v41(1) - vi(2)*v41(2) - vi(3)*v41(3)
 y4 =-vj(1)*v41(1) - vj(2)*v41(2) - vj(3)*v41(3)
 
!     CHECK TO SEE IF INTERIOR ANGLES ARE LESS THAN 180 DEGREES. IF NOT,
!     CALL FATAL ERROR MESSAGE.
 
 IF (y3 <= 0.) GO TO 7780
 IF (y4 <= 0.) GO TO 7800
 IF (x3 <= y3*x4/y4) GO TO 7790
 IF (x4 >= x2-(x2-x3)*y4/y3) GO TO 7810
 
!     TEST FOR PARALLEL EFFECTS.
 
 cep1 = DABS(y3-y4)
 cepx = DABS(x3-x4)
 temp = x3 - x2
 cep2 = DABS(y4*temp-y3*x4)
 cepy = DABS(x4*temp+y4*y3)
 ep   = 0.01D0
 IF (cep1 < ep*cepx) GO TO 30
 IF (cep2 < ep*cepy) GO TO 40
 GO TO 70
 30 IF (cep2 < ep*cepy) GO TO 50
 
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
 z    = ((p(1)*p(2)*pa)/(p(3)*p(4)*2.*g*t))*(1.+c23*nuc* (sa**2+sa*c+c**2))
 GO TO 80
 
!     AT THIS POINT THE LINE CONNECTING POINTS 1 AND 4 IS -PARALLEL- TO
!     THE LINE CONNECTING POINTS 2 AND 3.
 
 40 d    = -.5*(x4/y4 + (x3-x2)/y3)
 xq   = x4 - y4*(x3-x4)/(y3-y4)
 temp = 1.d0/DSQRT(1.d0+d**2)
 p(1) = (xq-x1-d*y1)*temp
 p(2) = (xq-x2-d*y2)*temp
 p(3) = (xq-x3-d*y3)*temp
 p(4) = (xq-x4-d*y4)*temp
 temp = xq - x4
 b    = (temp*d+y4)/(temp-y4*d)
 z    = ((p(1)*p(2)*pa)/(p(3)*p(4)*2.*g*t))*(1.+c23*nuc*(b**2+b*d + d**2))
 GO TO 80
 
!     IN THIS CASE THE PANEL APPROXIMATES A PARALLELOGRAM.
 
 50 DO  i = 1,4
   p(i) = 1.
 END DO
 d = -.5D0*(x4/y4+(x3-x2)/y3+(y3-y4)/(x3-x4))
 z = pa/(2.*g*t)*(1.+2.*d**2*nuc)
 GO TO 80
 
!     IN THIS CASE NO PARALLEL EFFECTS EXIST.
 
 70 xq   = x4 - (x3-x4)/(y3-y4)*y4
 temp = y3*x4 - y4*(x3-x2)
 xp   = x2*y3*x4/temp
 yp   = x2*y3*y4/temp
 xl   = DSQRT((xq-xp)**2 + yp**2)
 d    = (xq-xp)/yp
 temp = yp/xl
 p(1) = temp*(xq-x1-d*y1)
 p(2) = temp*(xq-x2-d*y2)
 p(3) = temp*(xq-x3-d*y3)
 p(4) = temp*(xq-x4-d*y4)
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
 temp = .5D0*p(1)*p(2)*p(3)*p(4)/xl**2
 term = (a + b + c23*(a3+b3) + .2D0*(a5+b5))*DLOG(DABS(a+b))
 term1= (c + d + c23*(c3+d3) + .2D0*(c5+d5))*DLOG(DABS(c+d))
 term2= (b + c + c23*(b3+c3) + .2D0*(b5+c5))*DLOG(DABS(b+c))
 term3= (d + a + c23*(d3+a3) + .2D0*(d5+a5))*DLOG(DABS(d+a))
 term4= .1D0*((a2-c2)*(b3-d3)+ (b2-d2)*(a3-c3))
 term5= .2D0*((a-c)*(b4-d4) + (b-d)*(a4-c4))
 f    = temp*(term+term1-term2-term3+term4-term5)
 z    = p(1)*p(2)/(p(3)*p(4)*2.*g*t)*(pa+4.*nuc*(f-c23*pa))
 80 xl13 = DSQRT(x3**2 + y3**2)
 xl24 = DSQRT((x4-x2)**2 + y4**2)
 smallu(1) = x3/xl13
 smallu(2) = (x4-x2)/xl24
 smallu(3) = smallu(1)
 smallu(4) = smallu(2)
 smallv(1) = y3/xl13
 smallv(2) = y4/xl24
 smallv(3) = smallv(1)
 smallv(4) = smallv(2)
 temp    = x4*y3 - x3*y4
 avec(1) = -.5*x2*y4*xl13/temp
 avec(2) =  .5*x2*y3 *xl24/(temp -x2*(y3-y4))
 avec(3) = -avec(1)
 avec(4) = -avec(2)
 
 DO  i = 1,144
   ke(i) = 0.
 END DO
 DO  ipvt = 1,4
   con = avec(ipvt)/(2.*z)
   
!     COMPUTE THE -VLEFT- VECTOR
   
   ivlbeg = 1
   vleft(1) = vi(1)*smallu(ipvt) + vj(1)*smallv(ipvt)
   vleft(2) = vi(2)*smallu(ipvt) + vj(2)*smallv(ipvt)
   vleft(3) = vi(3)*smallu(ipvt) + vj(3)*smallv(ipvt)
   IF (iecpt(4*ipvt+5) == 0) GO TO 150
   CALL transd (iecpt(4*ipvt+5),ti)
   ivlbeg = 4
   CALL gmmatd (ti,3,3,1, vleft(1),3,1,0, vleft(4))
   
!     COMPUTE THE 6 X 6 -S
   
   150 DO  j = 1,4
     ivrbeg = 1
     vright(1) = smallu(j)*vi(1) + smallv(j)*vj(1)
     vright(2) = smallu(j)*vi(2) + smallv(j)*vj(2)
     vright(3) = smallu(j)*vi(3) + smallv(j)*vj(3)
     IF (iecpt(4*j+5) == 0) GO TO 170
     CALL transd (iecpt(4*j+5),ti)
     CALL gmmatd (vright(1),1,3,0, ti,3,3,0, vright(4))
     ivrbeg = 4
     170   jt = (ipvt-1)*36 + (j-1)*9 + 1
     CALL gmmatd (vleft(ivlbeg),3,1,0, vright(ivrbeg),1,3,0, ke(jt))
     jt8 = jt + 8
     DO  k = jt,jt8
       ke(k) = con*ke(k)*avec(j)
     END DO
   END DO
 END DO
 
!     NOW REARRANGE KE BY INCREASING SIL THEN OUTPUT IT VIA EMGOUT
!     FIRST DETERMINE WHAT INCREASING SIL ORDER WILL BE
 
 ASSIGN 290 TO k OR m
 275 CONTINUE
 DO  i = 1,3
   ip1 = i + 1
   it  = ipart(i)
   DO  j = ip1,4
     jt = ipart(j)
     IF (isilno(it) <= isilno(jt)) CYCLE
     ipart(i) = jt
     ipart(j) = it
     it = jt
     GO TO 275
   END DO
 END DO
 isort = 1
 GO TO korm, (290,420)
 
!     NOW REARRANGE TERMS IN THE STIFFNESS MATRIX KE AND STORE IN KOUT
 
!     KE = (K  ,K  ,K  ,K  ,K  ,...,K  ,K  ,...,K  )
!            11  12  13  14  21      24  31      44
 
!     WHERE  K  IS A 3X3 SUBMATRIX AND  SILS ARE IN GRID POINT ORDER
!             IJ
 
!     AND    *****                  ****
!            * K     K     K     K     *
!            *  L1L1  L1L2  L1L3  L1L4 *
!            *                         *
!            * K     K     K     K     *
!     KOUT = *  L2L1  L2L2  L2L3  L2L4 *
!            *                         *
!            * K     K     K     K     *
!            *  L3L1  L3L2  L3L3  L3L4 *
!            *                         *
!            * K     K     K     K     *
!            *  L4L1  L4L2  L4L3  L4L4 *
!            ****                   ****
 
!     WHERE  KOUT     IS A   3X3    MATRIX AND SILS ARE IN INCREASING
!                LILJ
!     ORDER
 
 290  CONTINUE
 DO  i = 1,4
   is = ipart(i)
   DO  j = 1,4
     js = ipart(j)
     DO  k = 1,3
       DO  l = 1,3
         iout = (i -1)*36 + (j -1)*3 + (k-1)*12 + l
         ike  = (is-1)*36 + (js-1)*9 + (k-1)* 3 + l
         kout(iout) = ke(ike)
       END DO
     END DO
   END DO
 END DO
 
!     OUTPUT THE STIFFNESS MATRIX
 
 CALL emgout (kout,kout,144,1,dict,1,ip)
 
!     HERE WE CALCULATE THE MASS MATRIX VIA SUBROUTINE EMASTQ
 
 
 400 IF (ismb(2) == 0) RETURN
 
 CALL emadtq (6,me)
 IF (isort == 1) GO TO 420
 ASSIGN 420 TO korm
 GO TO 275
 
!     RETURN WITH A GRID POINT SORT ARRAY IN IPART
 
 420 DO  i = 1,4
   it = 1 + (ipart(i)-1)*3
   ij = (i-1)*3 + 1
   mout(ij  ) = me(it  )
   mout(ij+1) = me(it+1)
   mout(ij+2) = me(it+2)
 END DO
 
 dict(1) = estid
 dict(2) = 2
 dict(3) = 12
 dict(4) =  7
 dict5   = 0.
 
 CALL  emgout (kout,kout,12,1,dict,2,ip)
 RETURN
 
!     ERROR EXITS
 
 7770 CALL mesage (30,26,iecpt(1))
 7777 nogo = .true.
 RETURN
 
 7780 iecpt(2) = 2
 GO TO 7820
 7790 iecpt(2) = 4
 GO TO 7820
 7800 iecpt(2) = 1
 GO TO 7820
 7810 iecpt(2) = 3
 7820 CALL mesage (30,27,iecpt(1))
 GO TO 7777
END SUBROUTINE sheard
