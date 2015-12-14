SUBROUTINE dis2d8
     
!     2-D, 8 GRID POINT ISOPARAMETRIC STRUCTURAL ELEMENT
!     DIFFERENTIAL STIFFNESS MATRIX ROUTINE
 
 REAL :: kx,ky
 DOUBLE PRECISION :: kij,g,b,xi,eta,dnxi,dneta,xx,tb,dnl,pt,h,xjb,  &
     xxjb,determ,dnc,gsube,dumarg, bt,tempar,temp,dnx,dny,SAVE,tsave,  &
     kwd(36),cid(18),cjd(18),kmult(18),dhh, thick,pstmul(9),premul(9),e1t
 DIMENSION        sig(3),alphas(3),st(3),semp(9),ttb(9),stb(9),  &
     bb(72),db(72),s(6),r(9),se1t(6),dn(8),  &
     g(9),qq(15),xi(8),eta(8),tb(9),xy1(3),xy2(3),  &
     b(12),bt(12),ecpt(1),temp(9),tempar(1),dnx(1),  &
     dny(1),dnxi(1),dneta(1),SAVE(72),tsave(72),  &
     vec(3),vvec(3),veci(3),vecj(3),veck(3),e1t(9), iws(2,3)
 COMMON /ds1aaa/  npvt,icstm,ncstm
 COMMON /ds1aet/  necpt(1),ngrid(8),id1,th,matid1,t,isys1,x1,y1,z1,  &
     isys2,x2,y2,z2,isys3,x3,y3,z3,isys4,x4,y4,z4,  &
     isys5,x5,y5,z5,isys6,x6,y6,z6,isys7,x7,y7,z7,  &
     isys8,x8,y8,z8,ttemp,edt,isetno,tgrid(8),disp(24)
 COMMON /matin /  matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/  g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     tref,GE,kx,ky,c
 COMMON /ds1adp/  kij(36),xx(16),dnc(16),dnl(16),xxjb(2,2),xjb(4),  &
     pt(3),h(3),g,b,bt,tb,determ,gsube,dumarg,tsave
 EQUIVALENCE      (alphas(1),alpha1), (ecpt(1),necpt(1)),(temp(1),b(1)),  &
     (dnc(1),dnxi(1)),(dnc(9),dneta(1)),  &
     (dnl(1),dnx(1)),(dnl(9),dny(1)),(qq(1),g11),  &
     (tempar(1),bt(1)),(xy1(1),x1),(xy2(1),x2)
 DATA     xi   /  -1.d0, 1.d0, 1.d0,-1.d0, 0.d0, 1.d0, 0.d0,-1.d0/
 DATA     eta  /  -1.d0,-1.d0, 1.d0, 1.d0,-1.d0, 0.d0, 1.d0, 0.d0/
 
!     ECPT LIST
!                                      IN
!                                      THIS
!     ECPT       DESCRIPTION           ROUTINE        TYPE
!     ******************************************************************
!     ECPT( 1) = ELEMENT ID            NECPT(1)       INTEGER
!     ECPT( 2) = GRID POINT 1          NGRID(1)       INTEGER
!     ECPT( 3) = GRID POINT 2          NGRID(2)       INTEGER
!     ECPT( 4) = GRID POINT 3          NGRID(3)       INTEGER
!     ECPT( 5) = GRID POINT 4          NGRID(4)       INTEGER
!     ECPT( 6) = GRID POINT 5          NGRID(5)       INTEGER
!     ECPT( 7) = GRID POINT 6          NGRID(6)       INTEGER
!     ECPT( 8) = GRID POINT 7          NGRID(7)       INTEGER
!     ECPT( 9) = GRID POINT 8          NGRID(8)       INTEGER
!     ECPT(10) = COORD SYS ID-STRESS   ID1            INTEGER
!     ECPT(11) = ANIS. MATERIAL ANGLE  TH             REAL
!     ECPT(12) = MATERIAL ID           MATID1         INTEGER
!     ECPT(13) = THICKNESS             T              REAL
!     ECPT(14) = COORD SYS ID 1        ISYS1          INTEGER
!     ECPT(15) = X1                    X1             REAL
!     ECPT(16) = Y1                    Y1             REAL
!     ECPT(17) = Z1                    Z1             REAL
!     ECPT(18) = COORD SYS ID 2        ISYS2          INTEGER
!     ECPT(19) = X2                    X2             REAL
!     ECPT(20) = Y2                    Y2             REAL
!     ECPT(21) = Z2                    Z2             REAL
!     ECPT(22) = COORD SYS ID 3        ISYS3          INTEGER
!     ECPT(23) = X3                    X3             REAL
!     ECPT(24) = Y3                    Y3             REAL
!     ECPT(25) = Z3                    Z3             REAL
!     ECPT(26) = COORD SYS ID 4        ISYS4          INTEGER
!     ECPT(27) = X4                    X4             REAL
!     ECPT(28) = Y4                    Y4             REAL
!     ECPT(29) = Z4                    Z4             REAL
!     ECPT(30) = COORD SYS ID 5        ISYS5          INTEGER
!     ECPT(31) = X5                    X5             REAL
!     ECPT(32) = Y5                    Y5             REAL
!     ECPT(33) = Z5                    Z5             REAL
!     ECPT(34) = COORD SYS ID 6        ISYS6          INTEGER
!     ECPT(35) = X6                    XL             REAL
!     ECPT(36) = Y6                    Y6             REAL
!     ECPT(37) = Z6                    Z6             REAL
!     ECPT(38) = COORD SYS ID 7        ISYS7          INTEGER
!     ECPT(39) = X7                    X7             REAL
!     ECPT(40) = Y7                    Y7             REAL
!     ECPT(41) = Z7                    Z7             REAL
!     ECPT(42) = COORD SYS ID 8        ISYS8          INTEGER
!     ECPT(43) = X8                    X8             REAL
!     ECPT(44) = Y8                    Y8             REAL
!     ECPT(45) = Z8                    Z8             REAL
!     ECPT(46) = ELEMENT TEMP          TTEMP          REAL
!     ECPT(47) = 0.                    EDT            REAL
!     ECPT(48) = TEMPERATURE SET       ISETNO         INTEGER
!     ECPT(49) = *
!     ECPT(. ) = *  GRID POINT TEMPERATURES
!     ECPT(56) = *
!     ECPT(57) = *
!     ECPT(. ) = *  TRANSLATIONAL DOF-S OF GRIDS FOR THIS ELEMENT
!     ECPT(80) = *
 
 
!     TEST FOR PIVOT POINT
 
 DO  kk = 1,8
   IF (ngrid(kk) == npvt) GO TO 20
 END DO
 
!     IF FALL HERE NO ELEMENT GRID POINT IS THE PIVOT POINT
 
 CALL mesage (-30,34,ecpt(1))
 
!     UNIT I VECTOR IS FROM GRID POINT 1 TO GRID POINT 2
 
 20 DO  i = 1,3
   veci(i) = xy2(i) - xy1(i)
 END DO
 vecil = SQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 IF (vecil == 0.0) GO TO 60
 veci(1) = veci(1)/vecil
 veci(2) = veci(2)/vecil
 veci(3) = veci(3)/vecil
 
!     K VECTOR IS OBTAINED BY CROSSING I INTO VECTOR FROM GRID PT. 1 TO
!     GRID
 
 veck(1) = veci(2)*(z4-z1) - veci(3)*(y4-y1)
 veck(2) = veci(3)*(x4-x1) - veci(1)*(z4-z1)
 veck(3) = veci(1)*(y4-y1) - veci(2)*(x4-x1)
 veckl=SQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 IF (veckl == 0.0) GO TO 60
 veck(1) = veck(1)/veckl
 veck(2) = veck(2)/veckl
 veck(3) = veck(3)/veckl
 
!     J VECTOR IS OBTAINED BY CROSSING K INTO I
 
 vecj(1) = veck(2)*veci(3) - veck(3)*veci(2)
 vecj(2) = veck(3)*veci(1) - veck(1)*veci(3)
 vecj(3) = veck(1)*veci(2) - veck(2)*veci(1)
 
 e1t(1) = veci(1)
 e1t(2) = veci(2)
 e1t(3) = veci(3)
 e1t(4) = vecj(1)
 e1t(5) = vecj(2)
 e1t(6) = vecj(3)
 e1t(7) = veck(1)
 e1t(8) = veck(2)
 e1t(9) = veck(3)
 DO  i = 1,6
   se1t(i) = e1t(i)
 END DO
 
!     STORE ELEMENT COORDS FOR GRIDS 1 AND 2
 
 xx(1) = 0.
 xx(2) = 0.
 xx(3) = vecil
 xx(4) = 0.
 
!     FOR GRIDS 3-8, THE X COORDINATE IS THE DOT PRODUCT OF HTE VECTOR
!     FROM GRID POINT 1 TO THE GRID POINT AND THE I VECTOR. THE Y COORD.
!     IS THE L OF THE I VECTOR CROSSED INTO THE VECTOR FROM GRID 1 TO
!     THE GRID POINT.
 
 DO  i = 3,8
   ixx  = 2*i - 1
   isub = 4*i + 11
   vec(1)  = ecpt(isub  ) - x1
   vec(2)  = ecpt(isub+1) - y1
   vec(3)  = ecpt(isub+2) - z1
   xx(ixx) = vec(1)*veci(1) + vec(2)*veci(2) + vec(3)*veci(3)
   vvec(1) = veci(2)*vec(3) - veci(3)*vec(2)
   vvec(2) = veci(3)*vec(1) - veci(1)*vec(3)
   vvec(3) = veci(1)*vec(2) - veci(2)*vec(1)
   xx(ixx+1) = SQRT(vvec(1)**2 + vvec(2)**2 + vvec(3)**2)
 END DO
 GO TO 70
 
!     INAPPROPRIATE GEOMETRY
 
 60 CALL mesage (30,31,ecpt(1))
 nogo = 1
 RETURN
 
!     COMPUTE MATERIAL PROPERTIES
 
 70 tth   = th*3.1415927/180.
 sinth = SIN(tth)
 costh = COS(tth)
 eltemp= ttemp
 inflag= 2
 matid = matid1
 CALL mat (ecpt(1))
 DO  i = 1,3
   g(i)  = qq(i)
 END DO
 g(4)  = qq(2)
 g(5)  = qq(4)
 g(6)  = qq(5)
 g(7)  = qq(3)
 g(8)  = qq(5)
 g(9)  = qq(6)
 thick = t
 DO  i = 1,9
   r(i)  = g(i)
 END DO
 IF (isetno /= 0) CALL gmmats (r,3,3,0,alphas,3,1,0,st)
 
!     ZERO OUT THE KIJ AND SAVE MATRICES
 
 DO  i = 1,36
   kwd(i) = 0.d0
   kij(i) = 0.d0
 END DO
 DO  i = 1,72
   SAVE(i) = 0.d0
 END DO
 
 pt(1) =-0.57735027D0
 pt(2) =-pt(1)
 h(1)  = 1.d0
 h(2)  = 1.d0
 IF (id1 == 2) GO TO 120
 pt(1) =-0.77459667D0
 pt(2) = 0.d0
 pt(3) =-pt(1)
 h(1)  = 5.d0/9.d0
 h(2)  = 8.d0/9.d0
 h(3)  = h(1)
 
!     2 OR 3 QUADRATURE POINTS
 
 120 DO  iii = 1,id1
   DO  jjj = 1,id1
     
!     COMPUTE GAUSS POINT STRESSES
     
     
!     COMPUTE DERIVATIVES WITH RESPECT TO XI AND ETA
!     EACH GRID POINT
     
     DO  n = 1,4
       dnxi(n)  = .25D0*xi(n)*(1.d0+pt(jjj)*eta(n))*  &
           (2.d0*pt(iii)*xi(n)+pt(jjj)*eta(n))
       dneta(n) = .25D0*eta(n)*(1.d0+pt(iii)*xi(n))*  &
           (pt(iii)*xi(n)+2.d0*pt(jjj)*eta(n))
     END DO
     DO  n = 5,7,2
       
       dnxi(n)  = -pt(iii)*(1.d0+pt(jjj)*eta(n))
       dneta(n) = .5D0*(1.d0-pt(iii)*pt(iii))*eta(n)
     END DO
     
     DO  n = 6,8,2
       dnxi(n)  = .5D0*xi(n)*(1.d0-pt(jjj)*pt(jjj))
       dneta(n) = -pt(jjj)*(1.d0+pt(iii)*xi(n))
     END DO
     
!     COMPUTE JACOBEAN
     
!           N1XI   N2XI   N3XI   N4XI   N5XI   N6XI   N7XI   N8XI
!     DNC = N1ETA  N2ETA  N3ETA  N4ETA  N5ETA  N6ETA  N7ETA  N8ETA
     
!          X1  Y1
!          X2  Y2
!          X3  Y3
!     XX = X4  Y4
!          X5  Y5
!          X6  Y6
!          X7  Y7
!          X8  Y8
     
     CALL gmmatd (dnc,2,8,0,xx,8,2,0,xjb)
     
!     XJB IS ROW-STORED-IT MUST BE COLUMN-STORED AND DOUBLY DIMENSIONED
!     FOR INVERSION
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xxjb(i,j) = xjb(k)
       END DO
     END DO
     
!     COMPUTE INVERSE AND DETERMINANT OF JACOBEAN
     
     CALL inverd (2,xxjb,2,dumarg,0,determ,ising,iws)
     IF (ising == 2) CALL mesage (-30,143,ecpt(1))
     dhh = determ*h(iii)*h(jjj)
     
!     COMPUTE DERIVATIVES WITH RESPECT TO X AND Y
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xjb(k) = xxjb(i,j)
       END DO
     END DO
     CALL gmmatd (xjb,2,2,0,dnc,2,8,0,dnl)
     
!           N1X N2X N3X N4X N5X N6X N7X N8X
!     DNL = N1Y N2Y N3Y N4Y N5Y N6Y N7Y N8Y
     
     DO  i = 1,72
       bb(i) = 0.
     END DO
     
!     SET UP INDICATOR FOR GRID POINT TEMPERATURES
     
     idtemp = 0
     DO  i = 1,8
       IF (tgrid(i) /= 0.) GO TO 200
     END DO
     GO TO 210
     200 idtemp = 1
     
     210 DO  n = 1,8
       
       DO  i = 1,9
         semp(i) = 0.
       END DO
       DO  i = 1,6
         s(i) = 0.
       END DO
       s(1) = dnx(n)
       s(4) = dny(n)
       s(5) = dny(n)
       s(6) = dnx(n)
       
!     TRANSFORM TO ELEMENT COORDINATES
       
       IF (necpt(4*n+10) == 0) GO TO 240
       CALL transs (necpt(4*n+10),ttb)
       CALL gmmats (se1t,2,3,0,ttb,3,3,0,stb)
       GO TO 260
       240 DO  i = 1,6
         stb(i) = se1t(i)
       END DO
       260 CALL gmmats (s,3,2,0,stb,2,3,0,semp(1))
       n3 = 3*n
       bb(n3- 2) = semp(1)
       bb(n3- 1) = semp(2)
       bb(n3   ) = semp(3)
       bb(n3+22) = semp(4)
       bb(n3+23) = semp(5)
       bb(n3+24) = semp(6)
       bb(n3+46) = semp(7)
       bb(n3+47) = semp(8)
       bb(n3+48) = semp(9)
     END DO
     
!     BRING IN G MATRIX
     
     CALL gmmats (r,3,3,0,bb,3,24,0,db)
     
!     COMPUTE STRESSES
     
     CALL gmmats (db,3,24,0,disp,24,1,0,sig)
     
     
!     COMPUTE GAUSS POINT  TEMPERATURES
     
     IF (isetno == 0) GO TO 350
     IF (idtemp == 1) GO TO 280
     rgtemp = eltemp - tref
     GO TO 330
     
!     ALL TEMPERATURES ARE DEFAULT VALUE
     
     280 DO  n = 1,4
       dn(n) = .25*(1.+pt(iii)*xi(n))*(1.+pt(jjj)*eta(n))  &
           *(pt(iii)*xi(n)+pt(jjj)*eta(n)-1.)
     END DO
     DO  n = 5,7,2
       dn(n) = .5*(1.-pt(iii)*pt(iii))*(1.+pt(jjj)*eta(n))
     END DO
     DO  n = 6,8,2
       dn(n) = .5*(1.+pt(iii)*xi(n))*(1.-pt(jjj)*pt(jjj))
     END DO
     gstemp = 0.
     DO  n = 1,8
       gstemp = gstemp + dn(n)*tgrid(n)
     END DO
     rgtemp = gstemp - tref
     330 DO  i = 1,3
       sig(i) = sig(i) - st(i)*rgtemp
     END DO
     
     350 CONTINUE
     
!     FORM KWD MATRIX
     
     kwd( 1) = sig(2)
     kwd( 2) =-sig(3)
     kwd( 7) =-sig(3)
     kwd( 8) = sig(1)
     kwd(15) = sig(1) + sig(2)
     kwd(16) =-sig(3)
     kwd(17) = sig(3)
     kwd(18) = sig(1) - sig(2)
     kwd(21) =-sig(3)
     kwd(27) = sig(3)
     kwd(33) = sig(1) - sig(2)
     
!     FORM CID FOR I = NPVT
     
     DO  i = 1,18
       cid( i) = 0.d0
     END DO
     cid( 3) = dny(kk)
     cid( 6) =-dnx(kk)
     cid( 7) =-.5*dny(kk)
     cid( 8) = .5*dnx(kk)
     cid(10) = dnx(kk)
     cid(14) = dny(kk)
     cid(16) =.5*dny(kk)
     cid(17) =.5*dnx(kk)
     
     CALL gmmatd (cid,6,3,1,kwd,6,6,0,kmult)
     
!     LOOP FOR THE 8 6X6 PARTITIONS CORRESPONDING TO THE PRESENT
!     PIVOT POINT
     
     DO  n = 1,8
       
       DO  i = 1,18
         cjd(i) = 0.d0
       END DO
       
       cjd( 3) = dny(n)
       cjd( 6) =-dnx(n)
       cjd( 7) =-.5*dny(n)
       cjd( 8) = .5*dnx(n)
       cjd(10) = dnx(n)
       cjd(14) = dny(n)
       cjd(16) =.5*dny(n)
       cjd(17) =.5*dnx(n)
       
       CALL gmmatd (kmult,3,6,0,cjd,6,3,0,tempar(1))
       
!     THROW IN JACOBEAN DETERMINANT AND WEIGHT FACTORS
       
       DO  i = 1,9
         tempar(i) = tempar(i)*dhh
       END DO
       
!     ADD THE RESULTS OF THIS INTEGRATION TO THE PREVIOUS RESULTS
       
       ll = 9*(n-1)
       DO  i = 1,9
         l = ll + i
         SAVE(l) = SAVE(l) + tempar(i)
       END DO
       
!     LOOP FOR MORE PARTITIONS
       
     END DO
     
!     LOOP FOR MORE GAUSS POINTS
     
   END DO
 END DO
 
!     CHECK ON NECESSITY OF PRE-MULTIPLYING COORDINATE TRANSFORMATIONS
 
 isub = 4*kk + 10
 IF (necpt(isub) == 0) GO TO 420
 
!     ELEMENT TO GLOBAL
 
 CALL transd (necpt(isub),tb)
 CALL gmmatd (e1t,3,3,0,tb,3,3,0,premul)
 GO TO 440
 420 DO  i = 1,9
   premul(i) = e1t(i)
 END DO
 440 DO  n = 1,8
   ll = 9*n - 8
   CALL gmmatd (premul,3,3,1,SAVE(ll),3,3,0,temp)
   
!     STORE THE 3 X 3 IN TSAVE
   
   DO  i = 1,9
     l = 9*n + i - 9
     tsave(l) = temp(i)
   END DO
   
 END DO
 
!     NOW CHECK ON THE NECESSITY FOR POST-MULTIPLYING TRANSFORMATIONS
 
 DO  n = 1,8
   isub = 4*n + 10
   ll = 9*n - 8
   IF (necpt(isub) == 0) GO TO 470
   
!     GLOBAL TO ELEMENT
   
   CALL transd (necpt(isub),tb)
   CALL gmmatd (e1t,3,3,0,tb,3,3,0,pstmul)
   GO TO 490
   470 DO  i = 1,9
     pstmul(i) = e1t(i)
   END DO
   
!     POST-MULTIPLY
   
   490 CALL gmmatd (tsave(ll),3,3,0,pstmul,3,3,0,temp)
   
!     FILL OUT THE 6 X 6 PARTITION
   
   kij( 1) = temp(1)*thick
   kij( 2) = temp(2)*thick
   kij( 3) = temp(3)*thick
   kij( 7) = temp(4)*thick
   kij( 8) = temp(5)*thick
   kij( 9) = temp(6)*thick
   kij(13) = temp(7)*thick
   kij(14) = temp(8)*thick
   kij(15) = temp(9)*thick
   
!     INSERT INTO THE OVERALL STIFFNESS MATRIX
   
   CALL ds1b (kij,necpt(n+1))
   
!     LOOP FOR MORE PARTITIONS
   
 END DO
 
 RETURN
END SUBROUTINE dis2d8
