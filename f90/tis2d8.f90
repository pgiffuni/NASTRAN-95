SUBROUTINE tis2d8 (temp,pg)
     
!     THIS ROUTINE COMPUTES EQUIVALENT LOADS DUE TO GRID POINT
!     TEMPERATURES FOR THE 2-D, 8 GRID POINT ISOPARAMETRIC ELEMENT
 
!     ECPT LIST                        IN
!                                      THIS
!     ECPT       DESCRIPTION           ROUTINE         TYPE
!     --------   --------------------  --------       -------
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
 
 
 REAL, INTENT(IN OUT)                     :: temp(8)
 REAL, INTENT(OUT)                        :: pg(1)
 REAL :: kx,ky
 DIMENSION       g(9),qq(15),xi(8),eta(8),xy1(3),xy2(3),bt(12),  &
     ecpt(1),tempar(8),dnx(1),dny(1),iz(1),dnxi(1), dneta(1), dn(8),iws(2,3),  &
     vec(3),vvec(3),veci(3),vecj(3),veck(3),e1t(6)
 COMMON /tranx / z(14)
 COMMON /trimex/ necpt(1),ngrid(8),id1,th,matid1,t,isys1,x1,y1,z1,  &
     isys2,x2,y2,z2,isys3,x3,y3,z3,isys4,x4,y4,z4,  &
     isys5,x5,y5,z5,isys6,x6,y6,z6,isys7,x7,y7,z7,  &
     isys8,x8,y8,z8,ttemp,SAVE(16),rtside(3),tempav,  &
     alphas(3),xx(16),dnc(16),dnl(16),xxjb(2,2),  &
     xjb(4),pt(3),h(3),g,bt,determ,dumarg
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     tref,GE,kx,ky,c
 EQUIVALENCE     (ecpt(1),necpt(1)),(z(1),iz(1)),(qq(1),g11),  &
     (dnc(1),dnxi(1)),(dnc(9),dneta(1)), (dnl(1),dnx(1)),(dnl(9),dny(1)),  &
     (tempar(1),bt(1)),(xy1(1),x1),(xy2(1),x2)
 DATA    xi    / -1., 1., 1., -1., 0., 1., 0., -1./
 DATA    eta   / -1.,-1., 1.,  1.,-1., 0., 1.,  0./
 
!     UNIT I VECTOR IS FROM GRID POINT 1 TO GRID POINT 2
 
 DO  i = 1,3
   veci(i) = xy2(i) - xy1(i)
 END DO
 vecil = SQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 IF (vecil == 0.0) GO TO 40
 veci(1) = veci(1)/vecil
 veci(2) = veci(2)/vecil
 veci(3) = veci(3)/vecil
 
!     K VECTOR IS OBTAINED BY CROSSING I INTO VECTOR FROM GRID PT. 1 TO
!     GRID
 
 veck(1) = veci(2)*(z4-z1) - veci(3)*(y4-y1)
 veck(2) = veci(3)*(x4-x1) - veci(1)*(z4-z1)
 veck(3) = veci(1)*(y4-y1) - veci(2)*(x4-x1)
 veckl   = SQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 IF (veckl == 0.0) GO TO 40
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
 
!     STORE ELEMENT COORDS FOR GRIDS 1 AND 2
 
 xx(1) = 0.
 xx(2) = 0.
 xx(3) = vecil
 xx(4) = 0.
 
!     FOR GRIDS 3-8, THE X COORDINATE IS THE DOT PRODUCT OF HTE VECTOR
!     FROM THE GRID POINT TO
!     GRID POINT 1 TO THE GRID POINT AND THE I VECTOR. THE Y COORD. IS
!     THE L OF THE I VECTOR CROSSED INTO THE VECTOR FROM GRID 1 TO THE
!     GRID POINT.
 
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
 GO TO 150
 
!     INAPPROPRIATE GEOMETRY
 
 40 CALL mesage (-30,31,ecpt(1))
 
!     COMPUTE MATERIAL PROPERTIES
 
 150 tth   = th*3.1415927/180.
 sinth = SIN(tth)
 costh = COS(tth)
 inflag= 2
 matid = matid1
 
!     ZERO OUT SOME MATRICES
 
 DO  i = 1,16
   SAVE(i) = 0.
 END DO
 
 pt(1) =-0.57735027
 pt(2) =-pt(1)
 h(1)  = 1.
 h(2)  = 1.
 IF (id1 == 2) GO TO 226
 pt(1) =-0.77459667D0
 pt(2) = 0.d0
 pt(3) =-pt(1)
 h(1)  = 5.0/9.0
 h(2)  = 8.0/9.0
 h(3)  = h(1)
 
!     2 OR 3 QUADRATURE POINTS
 
 226 DO  iii = 1,id1
   DO  jjj = 1,id1
     
!     COMPUTE DERIVATIVES WITH RESPECT TO XI AND ETA
!     EACH GRID POINT
     
     DO  n = 1,4
       dn(n)   = 0.25*(1.+pt(iii)*xi(n))*(1.+pt(jjj)*eta(n))*  &
           (pt(iii)*xi(n)+pt(jjj)*eta(n)-1.)
       dnxi(n) = 0.25*xi(n)*(1.+pt(jjj)*eta(n))*  &
           (2.*pt(iii)*xi(n)+pt(jjj)*eta(n))
       dneta(n)= 0.25*eta(n)*(1.+pt(iii)*xi(n))*  &
           (pt(iii)*xi(n)+2.*pt(jjj)*eta(n))
     END DO
     
     DO  n = 5,7,2
       dn(n)   = 0.5*(1.-pt(iii)*pt(iii))*(1.+pt(jjj)*eta(n))
       dnxi(n) =-pt(iii)*(1.+pt(jjj)*eta(n))
       dneta(n)= 0.5*(1.-pt(iii)*pt(iii))*eta(n)
     END DO
     
     DO  n = 6,8,2
       dn(n)   = 0.5*(1.+pt(iii)*xi(n))*(1.-pt(jjj)*pt(jjj))
       dnxi(n) = 0.5*xi(n)*(1.-pt(jjj)*pt(jjj))
       dneta(n)=-pt(jjj)*(1.+pt(iii)*xi(n))
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
     
     CALL gmmats (dnc,2,8,0, xx,8,2,0, xjb)
     
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
     
     CALL invers (2,xxjb,2,dumarg,0,determ,ising,iws)
     IF (ising == 2) CALL mesage (-30,143,ecpt(1))
     
!     COMPUTE DERIVATIVES WITH RESPECT TO X,Y,AND Z
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xjb(k) = xxjb(i,j)
       END DO
     END DO
     CALL gmmats (xjb,2,2,0, dnc,2,8,0, dnl)
     
!           N1X N2X N3X N4X N5X N6X N7X N8X
!     DNL = N1Y N2Y N3Y N4Y N5Y N6Y N7Y N8Y
     
     coef = determ*h(iii)*h(jjj)
     
!     COMPUTE GAUSS POINT TEMPERATURE
     
     gstemp = 0.
     DO  n = 1,8
       gstemp = gstemp + dn(n)*(temp(n))
     END DO
     
!     GSTEMP IS THE GAUSS POINT TEMPERATURE. FIND MATERIAL PROPERTIES
!     BASED ON THIS TEMPERATURE. IF SAME AS PREVIOUS TEMPERATURE,DO NOT
!     RECOMPUTE.
     
     lll = iii*jjj
     IF (lll == 1) GO TO 236
     IF (gstemp == eltemp) GO TO 237
     236 eltemp = gstemp
     CALL mat (ecpt(1))
     DO  i = 1,3
       g(i) = qq(i)
     END DO
     g(4) = qq(2)
     g(5) = qq(4)
     g(6) = qq(5)
     g(7) = qq(3)
     g(8) = qq(5)
     g(9) = qq(6)
     alphas(1) = alpha1
     alphas(2) = alpha2
     alphas(3) = alp12
     
     CALL gmmats (g,3,3,0, alphas,3,1,0, rtside)
     
!     COMPUTE RELATIVE GAUSS POINT TEMPERATURE
     
     rgtemp = gstemp - tref
     
     237 CONTINUE
     
     coef = coef*rgtemp
     
!     SET UP BT
     
     DO  kk = 1,8
       
       DO  i = 1,6
         bt(i) = 0.
       END DO
       
       bt(1) = dnx(kk)
       bt(3) = dny(kk)
       bt(5) = bt(3)
       bt(6) = bt(1)
       
       CALL gmmats (bt,2,3,0, rtside,3,1,0, tempar(7))
       
!     ADD TO PREVIOUS RESULTS
       
       SAVE(2*kk-1) = SAVE(2*kk-1) + tempar(7)*coef
       SAVE(2*kk  ) = SAVE(2*kk  ) + tempar(8)*coef
       
!     CONTINUE FOR MORE GRID POINTS
       
     END DO
     
!     CONTINUE FOR MORE GAUSS POINTS
     
   END DO
 END DO
 
!     TRANSFORMATIONS AND ADD TO OVERALL VECTOR
 
 DO  kk = 1,8
   tempar(7) = SAVE(2*kk-1)
   tempar(8) = SAVE(2*kk  )
   
!     CONVERT FROM ELEMENT COORDINATES TO BASIC
   
   CALL gmmats (e1t,2,3,1, tempar(7),2,1,0, tempar(1))
   isub = 4*kk + 10
   IF (necpt(isub) == 0) GO TO 300
   
!     MUST TRANSFORM FROM BASIC COORDS TO GLOBAL
   
   CALL basglb (tempar(1),tempar(1),ecpt(isub+1),necpt(isub))
   
!     ADD THIS VECTOR TO OVERALL LOAD VECTOR
   
   300 DO  i = 1,3
     l = ngrid(kk) + i - 1
     pg(l) = pg(l) + tempar(i)*t
   END DO
   
!     GET ANOTHER PARTITION
   
 END DO
 RETURN
END SUBROUTINE tis2d8
