SUBROUTINE is2d8s
     
!     2-D, 8 GRID POINT ISOPARAMETRIC STRUCTURAL  ELEMENT STIFFNESS,
!     MASS, CONDUCTIVITY, AND CAPACITANCE ROUTINE
 
!     SINGLE PRECISION VERSION
 
 LOGICAL :: error,    heat
 INTEGER :: GE,       dict(13), ind6(36), isil(8),  nam(2), scr4,     bcd(2)
 REAL :: kx,       ky,       kxy,      mij,      kij
 DIMENSION       vec(3),   vvec(3),  veci(3),  vecj(3),  veck(3),  &
     xy1(3),   xy2(3),   iz(1),    ecpt(46), qq(15), iws(2,3)
 DIMENSION       g(9),     xi(8),    eta(8),   tb(9),  &
     b(12),    bt(12),   temp(9),  tempar(7),dnx(1),  &
     dny(1),   dnxi(1),  dneta(1), SAVE(144),tsave(216)  &
     ,               e1t(6),   temp1(9), temp2(9), temp3(9), pstmul(9),  &
     premul(6),dn(8),    xx(16),   dnc(16),  dnl(16),  &
     xjb(4),   xxjb(2,2),pt(3),    h(3),     savm(36)
 COMMON /BLANK / skip(16), volume,   surfac
 COMMON /emgest/ necpt(1), ngrid(8), id1,      th,       matid1,  &
     t,        isys1,    x1,       y1,       z1,  &
     isys2,    x2,       y2,       z2,       isys3,  &
     x3,       y3,       z3,       isys4,    x4,  &
     y4,       z4,       isys5,    x5,       y5,  &
     z5,       isys6,    x6,       y6,       z6,  &
     isys7,    x7,       y7,       z7,       isys8,  &
     x8,       y8,       z8,       ttemp,    dumb(119)
 COMMON /matin / matid,    inflag,   eltemp,   stress,   sinth, costh
 COMMON/matout / g11,      g12,      g13,      g22,      g23,  &
     g33,      rho,      alpha1,   alpha2,   alph12, tref,     GE,       dum3(3)
 COMMON /hmtout/ kx,       kxy,      ky,       c
 COMMON /emgdic/ dum2(2),  nlocs,    elid,     iestid
 COMMON /zzzzzz/ z(1)
 COMMON /emgprm/ idum,     jcore,    ncore,    dum12(12),kmbgg(3),  &
     iprec,    error,    heat,     coup
 EQUIVALENCE     (ecpt(1),necpt(1)), (z(1),iz(1)), (temp(1),b(1)),  &
     (dnc(1),dnxi(1))  , (dnc(9),dneta(1)),  &
     (dnl(1),dnx(1))   , (dnl(9),dny(1)), (qq(1),g11),  &
     (tempar(1),bt(1)) , (xy1(1),x1)    , (xy2(1),x2), (isil(1),ngrid(1))
 DATA    xi    / -1.00, 1.00, 1.00,-1.00, 0.00, 1.00, 0.00,-1.00/
 DATA    eta   / -1.00,-1.00, 1.00, 1.00,-1.00, 0.00, 1.00, 0.00/
 DATA    ind6  / 1,7,49,13,55,91,19,61,97,127,25,67,103,133,157,31,  &
     73,109,139,163,181,37,79,115,145,169,187,199,43, 85,121,151,175,193,205,211/
 DATA    nam   , bcd/ 4HIS2D,4H8S  , 4HCIS2,4HD8   /
 DATA    scr4  / 304 /
 
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
!     ECPT(10) = NO. OF GAUSS POINTS   ID1            INTEGER
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
 
 IF (jcore+576 > ncore) CALL mesage (-8,0,nam)
 dict(1) = iestid
 dict(2) = 1
 IF (heat) GO TO 1
 dict(3) = 24
 dict(4) = 7
 nsq     = 576
 GO TO 2
 1 dict(3) = 8
 dict(4) = 1
 nsq     = 64
 
!     SAVE NGRID IN DUMB
 
 2 DO  i = 1,9
   dumb(i) = ecpt(i)
 END DO
 area = 0.0
 
!     SET UP SIL ARRAY SO THAT MATRICES ARE SET UP IN INCREASING SIL
!     ORDER SIL(I)=PARTITION NUMBER OF ITH GRID POINT
 
 i =-8
 5 j = 0
 DO  k = 1,8
   IF (isil(k) < j) CYCLE
   j = isil(k)
   l = k
 END DO
 isil(l) = i
 i = i + 1
 IF (i < 0) GO TO 5
 DO  i = 1,8
   isil(i) =-isil(i)
 END DO
 
 DO  i = 1,nsq
   z(jcore+i) = 0.0
 END DO
 
!     UNIT I VECTOR IS FROM GRID POINT 1 TO GRID POINT 2
 
 DO  i = 1,3
   veci(i) = xy2(i)-xy1(i)
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
 
 e1t(1)  = veci(1)
 e1t(2)  = veci(2)
 e1t(3)  = veci(3)
 e1t(4)  = vecj(1)
 e1t(5)  = vecj(2)
 e1t(6)  = vecj(3)
 
!     STORE ELEMENT COORDS FOR GRIDS 1 AND 2
 
 xx(1) = 0.0
 xx(2) = 0.0
 xx(3) = vecil
 xx(4) = 0.0
 
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
 GO TO 150
 
!     INAPPROPRIATE GEOMETRY
 
 40 CALL mesage (30,31,ecpt(1))
 error = .true.
 
 150 IF (error) RETURN
 
!     SET UP QUADRATURE POINTS AND WEIGHTS
 
 pt(1) =-0.57735027
 pt(2) =-pt(1)
 h(1)  = 1.0
 h(2)  = 1.0
 IF (id1 == 2) GO TO 155
 pt(1) =-0.77459667
 pt(2) = 0.0
 pt(3) =-pt(1)
 h(1)  = 5.0/9.0
 h(2)  = 8.0/9.0
 h(3)  = h(1)
 
 155 IF (heat) GO TO 700
 
!     COMPUTE MATERIAL PROPERTIES
 
 tth   = th*3.1415927/180.
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
 drho  = rho*t
 
!     ZERO THE SAVE MATRICES TO COLLECT INTEGRATIONS
 
 DO  i = 1,36
   savm(i) = 0.0
 END DO
 DO  i = 1,144
   SAVE(i) = 0.0
 END DO
 
!     2 OR 3 QUADRATURE POINTS
 
 DO  iii = 1,id1
   DO  jjj = 1,id1
     
!     COMPUTE DERIVATIVES WITH RESPECT TO XI AND ETA
!     EACH GRID POINT
     
     DO  n = 1,4
       IF (kmbgg(2) /= 0) dn(n) = .25*(1.0+pt(iii)*xi(n))*  &
           (1.0+pt(jjj)*eta(n))*(pt(iii)*xi(n)+pt(jjj)*eta(n)-1.0)
       dnxi(n)  = .25*xi(n)*(1.0+pt(jjj)*eta(n))*  &
           (2.0*pt(iii)*xi(n)+pt(jjj)*eta(n))
       dneta(n) = .25*eta(n)*(1.0+pt(iii)*xi(n))*  &
           (pt(iii)*xi(n)+2.0*pt(jjj)*eta(n))
     END DO
     
     DO  n = 5,7,2
       IF (kmbgg(2) /= 0) dn(n) = .50*(1.0-pt(iii)*pt(iii))*  &
           (1.0+pt(jjj)*eta(n))
       dnxi(n)  = -pt(iii)*(1.0+pt(jjj)*eta(n))
       dneta(n) = .50*(1.0-pt(iii)*pt(iii))*eta(n)
     END DO
     
     DO  n = 6,8,2
       IF (kmbgg(2) /= 0) dn(n) = .50*(1.0+pt(iii)*xi(n))*  &
           (1.0-pt(jjj)*pt(jjj))
       dnxi(n)  = .50*xi(n)*(1.0-pt(jjj)*pt(jjj))
       dneta(n) =-pt(jjj)*(1.0+pt(iii)*xi(n))
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
     
     CALL gmmats (dnc,2,8,0,xx,8,2,0,xjb)
     
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
     IF (ising /= 2) GO TO 241
     CALL mesage (30,143,ecpt(1))
     error =.true.
     RETURN
     
     241 CONTINUE
     dhh  = determ*h(iii)*h(jjj)
     area = area + dhh
     
!     COMPUTE DERIVATIVES WITH RESPECT TO X AND Y
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xjb(k) = xxjb(i,j)
       END DO
     END DO
     CALL gmmats (xjb,2,2,0,dnc,2,8,0,dnl)
     
!           N1X N2X N3X N4X N5X N6X N7X N8X
!     DNL = N1Y N2Y N3Y N4Y N5Y N6Y N7Y N8Y
     
!     SET UP THE BT MATRIX
     
     ic = 0
     DO  kk = 1,8
       IF (kmbgg(1) == 0) GO TO 256
       
       DO  i = 1,12
         bt(i) = 0.0
       END DO
       bt(1) = dnx(kk)
       bt(3) = dny(kk)
       bt(5) = dny(kk)
       bt(6) = dnx(kk)
       
       CALL gmmats (tempar(1),2,3,0,g,3,3,0,tempar(7))
       
!     MULTIPLY G MATRIX BY PRESENT RESULTS
       
!     LOOP FOR THE 8 6X6 PARTITIONS CORRESPONDING TO THE PRESENT
!     PIVOT POINT
       
       256 CONTINUE
       DO  n = kk,8
         ic = ic + 1
         IF (kmbgg(1) == 0) GO TO 281
         
!     SET UP THE B MATRIX
         
         DO  i = 1,12
           b( i) = 0.0
         END DO
         b( 1) = dnx(n)
         b( 4) = dny(n)
         b( 5) = dny(n)
         b( 6) = dnx(n)
!                                    T
!     PERFORM MULTIPLICATION TO GET B *D*B
         
         CALL gmmats (tempar(7),2,3,0,b,3,2,0,tempar(1))
         
!     THROW IN JACOBEAN DETERMINANT AND WEIGHT FACTORS
         
         DO  i = 1,4
           tempar(i) = tempar(i)*dhh
         END DO
         
!     ADD THE RESULTS OF THIS INTEGRATION TO THE PREVIOUS RESULTS
         
         ll = 4*(ic-1)
         DO  i = 1,4
           l = ll + i
           SAVE(l) = SAVE(l) + tempar(i)
         END DO
         281 CONTINUE
         
         IF (kmbgg(2) == 0) GO TO 289
         
         mij = dn(kk)*dn(n)*dhh
         savm(ic) = savm(ic) + mij
         289 CONTINUE
         
!     LOOP FOR MORE PARTITIONS
         
       END DO
     END DO
     
!     LOOP FOR MORE GAUSS POINTS
     
   END DO
 END DO
 IF (kmbgg(2) == 0) GO TO 306
 DO  i = 1,36
   savm(i) = savm(i)*drho
 END DO
 306 CONTINUE
 
!     CHECK ON NECESSITY OF PRE-MULTIPLYING COORDINATE TRANSFORMATIONS
 
 IF (kmbgg(1) == 0) GO TO 500
 ic = 0
 DO  kk = 1,8
   isub = 4*kk + 10
   IF (necpt(isub) == 0) GO TO 310
   
!     ELEMENT TO GLOBAL
   
   CALL transs (necpt(isub),tb)
   CALL gmmats (e1t,2,3,0,tb,3,3,0,premul)
   GO TO 350
   310 DO  i = 1,6
     premul(i) = e1t(i)
   END DO
   350 DO  n = kk,8
     ic = ic + 1
     ll = 4*ic - 3
     CALL gmmats (premul,2,3,1,SAVE(ll),2,2,0,temp)
     
!     STORE THE 3 X 2 IN TSAVE
     
     DO  i = 1,6
       l = 6*ic + i - 6
       tsave(l) = temp(i)
     END DO
     
   END DO
 END DO
 
!     NOW CHECK ON THE NECESSITY FOR POST-MULTIPLYING TRANSFORMATIONS
 
 ic = 0
 DO  n = 1,8
   isub = 4*n + 10
   IF (necpt(isub) == 0) GO TO 410
   
!     GLOBAL TO ELEMENT
   
   CALL transs (necpt(isub),tb)
   CALL gmmats (e1t,2,3,0,tb,3,3,0,pstmul)
   GO TO 450
   410 DO  i = 1,6
     pstmul(i) = e1t(i)
   END DO
   
!     POST-MULTIPLY
   
!     IND6 GIVES STARTING POSITIONS OF VERTICAL 3X3 PARTITIONS, SINCE
!     THE NTH COLUMN MULTIPLIES INTO THE NTH POST-MULTIPLIER
   
   450 DO  m = 1,n
     ic = ic + 1
     ll = ind6(ic)
     CALL gmmats (tsave(ll),3,2,0,pstmul,2,3,0,temp)
     DO  i =  1,9
       temp(i) = temp(i)*thick
     END DO
     
!     PICK UP ROW AND COLUMN PARTITION NUMBERS AND CONVERT TO STARTING
!     POINTS IN OPEN CORE FOR THIS PARTITION  AND ITS TRANSPOSE.
!     TEMP IS PUT INTO ONE PARTITION AND TEMP-TRANSPOSE INTO THE OTHER
     
     ncol = isil(n)
     nrow = isil(m)
     CALL insert (ncol,nrow,3,8,jcore,z,z,temp,temp,iprec)
     
!     LOOP FOR ANOTHER PARTITION FOR THIS POST-MULTIPLIER
     
   END DO
   
!     LOOP FOR ANOTHER POST-MULTIPLIER
   
 END DO
 
!     ADD TO DICTIONARY
 
 dict(5) = GE
 CALL emgout (z(jcore),z(jcore),nsq,1,dict,1,iprec)
 
 500 IF (kmbgg(2) == 0) GO TO 1000
 
 ic = 0
 DO  kk = 1,8
   DO  n = kk,8
     ic = ic + 1
     DO  i = 1,9
       temp(i) = 0.0
     END DO
     
!     CHECK ON TEH NECESSITY OF COORDINATE TRANSFORMATIONS.
!     SINCE EACH PARTITION IS A MULTIPLE OF A 3X3 IDENTITY AND SINCE
!     THE TRANSFORAMATION MATRICES ARE ORTHOGONAL, NO EXPLICIT
!     TRANSFORMA-TIONS FROM THE ELEMENT COORDINATE SYSTEM ARE REQUIRED.
!     ALSO, NO TRANSFORAMTION IS REQUIRED IF TRANSFORMATION MATRICES ARE
!     THE SAME FOR THE GRIDS CORRESPONDING TO THE THE ROW AND COLUMN
     
     term  = savm(ic)
     IF (kk == n) GO TO 570
     isub  = 4*kk + 10
     isub1 =4*n + 10
     IF (necpt(isub) == 0 .AND. necpt(isub1) == 0) GO TO 570
     IF (necpt(isub) == 0) GO TO 520
     CALL transs (necpt(isub),temp1)
     IF (necpt(isub1) == 0) GO TO 530
     520 CALL transs (necpt(isub1),temp2)
     IF (necpt(isub) == 0) GO TO 550
     
!     MULTIPLY THE TRANSFORMATION MATRICES
     
     CALL gmmats (temp1,3,3,1,temp2,3,3,0,temp3)
     GO TO 580
     530 temp3(1) = temp1(1)
     temp3(2) = temp1(4)
     temp3(3) = temp1(7)
     temp3(4) = temp1(2)
     temp3(5) = temp1(5)
     temp3(6) = temp1(8)
     temp3(7) = temp1(3)
     temp3(8) = temp1(6)
     temp3(9) = temp1(9)
     GO TO 580
     550 DO  i = 1,9
       temp3(i) = temp2(i)
     END DO
     GO TO 580
     570 temp(1) = term
     temp(5) = term
     temp(9) = term
     GO TO 600
     
     580 DO  i = 1,9
       temp(i) = term*temp3(i)
     END DO
     
     600 nrow = isil(kk)
     ncol = isil(n)
     CALL insert (ncol,nrow,3,8,jcore,z,z,temp,temp,iprec)
   END DO
 END DO
 
 CALL emgout (z(jcore),z(jcore),nsq,1,dict,2,iprec)
 GO TO 1000
 
!     HEAT FORMULATION
 
!     COMPUTE MATERIAL PROPERTIES
 
 700 sinth  = 0.
 costh  = 0.
 eltemp = ttemp
 inflag = 2
 matid  = matid1
 CALL hmat (ecpt(1))
 thick  = t
 dc     = c*t
 
!     ZERO OUT THE SAVE MATRIX
 
 DO  i = 1,36
   savm(i) = 0.0
   SAVE(i) = 0.0
 END DO
 
 DO  iii = 1,id1
   DO  jjj = 1,id1
     
!     COMPUTE DERIVATIVES WITH RESPECT TO XI AND ETA
!     EACH GRID POINT
     
     DO  n = 1,4
       IF (kmbgg(3) /= 0) dn(n) = .25*(1.0+pt(iii)*xi(n))*  &
           (1.0+pt(jjj)*eta(n))*(pt(iii)*xi(n)+pt(jjj)*eta(n)-1.0)
       dnxi(n)  = .25*xi(n)*(1.0+pt(jjj)*eta(n))*  &
           (2.0*pt(iii)*xi(n)+pt(jjj)*eta(n))
       dneta(n) = .25*eta(n)*(1.0+pt(iii)*xi(n))*  &
           (pt(iii)*xi(n)+2.0*pt(jjj)*eta(n))
     END DO
     
     DO  n = 5,7,2
       IF (kmbgg(3) /= 0) dn(n) = .50*(1.0-pt(iii)*pt(iii))*  &
           (1.0+pt(jjj)*eta(n))
       dnxi(n)  = -pt(iii)*(1.0+pt(jjj)*eta(n))
       dneta(n) = .50*(1.0-pt(iii)*pt(iii))*eta(n)
     END DO
     
     DO  n = 6,8,2
       IF (kmbgg(3) /= 0) dn(n) = .50*(1.0+pt(iii)*xi(n))*  &
           (1.0-pt(jjj)*pt(jjj))
       dnxi(n)  = .50*xi(n)*(1.0-pt(jjj)*pt(jjj))
       dneta(n) = -pt(jjj)*(1.0+pt(iii)*xi(n))
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
     
     CALL gmmats (dnc,2,8,0,xx,8,2,0,xjb)
     
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
     IF (ising /= 2) GO TO 741
     CALL mesage (30,143,ecpt(1))
     error =.true.
     RETURN
     
     741 CONTINUE
     dhh = determ*h(iii)*h(jjj)
     
!     COMPUTE DERIVATIVES WITH RESPECT TO X,Y,AND Z
     
     k = 0
     DO  i = 1,2
       DO  j = 1,2
         k = k + 1
         xjb(k) = xxjb(i,j)
       END DO
     END DO
     CALL gmmats (xjb,2,2,0,dnc,2,8,0,dnl)
     
!           N1X N2X N3X N4X N5X N6X N7X N8X
!     DNL = N1Y N2Y N3Y N4Y N5Y N6Y N7Y N8Y
     
!     SET UP THE BT MATRIX
     
     ic = 0
     DO  kk = 1,8
       IF (kmbgg(1) == 0) GO TO 800
       bt(1) = kx *dnx(kk) + kxy*dny(kk)
       bt(2) = kxy*dnx(kk) + ky*dny(kk)
       
!     DO NOT TRANSFORM FROM MATERIAL COORD SYSTEM TO BASIC AND GLOBAL
!     SINCE THIS IS A SCALAR PROBLEM
       
       800 DO  n = kk,8
         ic = ic + 1
         
!     SET UP THE B MATRIX
         
         IF (kmbgg(1) == 0) GO TO 810
         b(1) = dnx(n)
         b(2) = dny(n)
         
!     O.K. NOW PERFORM FINAL MULTIPLICATION
         
         kij = bt(1)*b(1) + bt(2)*b(2)
         
!     THROW IN JACOBEAN DETERMINANT AND WEIGHT FACTORS
         
         kij = kij*dhh
         
!     ADD THE RESULTS OF THIS INTEGRATION TO PREVIOUS RESULTS
         
         SAVE(ic) = SAVE(ic) + kij
         810 IF (kmbgg(3) == 0) CYCLE
         bij = dn(kk)*dn(n)*dhh
         savm(ic) = savm(ic) + bij
         
!     LOOP FOR MORE PARTITIONS
         
       END DO
     END DO
     
!     LOOP FOR ADDITIONAL GAUSS POINTS
     
   END DO
 END DO
 
 DO  i = 1,36
   IF (kmbgg(1) /= 0) SAVE(i) = SAVE(i)*thick
   IF (kmbgg(3) /= 0) savm(i) = savm(i)*dc
 END DO
 
!     INSERT INTO OVERALL STIFFNESS MATRIX
 
 ic = 0
 DO  i = 1,8
   DO  j = i,8
     ic = ic + 1
     nrow = isil(i)
     ncol = isil(j)
     IF (kmbgg(1) /= 0) CALL insert (ncol,nrow,1,8,jcore,z,z,  &
         SAVE(ic),SAVE(ic),iprec)
     IF (kmbgg(3) /= 0) CALL insert (ncol,nrow,1,8,jcore+64,z,z,  &
         savm(ic),savm(ic),iprec)
   END DO
 END DO
 
 IF (kmbgg(1) /= 0) CALL emgout (z(jcore),z(jcore),nsq,1,dict,1, iprec)
 IF (kmbgg(3) /= 0) CALL emgout (z(jcore+64),z(jcore+64),nsq,1, dict,3,iprec)
 GO TO 5000
 
!     SAVE ELEMENT NAME, ID, THICKNESS, DENSITY, NO. OF GRID POINTS,
!     GRID POINT DATA, AND AREA IF USER REQUESTED VOLUME AND AREA
!     COMPUTATION
 
 1000 IF (volume <= 0.0 .AND. surfac <= 0.0) GO TO 5000
 ecpt(2) = ecpt(13)
 ecpt(3) = rho
 j = 4
 necpt(j) = 8
 ecpt(46) = area
 CALL WRITE (scr4,bcd,2,0)
 CALL WRITE (scr4,ecpt(1),4,0)
 CALL WRITE (scr4,dumb(2),8,0)
 CALL WRITE (scr4,ecpt(14),33,1)
 5000 RETURN
END SUBROUTINE is2d8s
