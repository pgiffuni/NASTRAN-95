SUBROUTINE qdmm1s
     
!     THIS SUBROUTINE COMPUTES THE STIFFNESS AND MASS MATRIX FOR THE
!     FIRST QUADRILATERAL MEMBRANE ELEMENT.
 
!     SINGLE PRECISION VERSION
 
!     ECPT LIST
!                                                   IN THIS
!        ECPT       DESCRIPTION                     ROUTINE    TYPE
!     ========   ================================   ========  =======
!     ECPT( 1) = ELEMENT ID                         NECPT(1)  INTEGER
!     ECPT( 2)   GRID POINT A                       NGRID(1)  INTEGER
!     ECPT( 3)   GRID POINT B                       NGRID(2)  INTEGER
!     ECPT( 4)   GRID POINT C                       NGRID(3)  INTEGER
!     ECPT( 5)   GRID POINT D                       NGRID(4)  INTEGER
!     ECPT( 6) = THETA = ANGLE OF MATERIAL          ANGLE     REAL
!     ECPT( 7)   MATERIAL ID                        MATID     INTEGER
!     ECPT( 8) = THICKNESS                          T         REAL
!     ECPT( 9) = NON-STRUCTURAL MASS                FMU       REAL
!     ECPT(10)   COORD. SYSTEM ID 1                 NECPT(10) INTEGER
!     ECPT(11) = X1                                 X1        REAL
!     ECPT(12) = Y1                                 Y1        REAL
!     ECPT(13) = Z1                                 Z1        REAL
!     ECPT(14)   COORD. SYSTEM ID 2                 NECPT(14) INTEGER
!     ECPT(15) = X2                                 X2        REAL
!     ECPT(16) = Y2                                 Y2        REAL
!     ECPT(17) = Z2                                 Z2        REAL
!     ECPT(18)   COORD. SYSTEM ID 3                 NECPT(18) INTEGER
!     ECPT(19) = X3                                 X3        REAL
!     ECPT(20) = Y3                                 Y3        REAL
!     ECPT(21) = Z3                                 Z3        REAL
!     ECPT(22)   COORD. SYSTEM ID 4                 NECPT(22) INTEGER
!     ECPT(23) = X4                                 X4        REAL
!     ECPT(24) = Y4                                 Y4        REAL
!     ECPT(25)   Z4                                 Z4        REAL
!     ECPT(26) = ELEMENT TEMPERATURE                ELTEMP    REAL
 
 LOGICAL :: nogo,     heat,     planar
 INTEGER :: outpt,    dict(9),  map(2,4), elid,     estid
 REAL :: kij,      la,       lb,       lc,       ld,  &
     lbd1,     lcd1,     lcd2,     ldd2,     mgg(4),  &
     magi,     magj,     magk,     eta01(2), ecpt(26), tempar(144)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,      uwm,      uim
 COMMON /system/ ksystm,   outpt
 COMMON /condas/ consts(4),degra
 COMMON /emgest/ necpt(1), ngrid(4), angle,    matid1,   thick,  &
     fmu,      dummy1,   x1,       y1,       z1,  &
     dummy2,   x2,       y2,       z2, dummy3,   x3,       y3,       z3,  &
     dummy4,   x4,       y4,       z4, dum(75)
 COMMON /emgprm/ dum2(16), mass,     dum3,     iprec,    nogo, heat
 COMMON /matin / matid,    inflag,   eltemp,   stress,   sinth, costh
 COMMON /matout/ g11,      g12,      g13,      g22,      g23,  &
     g33,      rho,      alpha1,   alpha2,   alp12,  &
     tsub0,    gsube,    sigten,   sigcom,   sigshe, g2x211,   g2x212,   g2x222
 COMMON /sma1dp/ tie(9,4), kij(3,3), b(144),   e(9),     etj(9,4)
 COMMON /sma2dp/ u(64),    c(6),     aq(24),   bq(24),   cq(30),  &
     ti(9),    btxk(96)
 COMMON /emgdic/ dmmm(2),  nlocs,    elid,     estid
 EQUIVALENCE     (dict5,dict(5)),    (ecpt(1),necpt(1)), (u(1),tempar(1))
 
 o(d,v,f,h,p,q,y4a,x12,y34,y3a,x23,x14,eta,tea) = (d+(v*tea) +  &
     (f*eta) + (h*tea*eta) + (p*tea*tea) + (q*eta*eta))/  &
     ((-y4a*x12) + (-y34*x12*eta) + ((-y4a*x23) + (y3a*x14))*tea)
 
 eta = 1.0
 tea = 1.0
 IF (heat) GO TO 450
 eta01(1) = 0.211324865
 eta01(2) = 0.788675135
 
!     COMPUTE DIFFERENCES OF COORDINATES OF ACTUAL GRID POINTS
 
 x21 = x2 - x1
 y21 = y2 - y1
 z21 = z2 - z1
 x31 = x3 - x1
 y31 = y3 - y1
 z31 = z3 - z1
 x41 = x4 - x1
 y41 = y4 - y1
 z41 = z4 - z1
 x42 = x4 - x2
 y42 = y4 - y2
 z42 = z4 - z2
 
!     COMPUTE ELEMENTS OF THE E MATRIX
 
 pk1 = y31*z42 - z31*y42
 pk2 = z31*x42 - x31*z42
 pk3 = x31*y42 - y31*x42
 magk= SQRT(pk1**2 + pk2**2 + pk3**2)
 IF (magk <= 1.0E-6) GO TO 410
 pk1 = pk1/magk
 pk2 = pk2/magk
 pk3 = pk3/magk
 
!     HH IS THE MEASURE OF NON-PLANARITY OF THE ELEMENT
 
 hh  = x21*pk1 + y21*pk2 + z21*pk3
 pi1 = x21 - hh*pk1
 pi2 = y21 - hh*pk2
 pi3 = z21 - hh*pk3
 magi= SQRT(pi1**2 + pi2**2 + pi3**2)
 IF (magi <= 1.0E-6) GO TO 420
 pi1 = pi1/magi
 pi2 = pi2/magi
 pi3 = pi3/magi
 hh  =-hh/2.0
 
!     THIS SIGN CHANGE MADE BECAUSE SIGN OF H AS DEFINED ON
!     PAGE 4.87-105 OF PROGRAMMERS MANUAL IS WRONG
 
 temp  = SQRT(x31**2 + y31**2 + z31**2)
 ysub4 = SQRT(x42**2 + y42**2 + z42**2)
 h1    = (2.0*hh)/(temp + ysub4)
 planar= .true.
 IF (h1 > 1.0E-6) planar = .false.
 IF (h1 >= 1.0E-2) WRITE (outpt,28) uim,h1,necpt(1)
 28 FORMAT (a29,' 3061, THE MEASURE OF NON-PLANARITY IS ',e13.5,  &
     ' FOR ELEMENT NUMBER',i9)
 pj1 = pk2*pi3 - pk3*pi2
 pj2 = pk3*pi1 - pk1*pi3
 pj3 = pk1*pi2 - pk2*pi1
 magj= SQRT(pj1**2 + pj2**2 + pj3**2)
 IF (magj <= 1.0E-6) GO TO 430
 pj1 = pj1/magj
 pj2 = pj2/magj
 pj3 = pj3/magj
 
!  *  SET UP E MATRIX (3X3) FOR QUAD-MEMBRANE PROJECTION ONTO MEAN PLANE
!     E IS TRANSPOSE OF E MATRIX IN THEORETICAL MANUAL
 
!     E(1),E(4),E(7) IS I-VECTOR
!     E(2),E(5),E(8) IS J-VECTOR
!     E(3),E(6),E(9) IS K-VECTOR
 
 e(1) = pi1
 e(2) = pj1
 e(3) = pk1
 e(4) = pi2
 e(5) = pj2
 e(6) = pk2
 e(7) = pi3
 e(8) = pj3
 e(9) = pk3
 
!     COMPUTE DIFFERENCES OF COORDINATES OF GRID POINTS IN THE MEAN PLAN
 
 x12 =-(x21*e(1) + y21*e(4) + z21*e(7))
 x13 =-(x31*e(1) + y31*e(4) + z31*e(7))
 x24 =-(x42*e(1) + y42*e(4) + z42*e(7))
 x14 = x12 + x24
 y3a = x31*e(2) + y31*e(5) + z31*e(8)
 y4a = x42*e(2) + y42*e(5) + z42*e(8)
 x34 = x14 - x13
 y34 = y3a - y4a
 x23 = x13 - x12
 IF (y3a <= 0. .OR. y4a <= 0.) GO TO 430
 temp  = x12 + x23*(y4a/y3a)
 ysub4 = (y3a/y4a)*x14
 
!                                              0
!     CHECK FOR INTERNAL ANGLE GREATER THAN 180
 
 IF (x13 >= ysub4 .OR. x14 <= temp) GO TO 430
 
!     GET MASS MATRIX DIAGONALS
 
 IF( mass == 0) GO TO 60
 inflag = 4
 matid  = matid1
 CALL mat (ecpt(1))
 
!     GET TRIANGULAR AREA TIMES TWO
 
 at1 = -x12*y4a
 at2 = -x12*y3a
 at3 = -x23*y4a + x24*y3a
 at4 = -x13*y4a + x14*y3a
 
 fact   = (fmu + g11*thick)/12.0
 mgg(1) = (at4 + at1 + at2)*fact
 mgg(2) = (at1 + at2 + at3)*fact
 mgg(3) = (at2 + at3 + at4)*fact
 mgg(4) = (at3 + at4 + at1)*fact
 
!     COMPUTE LENGTHS OF SIDES OF ELEMENT IN THE MEAN PLANE
 
 60 la = ABS(x12)
 lb = SQRT(x23**2 + y3a**2)
 lc = SQRT(x34**2 + y34**2)
 ld = SQRT(x14**2 + y4a**2)
 IF (la == 0. .OR. lb == 0. .OR. lc == 0. .OR. ld == 0.) GO TO 430
 
!     COMPUTE THE CHARACTERISTIC ANGLES OF ELEMENT IN THE MEAN PLANE
 
 IF (planar) GO TO 75
 cth1  =-x14/ld
 sth1  = y4a/ld
 cth2  = x23/lb
 sth2  = y3a/lb
 cth31 = x34/lc
 sth31 =-y34/lc
 cth41 = cth1
 sth41 = sth1
 cth32 = sth2
 sth32 = cth2
 cth42 = sth31
 sth42 = cth31
 dlt1  = cth31*cth32 - sth31*sth32
 dlt2  = cth42*cth41 + sth41*sth42
 ldd2  = ld*dlt2
 lbd1  = lb*dlt1
 lcd1  = lc*dlt1
 lcd2  = lc*dlt2
 
!     SET UP THE (12X8)  TRANSFORMATION MATRIX B BETWEEN THE MEAN PLANE
!                        AND ACTUAL GRID POINTS
 
 do70 i = 2,92
 b(i) = 0.0
 70 CONTINUE
 
 b( 1) = 1.0
 b(10) = 1.0
 b(17) =-hh/la
 b(18) =-hh/(ld*sth1) + ((hh*cth1)/(la*sth1))
 b(19) = hh/la
 b(20) = (hh*cth2)/(la*sth2)
 b(23) = (hh*cth42)/ldd2
 b(24) = (hh*sth42)/ldd2
 b(27) = 1.0
 b(36) = 1.
 b(41) =-b(17)
 b(42) = (-hh*cth1)/(la*sth1)
 b(43) = b(17)
 b(44) = ((-hh*cth2)/(la*sth2)) + (hh/(lb*sth2))
 b(45) = (-hh*sth31)/lbd1
 b(46) = (-hh*cth31)/lbd1
 b(53) = 1.
 b(62) = 1.
 b(68) =-hh/(lb*sth2)
 b(69) = hh*((sth31/lbd1) + (cth32/lcd1))
 b(70) = hh*((cth31/lbd1) + (sth32/lcd1))
 b(71) = (-hh*sth41)/lcd2
 b(72) = (hh*cth41)/lcd2
 b(79) = 1.0
 b(88) = 1.0
 b(90) = hh/(ld*sth1)
 b(93) = (-hh*cth32)/lcd1
 b(94) = (-hh*sth32)/lcd1
 b(95) = hh*((-cth42/ldd2) + (sth41/lcd2))
 b(96) = hh*((-sth42/ldd2) - (cth41/lcd2))
 
 75 theta = angle*degra
 sinth = SIN(theta)
 costh = COS(theta)
 IF (ABS(sinth) < 1.0E-06) sinth = 0.0E0
 eltemp = ecpt(26)
 inflag = 2
 matid  = matid1
 
!                                                     T
!     COMPUTE TRANSFORMED MATRIX OF STIFFNESSES  C = P  * G * P
 
 CALL mat (ecpt(1))
 
!     STORE INTO G MATRIX
 
 c(1) = g11
 c(2) = g12
 c(3) = g22
 c(4) = g13
 c(5) = g23
 c(6) = 0.
 fact = g33*thick/(x24*y3a - x13*y4a)*2.
 
!     COMPUTE COEFFICIENTS OF THE GENERAL INTEGRAL
 
!                                            2         2
!     D + E*ETA + F*ZETA + H*ETA*ZETA + P*ETA  + Q*ZETA
!     --------------------------------------------------
!     Y *X   +Y  *X  *ZETA + (Y *X   - Y *X  ) * ETA
!      4  21   34  21          4  32    3  41
 
 aq( 1) =-y4a
 aq( 3) =-x24
 aq( 5) =-x24
 aq( 6) =-y4a
 aq( 7) = y4a
 aq( 9) = x14
 aq(11) = x14
 aq(12) = y4a
 aq(13) = 0.0
 aq(15) = 0.0
 aq(17) = 0.0
 aq(18) = 0.0
 aq(19) = 0.0
 aq(21) =-x12
 aq(23) =-x12
 aq(24) = 0.0
 
 bq( 1) = y3a
 bq( 3) = x23
 bq( 5) = x23
 bq( 6) = y3a
 bq( 7) =-y4a
 bq( 9) =-x14
 bq(11) =-x14
 bq(12) =-y4a
 bq(13) = y4a
 bq(15) = x14
 bq(17) = x14
 bq(18) = y4a
 bq(19) =-y3a
 bq(21) =-x23
 bq(23) =-x23
 bq(24) =-y3a
 
 cq( 1) =-y34
 cq( 3) = x34
 cq( 5) = x34
 cq( 6) =-y34
 cq( 7) = y34
 cq( 9) =-x34
 cq(11) =-x34
 cq(12) = y34
 cq(13) = 0.0
 cq(15) =-x12
 cq(17) =-x12
 cq(18) = 0.0
 cq(19) = 0.0
 cq(21) = x12
 cq(23) = x12
 cq(24) = 0.0
 
 nn = 0
 DO  i = 1,4
   DO  k = 1,2
     DO  j = 1,4
       DO   l = 1,2
         nn  = nn + 1
         im1 = i  - 1
         jm1 = j  - 1
         km1 = k  - 1
         lm1 = l  - 1
         k1  = 6*im1 + 4*km1 + 1
         k2  = 6*im1 + 3*km1 + 3
         l1  = 6*jm1 + 4*lm1 + 1
         l2  = 6*jm1 + 3*lm1 + 3
         kl  = k + l - 1
         k3  = k + 3
         l3  = l + 3
         d = c(kl)*aq(k1)*aq(l1)+c(k3)*aq(k1)*aq(l2)+c(l3)*aq(k2)*aq(l1)
         
         v = c(kl)*((aq(k1)*bq(l1))+(bq(k1)*aq(l1)))+c(k3)*((aq(k1)*bq(l2))  &
             + (bq(k1)*aq(l2)))+c(l3)*((aq(k2)*bq(l1))+(bq(k2)*aq(l1)))
         
         f = c(kl)*((aq(k1)*cq(l1))+(cq(k1)*aq(l1)))+c(k3)*((aq(k1)*cq(l2))  &
             + (cq(k1)*aq(l2)))+c(l3)*((aq(k2)*cq(l1))+(cq(k2)*aq(l1)))
         
         h = c(kl)*((bq(k1)*cq(l1))+(cq(k1)*bq(l1)))+c(k3)*((bq(k1)*cq(l2))  &
             + (cq(k1)*bq(l2)))+c(l3)*((bq(k2)*cq(l1))+(cq(k2)*bq(l1)))
         
         p = c(kl)*bq(k1)*bq(l1)+c(k3)*bq(k1)*bq(l2)+c(l3)*bq(k2)*bq(l1)
         
         q = c(kl)*cq(k1)*cq(l1)+c(k3)*cq(k1)*cq(l2)+c(l3)*cq(k2)*cq(l1)
         
!     USE GAUSSIAN INTEGRATION TO FIND THE PARTITIONS OF
!     THE STIFFNESS MATRIX FOR THE MEAN PLANE ELEMENT
         
         u(nn) = 0.0
         DO  ia01 = 1,2
           DO  ja01 = 1,2
             u(nn) = u(nn) +  &
                 o(d,v,f,h,p,q,y4a,x12,y34,y3a,x23,x14,eta01(ia01),eta01(ja01))
           END DO
         END DO
         u(nn) = u(nn)/4.0*thick
         
!     ADD SHEAR TERMS HERE
         u(nn) = u(nn) + fact*(aq(k2)+0.5*(bq(k2)+cq(k2)))  &
             *(aq(l2)+0.5*(bq(l2)+cq(l2)))
       END DO
     END DO
   END DO
 END DO
 
!     TRANSFORM FROM MEAN PLANE TO ACTUAL GRID POINTS
 
!                   T
!      K = B * K * B
 
!     EXPAND MATRIX TO INCLUDE Z COORDINATES
!     IF NON-PLANAR,
 
 IF (planar) GO TO 130
 CALL gmmats (b(1),12,8,0, u(1),8,8,0, btxk(1))
 CALL gmmats (btxk(1),12,8,0, b(1),12,8,1, tempar(1))
 GO TO 200
 
!  *  IF PLANAR, TEMPAR(12X12) .EQ. U(8X8)
 
 130 ij1 =-12
 i2  = 144
 DO  i = 1,64
   tempar(i2+i) = u(i)
 END DO
 DO  i = 1,12
   ij1 = ij1 + 12
   IF (MOD(i,3) /= 0) GO TO 160
   DO  j = 1,12
     ij = ij1 + j
     tempar(ij) = 0.0
   END DO
   CYCLE
   160 DO  j = 1,12
     ij = ij1 + j
     IF (MOD(j,3) /= 0) GO TO 170
     tempar(ij) = 0.0
     CYCLE
     170 i2 = i2 + 1
     tempar(ij) = tempar(i2)
   END DO
 END DO
 
!                T            T
!  *  GENERATE (T  * E) AND (E  * T )
!                I                 J
 
 200 DO  i = 1,4
   ka = 4*i + 6
   IF (necpt(ka) == 0) GO TO 210
   CALL transs (necpt(ka),ti)
   CALL gmmats (ti,3,3,1, e,3,3,0, tie(1,i))
   CALL gmmats (e,3,3,1, ti,3,3,0, etj(1,i))
   CYCLE
   210 DO  ii = 1,9
     tie(ii,i) = e(ii)
   END DO
   etj(1,i) = e(1)
   etj(2,i) = e(4)
   etj(3,i) = e(7)
   etj(4,i) = e(2)
   etj(5,i) = e(5)
   etj(6,i) = e(8)
   etj(7,i) = e(3)
   etj(8,i) = e(6)
   etj(9,i) = e(9)
 END DO
!                                      T              T
!     COMPUTE STIFFNESS MATRIX  K   = T  * E * S   * E  * T
!                                IJ    I        IJ         J
 
!     EXTRACT 3 BY 3 PARTITIONS, TRANSFORM TO GLOBAL, AND INSERT BY
!     ORDER OF SILS INTO A 12 BY 12 MATRIX
 
 DO  i = 1,4
   j = ngrid(i)
   DO  k = 2,5
     IF (necpt(k) == j) GO TO 250
   END DO
   CALL mesage (-30,34,ecpt(1))
   250 map(1,i) = j
   map(2,i) = i
 END DO
 CALL sort (0,0,2,1,map(1,1),8)
 
!     REPLACE SILS WITH INDICES
!     RESORT FOR ORIGINAL ORDER - WORD 1 WILL CONTAIN NEW LOCATION
 
 DO  i = 1,4
   map(1,i) = i
 END DO
 CALL sort (0,0,2,2,map(1,1),8)
 
!     MOVE AND TRANSFORM HERE
!     ROW LOOP
 
 DO  i = 1,4
   ior = 36*(i-1)
   inr = 36*(map(1,i) - 1)
   
!     COLUMN LOOP
   
   DO  j = 1,4
     iocl = ior + 3*(j-1)
     incl = inr + 3*(map(1,j) - 1)
     
!     INNER LOOPS
     
     DO  k = 1,3
       kl = iocl + 12*(k-1)
       DO  l = 1,3
         kij(l,k) = tempar(kl+l)
       END DO
     END DO
     
!     TRANSFORM 3 BY 3
     
     CALL gmmats (kij,3,3,0, etj(1,j),3,3,0, e)
     CALL gmmats (tie(1,i),3,3,0, e,3,3,0, kij)
     
!     INSERT
     
     DO  k = 1,3
       kl = incl + 12*(k-1)
       DO  l = 1,3
         b(kl+l) = kij(l,k)
       END DO
     END DO
   END DO
 END DO
 
!     INSERT WHOLE 12 BY 12 USING EMGOUT
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 12
 dict(4) = 7
 dict5   = gsube
 ldata   = 144
 ieoe    = 1
 ifile   = 1
 CALL emgout (b,b,ldata,ieoe,dict,ifile,iprec)
 
!     DO MASS IF NECESSARY
 
 IF (mass == 0) RETURN
 DO  i = 1,4
   kl = 3*(map(1,i) - 1)
   DO  j = 1,3
     b(kl+j) = mgg(i)
   END DO
 END DO
 dict(2) = 2
 dict(5) = 0
 ldata   = 12
 ifile   = 2
 CALL emgout (b,b,ldata,ieoe,dict,ifile,iprec)
 RETURN
 
!     ERROR EXITS
 
 410 j = 32
 GO TO 440
 420 j = 31
 GO TO 440
 430 j = 26
 440 CALL mesage (30,j,ecpt(1))
 nogo = .true.
 RETURN
 
 450 WRITE  (outpt,460) uwm,necpt(1)
 460 FORMAT (a25,' 3115, QDMM1S FINDS ELEMENT NO.',i9,' PRESENT IN A',  &
     ' HEAT FORMULATION AND IS IGNORING SAME.')
 RETURN
END SUBROUTINE qdmm1s
