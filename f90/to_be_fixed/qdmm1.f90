SUBROUTINE qdmm1 (tbar,pg)
     
!     QUADRILATERAL MEMBRANE ELEMENT
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
!           MAT    - MATERIAL DATA ROUTINE
!           MESAGE - ERROR MESSAGE WRITER
!           BASGLB - TRANSFER COORDINATES FROM BASIC TO GLOBAL
!           GMMATS - SINGLE PRECISION MATRIX MULTIPLY AND TRANSPOSE
!           TRANSS - SINGLE PRECISION TRANSFORMATION SUPPLIER
 
!     ECPT LIST
!                                                    IN THIS
!      ECPT       DESCRIPTION                        ROUTINE   TYPE
!     ========   =================================   ========  =======
!     ECPT( 1) = ELEMENT ID                          NECPT(1)  INTEGER
!     ECPT( 2)   GRID POINT A                        NGRID(1)  INTEGER
!     ECPT( 3)   GRID POINT B                        NGRID(2)  INTEGER
!     ECPT( 4)   GRID POINT C                        NGRID(3)  INTEGER
!     ECPT( 5)   GRID POINT D                        NGRID(4)  INTEGER
!     ECPT( 6) = THETA = ANGLE OF MATERIAL           ANGLE      REAL
!     ECPT( 7)   MATERIAL ID                         MATID1    INTEGER
!     ECPT( 8) = THICKNESS                           T          REAL
!     ECPT( 9) = NON-STRUCTURAL MASS                 FMU        REAL
!     ECPT(10)   COORD. SYSTEM ID 1                  NECPT(10) INTEGER
!     ECPT(11) = X1                                   X1        REAL
!     ECPT(12) = Y1                                   Y1        REAL
!     ECPT(13) = Z1                                   Z1        REAL
!     ECPT(14)   COORD. SYSTEM ID 2                  NECPT(14) INTEGER
!     ECPT(15) = X2                                   X2        REAL
!     ECPT(16) = Y2                                   Y2        REAL
!     ECPT(17) = Z2                                   Z2        REAL
!     ECPT(18)   COORD. SYSTEM ID 3                  NECPT(18) INTEGER
!     ECPT(19) = X3                                   X3        REAL
!     ECPT(20) = Y3                                   Y3        REAL
!     ECPT(21) = Z3                                   Z3        REAL
!     ECPT(22)   COORD. SYSTEM ID 4                  NECPT(22) INTEGER
!     ECPT(23) = X4                                   X4        REAL
!     ECPT(24) = Y4                                   Y4        REAL
!     ECPT(25)   Z4                                   Z4        REAL
 
 
 REAL, INTENT(IN)                         :: tbar
 REAL, INTENT(OUT)                        :: pg(1)
 REAL :: la,lb,lc,ld,ldd2,lbd1,lcd1,lcd2,magi,magj,magk
 DIMENSION       ecpt(26), g(9),e(9)
 COMMON /condas/ consts(5)
 COMMON /trimex/ necpt(1),ngrid(4),angle,matid1,t,fmu,  &
     dummy1,x1,y1,z1,dummy2,x2,y2,z2,dummy3,x3,y3,z3, dummy4,x4,y4,z4
 COMMON /ssgwrk/ ee(144),b(96),tempar(24),c(24),ti(9)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alphas(3),  &
     t_sub_0,g sub e,sigten,sigcom,sigshe, g2x211,g2x212,g2x222
 EQUIVALENCE     (consts(4),degra),(ecpt(1),necpt(1))
 
!     SET UP THE E MATRIX WHICH IS (12X12) FOR THE QUAD-MEMBRANE PROJECT
!                         ONTO THE MEAN PLANE
 
 DO  i = 1,144
   ee(i) = 0.
 END DO
 
!     E(1), E(4), E(7) WILL BE THE I-VECTOR
!     E(2), E(5), E(8) WILL BE THE J-VECTOR
!     E(3), E(6), E(9) WILL BE THE K-VECTOR
 
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
 
!     COMPUTE THE ELEMENTS OF THE (3X3) E MATRIX
 
 pk1  = y31*z42 - z31*y42
 pk2  = z31*x42 - x31*z42
 pk3  = x31*y42 - y31*x42
 magk = SQRT(pk1**2 + pk2**2 + pk3**2)
 IF (magk > 1.0E-06) GO TO 40
 CALL mesage (-30,32,ecpt(1))
 40 pk1  = pk1/magk
 pk2  = pk2/magk
 pk3  = pk3/magk
 
!     HH IS THE MEASURE OF NON-PLANARITY OF THE ELEMENT
 
 hh   = x21*pk1 + y21*pk2 + z21*pk3
 pi1  = x21 - hh*pk1
 pi2  = y21 - hh*pk2
 pi3  = z21 - hh*pk3
 magi = SQRT(pi1**2 + pi2**2 + pi3**2)
 IF (magi > 1.0E-06) GO TO 41
 CALL mesage (-30,31,ecpt(1))
 41 pi1  = pi1/magi
 pi2  = pi2/magi
 pi3  = pi3/magi
 hh   =-hh/2.
 
!     THIS SIGN CHANGE MADE BECAUSE SIGN OF H AS DEFINED ON
!     PAGE 4.87-105 OF PROGRAMMERS MANUAL IS WRONG
 
 pj1  = pk2*pi3 - pk3*pi2
 pj2  = pk3*pi1 - pk1*pi3
 pj3  = pk1*pi2 - pk2*pi1
 magj = SQRT(pj1**2 + pj2**2 + pj3**2)
 pj1  = pj1/magj
 pj2  = pj2/magj
 pj3  = pj3/magj
 
!     INSERT ELEMENTS INTO THE (3X3) E MATRIX
 
 e(1) = pi1
 e(2) = pj1
 e(3) = pk1
 e(4) = pi2
 e(5) = pj2
 e(6) = pk2
 e(7) = pi3
 e(8) = pj3
 e(9) = pk3
 
!     STORE FOUR (3X3) E MATRICES INTO (12X12) E MATRIX
 
 llct = -39
 DO  iict = 1,12,3
   llct = llct + 39
   nnct = 0
   mmct =-12
   DO  jjct = 1,3
     mmct = mmct + 12
     DO  kkct = 1,3
       nnct = nnct + 1
       ktot = kkct + llct + mmct
       ee(ktot) = e(nnct)
     END DO
   END DO
 END DO
 
!     COMPUTE DIFFERENCES OF COORDINATES OF GRID POINTS IN THE MEAN PLAN
 
 x12 = -x21*e(1) - y21*e(4) - z21*e(7)
 x13 = -x31*e(1) - y31*e(4) - z31*e(7)
 x14 = -x41*e(1) - y41*e(4) - z41*e(7)
 y3a =  x31*e(2) + y31*e(5) + z31*e(8)
 y4a =  x42*e(2) + y42*e(5) + z42*e(8)
 x24 =  x14 - x12
 x23 =  x13 - x12
 x34 =  x14 - x13
 y34 =  y3a - y4a
 
!     COMPUTE LENGTHS OF SIDES OF ELEMENT IN THE MEAN PLANE
 
 la = ABS(x12)
 lb = SQRT(x23**2 + y3a**2)
 lc = SQRT(x34**2 + y34**2)
 ld = SQRT(x14**2 + y4a**2)
 
!     COMPUTE THE CHARACTERISTIC ANGLES OF ELEMENT IN THE MEAN PLANE
 
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
 
!     SET UP THE (12X12) TRANSFORMATION MATRIX B BETWEEN THE MEAN PLANE
!                        AND ACTUAL GRID POINTS
 
 DO  i = 1,96
   b( i) = 0.
 END DO
 b( 1) = 1.
 b(10) = 1.
 b(17) =-hh/la
 b(18) =-hh/(ld*sth1)+((hh*cth1)/(la*sth1))
 b(19) = hh/la
 b(20) = (hh*cth2)/(la*sth2)
 b(23) = (hh*cth42)/ldd2
 b(24) = (hh*sth42)/ldd2
 b(27) = 1.
 b(36) = 1.
 b(41) =-b(17)
 b(42) =-(hh*cth1)/(la*sth1)
 b(43) = b(17)
 b(44) = ((-hh*cth2)/(la*sth2))+(hh/(lb*sth2))
 b(45) =-(hh*sth31)/lbd1
 b(46) =-(hh*cth31)/lbd1
 b(53) = 1.
 b(62) = 1.
 b(68) =-hh/(lb*sth2)
 b(69) = hh*((sth31/lbd1)+(cth32/lcd1))
 b(70) = hh*((cth31/lbd1)+(sth32/lcd1))
 b(71) =-(hh*sth41)/lcd2
 b(72) = (hh*cth41)/lcd2
 b(79) = 1.
 b(88) = 1.
 b(90) = hh/(ld*sth1)
 b(93) =-(hh*cth32)/lcd1
 b(94) =-(hh*sth32)/lcd1
 b(95) = hh*((-cth42/ldd2)+(sth41/lcd2))
 b(96) = hh*((-sth42/ldd2)-(cth41/lcd2))
 h     = ecpt( 8)
 eltemp= ecpt(26)
 
!     SET UP (3X8) C MATRIX (SEE FMMS)
 
 c( 1) =-(h*y4a)/2.
 c( 2) = 0.
 c( 3) =-(h*x24)/2.
 c( 4) = 0.
 c( 5) =-(h*x24)/2.
 c( 6) =-(h*y4a)/2.
 c( 7) = (h*y3a)/2.
 c( 8) = 0.
 c( 9) = (h*x13)/2.
 c(10) = 0.
 c(11) = (h*x13)/2.
 c(12) = (h*y3a)/2.
 c(13) = (h*y4a)/2.
 c(14) = 0.
 c(15) = (h*x24)/2.
 c(16) = 0.
 c(17) = (h*x24)/2.
 c(18) = (h*y4a)/2.
 c(19) =-(h*y3a)/2.
 c(20) = 0.
 c(21) =-(h*x13)/2.
 c(22) = 0.
 c(23) =-(h*x13)/2.
 c(24) =-(h*y3a)/2.
 theta = angle*degra
 sinth = SIN(theta)
 costh = COS(theta)
 IF (ABS(sinth) < 1.0E-06) sinth = 0.0
 matid  = matid1
 inflag = 2
!                                                     T
!     COMPUTE TRANSFORMED MATRIX OF STIFFNESSES  G = P  * G * P
 
 CALL mat (ecpt(1))
 
!     STORE INTO G MATRIX
 
 g(1) = g11
 g(2) = g12
 g(3) = g13
 g(4) = g12
 g(5) = g22
 g(6) = g23
 g(7) = g13
 g(8) = g23
 g(9) = g33
 
!                   T                            -
!     COMPUTE PG = T  * E * B * C * G * ALPHA * (T - T )
!                                                     0
 
 temp = tbar - tsub0
 tempar(1) = alphas(1)*temp
 tempar(2) = alphas(2)*temp
 tempar(3) = alphas(3)*temp
 CALL gmmats (g(1),3,3,0,tempar(1),3,1,0,tempar(13))
 CALL gmmats (c(1),8,3,0,tempar(13),3,1,0,tempar(1))
 CALL gmmats (b(1),12,8,0,tempar(1),8,1,0,tempar(13))
 CALL gmmats (ee(1),12,12,0,tempar(13),12,1,0,tempar(1))
 DO  i = 1,4
   
!     T-SUB-I WILL BE USED BELOW ONLY IF THE PIVOT COORDINATE SYSTEM ID
!     IS NOT ZERO, OTHERWISE IT IS ASSUMED TO BE THE IDENTITY MATRIX.
   
   ka = 4*i + 6
   
!     DO WE NEED TRANSFORMATION TI
   
   isw = 0
   jj  = 3*i - 2
   IF (necpt(ka) == 0) GO TO 9
   isw = 1
   CALL basglb (tempar(jj),tempar(20),necpt(ka+1),necpt(ka))
   
!     COMPUTE PG VECTOR
   
   9 DO  k = 1,3
     jjk = jj + k - 1
     k19 = k + 19
     IF (isw == 0) tempar(k19) = tempar(jjk)
     i1 = i + 1
     l  = necpt(i1) + k - 1
     pg(l) = pg(l) + tempar(k19)
   END DO
 END DO
 RETURN
END SUBROUTINE qdmm1
