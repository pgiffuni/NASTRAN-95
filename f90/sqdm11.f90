SUBROUTINE sqdm11
     
!     PHASE I OF STRESS DATA RECOVERY FOR THE  QUADRILATERAL MEMBRANE
!     ELEMENT
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
 
!            MAT    - MATERIAL DATA ROUTINE
!            MESAGE - ERROR MESSAGE WRITER
!            GMMATS - SINGLE MATRIX MULTIPLY AND TRANSPOSE
!            TRANSS - SINGLE PRECISION TRANSFORMATION SUPPLIER
 
 
 REAL :: la,lb,lc,ld,ldd2,lbd1,lcd1,lcd2,magi,magj,magk
 DIMENSION       ecpt(26),ee(144)
 COMMON /system/ dummy(39),nbpw
 COMMON /condas/ consts(5)
 COMMON /sdr2x5/ necpt(1),ngrid(4),angle,matid1,t,fmu,  &
     dummy1,x1,y1,z1,dummy2,x2,y2,z2, dummy3,x3,y3,z3,dummy4,x4,y4,z4,dumb(75),  &
     ph1out(100),forvec(25)
 COMMON /sdr2x6/ e(9),ti(9),theta,tempar(150),a(24),g(9),b(96)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alphas(3),tsub0,gsube,  &
     sigten,sigcom,sigshe,g2x211,g2x212,g2x222
 EQUIVALENCE     (consts(4),degra),(ecpt(1),necpt(1))
 
 
!      ECPT LIST
!                                                      IN
!                                                      THIS
!        ECPT       DESCRIPTION                        ROUTINE    TYPE
!     *******    *********************************     ********  ******
!     ECPT( 1) = ELEMENT ID                            NECPT(1)  INTEGER
!     ECPT( 2)   GRID POINT A                          NGRID(1)  INTEGER
!     ECPT( 3)   GRID POINT B                          NGRID(2)  INTEGER
!     ECPT( 4)   GRID POINT C                          NGRID(3)  INTEGER
!     ECPT( 5)   GRID POINT D                          NGRID(4)  INTEGER
!     ECPT( 6) = THETA = ANGLE OF MATERIAL             ANGLE      REAL
!     ECPT( 7)   MATERIAL ID                           MATID     INTEGER
!     ECPT( 8) = THICKNESS                             T          REAL
!     ECPT( 9) = NON-STRUCTURAL MASS                   FMU        REAL
!     ECPT(10)   COORD. SYSTEM ID 1                    NECPT(10) INTEGER
!     ECPT(11) = X1                                     X1        REAL
!     ECPT(12) = Y1                                     Y1        REAL
!     ECPT(13) = Z1                                     Z1        REAL
!     ECPT(14)   COORD. SYSTEM ID 2                    NECPT(14) INTEGER
!     ECPT(15) = X2                                     X2        REAL
!     ECPT(16) = Y2                                     Y2        REAL
!     ECPT(17) = Z2                                     Z2        REAL
!     ECPT(18)   COORD. SYSTEM ID 3                    NECPT(18) INTEGER
!     ECPT(19) = X3                                     X3        REAL
!     ECPT(20) = Y3                                     Y3        REAL
!     ECPT(21) = Z3                                     Z3        REAL
!     ECPT(22)   COORD. SYSTEM ID 4                    NECPT(22) INTEGER
!     ECPT(23) = X4                                     X4        REAL
!     ECPT(24) = Y4                                     Y4        REAL
!     ECPT(25)   Z4                                     Z4        REAL
!     ECPT(26) = ELEMENT TEMPERATURE                    ECPT(26)  REAL
 
 
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
 
!     COMPUTE ELEMENTS OF THE E MATRIX (3X3)
 
 pk1 = y31*z42 - z31*y42
 pk2 = z31*x42 - x31*z42
 pk3 = x31*y42 - y31*x42
 magk= SQRT(pk1**2 + pk2**2 + pk3**2)
 IF (magk > 1.0E-06) GO TO 40
 CALL mesage (-30,32,ecpt(1))
 40 pk1 = pk1/magk
 pk2 = pk2/magk
 pk3 = pk3/magk
 
!     HH IS THE MEASURE OF NON-PLANARITY OF THE ELEMENT
 
 hh  = x21*pk1 + y21*pk2 + z21*pk3
 pi1 = x21 - hh*pk1
 pi2 = y21 - hh*pk2
 pi3 = z21 - hh*pk3
 magi= SQRT(pi1**2 + pi2**2 + pi3**2)
 IF (magi > 1.0E-06) GO TO 41
 CALL mesage (-30,31,ecpt(1))
 41 pi1 = pi1/magi
 pi2 = pi2/magi
 pi3 = pi3/magi
 hh  =-hh/2.
 
!     THIS SIGN CHANGE MADE BECAUSE SIGN OF H AS DEFINED ON
!     PAGE 4.87-105 OF PROGRAMMERS MANUAL IS WRONG
 
 pj1 = pk2*pi3 - pk3*pi2
 pj2 = pk3*pi1 - pk1*pi3
 pj3 = pk1*pi2 - pk2*pi1
 magj= SQRT(pj1**2 + pj2**2 + pj3**2)
 pj1 = pj1/magj
 pj2 = pj2/magj
 pj3 = pj3/magj
 e(1)= pi1
 e(2)= pj1
 e(3)= pk1
 e(4)= pi2
 e(5)= pj2
 e(6)= pk2
 e(7)= pi3
 e(8)= pj3
 e(9)= pk3
 
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
 
 x12 =-(x21*e(1) + y21*e(4) + z21*e(7))
 x13 =-(x31*e(1) + y31*e(4) + z31*e(7))
 x14 =-(x41*e(1) + y41*e(4) + z41*e(7))
 y3a = (x31*e(2) + y31*e(5) + z31*e(8))
 y4a = (x42*e(2) + y42*e(5) + z42*e(8))
 x24 = x14 - x12
 x23 = x13 - x12
 x34 = x14 - x13
 y34 = y3a - y4a
 
 
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
 
!                                                  *       *
!     COMPUTE THE INTERSECTION OF THE DIAGONALS(ETA  AND XI ) OF
!     THE ELEMENTS IN THE MEAN PLANE
 
 tol  = 1.0E-3*(-x12)
 IF (nbpw >= 60) tol  = 1.0E-5*(-x12)
 tol2 = 1.0E-3*x12*x12
 IF (nbpw >= 60) tol2 = 1.0E-5*x12*x12
 IF (ABS(x34+x12) > tol .OR. ABS(y34) > tol) GO TO 6
 etas =.5
 xis  =.5
 GO TO 16
 6 IF (ABS(x24) < tol .OR. ABS(x13) < tol) GO TO 7
 xstar = (y4a*x13*x12)/((y3a*x24)-(y4a*x13))
 ystar = (-y4a/x24)*(xstar+x12)
 GO TO 9
 7 IF (ABS(x13) > tol) GO TO 8
 xstar = -x13
 ystar = (-y4a/x24)*x12
 GO TO 9
 8 xstar = -x12
 ystar = (y3a/x13)*x12
 9 IF (ABS(x34+x12) < tol) GO TO 13
 c1 = y34*xstar - ystar*(x34+x12)
 a2 =-y4a*x23 + y3a*x14
 b2 =-y4a*x12 + c1
 IF (ABS(a2) <= tol2) GO TO 10
 temp2 = SQRT(b2**2-(4.*a2*x12*ystar))/(2.*a2)
 temp1 =-b2/(2.*a2)
 etas  = temp1 - temp2
 IF (etas <= 0. .OR. etas >= 1.) etas = temp1 + temp2
 GO TO 11
 10 etas = (-x12*ystar)/b2
 11 IF (ABS(y34) < tol) GO TO 12
 xis = (-c1 + ((y4a*x23) - (y3a*x14))*etas)/(y34*x12)
 GO TO 16
 12 xis = (xstar + (x14*etas))/((etas*(x34+x12)) - x12)
 GO TO 16
 13 a3 = -x14*y34
 b3 =  x12*y4a - y34*xstar
 IF (ABS(a3) <= tol2) GO TO 14
 temp2 = SQRT(b3**2 + (4.*a3*x12*ystar))/(2.*a3)
 temp1 = -b3/(2.*a3)
 etas  = temp1 - temp2
 IF (etas <= 0. .OR. etas >= 1.) etas = temp1 + temp2
 GO TO 15
 14 etas = (x12*ystar)/b3
 15 xis  = (ystar - (y4a*etas))/(y34*etas)
 
!     SET UP THE (12X12) TRANSFORMATION MATRIX B BETWEEN THE MEAN PLANE
!                        AND ACTUAL GRID POINTS
 
 16 DO  i = 2,92
   b(i)  = 0.
 END DO
 b(1)  = 1.
 b(10) = 1.
 b(17) =-hh/la
 b(18) =-hh/(ld*sth1) + ((hh*cth1)/(la*sth1))
 b(19) = hh/la
 b(20) = (hh*cth2)/(la*sth2)
 b(23) = (hh*cth42)/ldd2
 b(24) = (hh*sth42)/ldd2
 b(27) = 1.
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
 b(79) = 1.
 b(88) = 1.
 b(90) = hh/(ld*sth1)
 b(93) = (-hh*cth32)/lcd1
 b(94) = (-hh*sth32)/lcd1
 b(95) = hh*((-cth42/ldd2) + (sth41/lcd2))
 b(96) = hh*((-sth42/ldd2) - (cth41/lcd2))
 DO  i = 1,24
   a(i)  = 0.
 END DO
!                                                     T
!     COMPUTE TRANSFORMED MATRIX OF STIFFNESSES  G = P  * G * P
 
 theta  = angle*degra
 sinth  = SIN(theta)
 costh  = COS(theta)
 IF (ABS(sinth) < 1.0E-06) sinth = 0.0E0
 matid  = matid1
 inflag = 2
 eltemp = ecpt(26)
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
 
!     COMPUTE MATRIX A TO RELATE DISPLACEMENTS TO STRAINS
 
 aj    = (-y4a*x12) + (-y34*x12*xis) + etas*((-y4a*x23)+(y3a*x14))
 a(1)  = (-y4a + (y3a*etas) - (y34*xis))/aj
 a(3)  = ( y4a - (y4a*etas) + (y34*xis))/aj
 a(5)  = ( y4a*etas)/aj
 a(7)  = (-y3a*etas)/aj
 a(10) = (-x24 + (x23*etas) + (x34*xis))/aj
 a(12) = ( x14 - (x14*etas) - (x34*xis))/aj
 a(14) = ((x14*etas) - (x12*xis))/aj
 a(16) = (-x12 - (x23*etas) + (x12*xis))/aj
 a(17) = (-x24 + (x23*etas) + (x34*xis))/aj
 a(18) = (-y4a + (y3a*etas) - (y34*xis))/aj
 a(19) = ( x14 - (x14*etas) - (x34*xis))/aj
 a(20) = ( y4a - (y4a*etas) + (y34*xis))/aj
 a(21) = ((x14*etas) - (x12*xis))/aj
 a(22) = (y4a*etas)/aj
 a(23) = (-x12 - (x23*etas) + (x12*xis))/aj
 a(24) = (-y3a*etas)/aj
 
!                          T    T
!     COMPUTE S = G * A * B  * E
 
 CALL gmmats (b(1),12,8,1,ee(1),12,12,1,tempar(1))
 CALL gmmats (a(1),3,8,0,tempar(1),8,12,0,tempar(100))
 CALL gmmats (g(1),3,3,0,tempar(100),3,12,0,tempar(1))
 DO  l = 1,4
   DO  n = 2,5
     IF (necpt(n) /= ngrid(l)) CYCLE
     ka = 4*n + 2
     GO TO 20
   END DO
   CALL mesage (-30,34,ecpt(1))
   20 IF (necpt(ka) == 0) GO TO 21
   CALL transs (necpt(ka),ti)
   GO TO 23
   21 DO  ii = 1,9
     ti(ii) = 0.
   END DO
   ti(1) = 1.
   ti(5) = 1.
   ti(9) = 1.
   23 lcnt  = 3*(l-1)
   irowct= -12
   nn    = 0
   DO  jj = 1,3
     irowct= irowct + 12
     DO  kk = 1,3
       nn   = nn + 1
       ktot = kk + irowct + lcnt
       nn49 = nn + 49
       tempar(nn49) = tempar(ktot)
     END DO
   END DO
   CALL gmmats (tempar(50),3,3,0,ti,3,3,0,tempar(60))
   
!                                                          TH
!     MATRICES S  RELATE DISPLACEMENTS TO STRESSES AT THE I   GRIDPOINT
!               I
   
   DO  il = 1,9
     ktot = il + 9*l
     il59 = il + 59
     ph1out(ktot) = tempar(il59)
   END DO
 END DO
 CALL gmmats (g(1),3,3,0,alphas(1),3,1,0,ph1out(7))
 ph1out(1) = ecpt(1)
 ph1out(2) = ecpt(2)
 ph1out(3) = ecpt(3)
 ph1out(4) = ecpt(4)
 ph1out(5) = ecpt(5)
 ph1out(6) = tsub0
 RETURN
END SUBROUTINE sqdm11
