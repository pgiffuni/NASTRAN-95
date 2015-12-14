SUBROUTINE stord1
     
 
!*****
! THIS ROUTINE IS PHASE  I OF STRESS DATA RECOVERY FOR AN AXI-SYMMETRIC
! TOROIDAL THIN SHELL RING
!*****
 
 
 
!                        ECPT FOR THE TOROIDAL RING
 
!                                                       TYPE
! ECPT( 1) ELEMENT IDENTIFICATION                         I
! ECPT( 2) SCALAR INDEX NO. FOR GRID POINT A              I
! ECPT( 3) SCALAR INDEX NO. FOR GRID POINT B              I
! ECPT( 4) ANGLE OF CURVATURE AT GRID POINT A             R
! ECPT( 5) ANGLE OF CURVATURE AT GRID POINT B(NOT USED)   R
! ECPT( 6) MATERIAL ORIENTATION (NOT USED)                R
! ECPT( 7) MATERIAL IDENTIFICATION                        I
! ECPT( 8) MEMBRANE THICKNESS                             R
! ECPT( 9) FLEXURE THICKNESS                              R
! ECPT(10) COOR. SYS. ID. FOR GRID POINT A                I
! ECPT(11) X-COOR. OF GRID POINT A (IN BASIC COOR.)       R
! ECPT(12) Y-COOR. OF GRID POINT A (IN BASIC COOR.)       R
! ECPT(13) Z-COOR. OF GRID POINT A (IN BASIC COOR.)       R
! ECPT(14) COOR. SYS. ID. FOR GRID POINT B                I
! ECPT(15) X-COOR. OF GRID POINT B (IN BASIC COOR.)       R
! ECPT(16) Y-COOR. OF GRID POINT B (IN BASIC COOR.)       R
! ECPT(17) Z-COOR. OF GRID POINT B (IN BASIC COOR.)       R
! ECPT(18) EL. TEMPERATURE FOR MATERIAL PROPERTIES        R
 
 
 DIMENSION          iecpt(18)
 DIMENSION          gambqf(72),    gambqm(48)
 DIMENSION          ee(4), gambq(144), gamrs(144)
 DIMENSION          aki(36),  delint(66)
 DIMENSION                    ics(2)
 DIMENSION          gambl(144)
 
 COMMON /condas/ consts(5)
 COMMON   /sdr2x5/ ecpt(18)  &
     ,                  dum5(82) ,                  idel,     igp(2),   tz  &
     ,                  sel(180), ts(30),   ak(144)
 COMMON   /matin/ matidc        ,matflg  &
     ,                  eltemp        ,stress  &
     ,                  sinth         ,costh
 COMMON   /matout/ e(3)          ,anu(3)  &
     ,                  rho           ,g(3)  &
     ,                  alf(3)        ,tzero,    gsube
 COMMON   /sdr2x6/ d(180),   r(2),     z(2),     alph(2)
 
 EQUIVALENCE ( consts(2) , twopi  )
 EQUIVALENCE ( consts(4) , degra  )
 EQUIVALENCE        (iecpt(1) , ecpt(1))
 EQUIVALENCE        (a1, alph(1)), (a2, alph(2))
 EQUIVALENCE        (r1, r(1)),    (r2, r(2))
 EQUIVALENCE        (z1, z(1)),    (z2, z(2))
 EQUIVALENCE        (gambqf(1), gambq(1))
 EQUIVALENCE        (gambqm(1), gambq(73))
 EQUIVALENCE        (delint(1), gambq(1))
 EQUIVALENCE        (gamrs(1),  gambq(1))
 EQUIVALENCE        (aki(1),    gambq(1))
 EQUIVALENCE        (gambl(1), gambq(1))
 
! ----------------------------------------------------------------------
 
! STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel   = iecpt(1)
 igp(1) = iecpt(2)
 igp(2) =  iecpt(3)
 matid  =  iecpt(7)
 ics(1) =  iecpt(10)
 ics(2) =  iecpt(14)
 alph(1)=  ecpt(4)
 alph(2)=  ecpt(5)
 tm     =  ecpt(8)
 tf     =  ecpt(9)
 r(1)   =  ecpt(11)
 d(1)   =  ecpt(12)
 z(1)   =  ecpt(13)
 r(2)   =  ecpt(15)
 d(2)   =  ecpt(16)
 z(2)   =  ecpt(17)
 tempe  =  ecpt(18)
 
 
! TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 DO  i = 1,2
   IF (r(i) < 0.0E0) CALL mesage(-30,37,idel)
   IF (d(i) /= 0.0E0) CALL mesage(-30,37,idel)
 END DO
 
 
! DETERMINE IF ELEMENT IS A TOROIDAL, CONICAL OR CYLINDRICAL RING
 
 itord = 0
 IF (ABS(a1-a2) <= .000001) itord = 1
 IF (itord == 1  .AND.  ABS(a1-90.0E0) <= .00001) itord = -1
 
 
! COMPUTE THE ELEMENT COORDINATES
 
 a1 = a1 * degra
 a2 = a2 * degra
 phib = a2 - a1
 sina1 =  SIN(a1)
 cosa1 =  COS(a1)
 sina2 =  SIN(a2)
 cosa2 =  COS(a2)
 
 IF (itord /= 0) GO TO 300
 
! FOR THE TOROIDAL RING
 
 rp =  SQRT( (r2-r1)**2 + (z2-z1)**2 ) / (2.0E0 *  SIN(phib/2.0E0))
 s = phib * rp
 GO TO 350
 
! FOR THE CONICAL OR CYLINDRICAL RING
 
 300 CONTINUE
 rp = 0.0D0
 s  =  SQRT( (r2-r1)**2 + (z2-z1)**2 )
 
 350 CONTINUE
 
 
! COMPUTE THE BASIC AND REQUIRED INTEGRALS
 
 
! SET UP ARRAY OF CONSTANTS FOR ROMBER INTEGRATION ROUTINE
 
 d(21) = 0.0E0
 d(22) = rp
 d(23) = r1
 d(24) = cosa1
 d(25) = sina1
 
! COMPUTE CONSTANTS NEEDED FOR INTEGRAL CALCULATIONS
 
 d(30) = r1 - rp * sina1
 d(31) = rp * cosa1
 d(32) = rp * sina1
 d(33) = cosa1 ** 2
 d(34) = sina1 * cosa1
 d(35) = sina1 ** 2
 d(36) = 0.5 - d(35)
 
! START LOOP  FOR CALCULATIONS OF INTEGRALS
 
 DO  jp1 = 1,11
   j = jp1 - 1
   k = (j * 6) + 1
   djp1 = jp1
   
! TEST FOR ELEMENT SHAPE
   
   IF (itord < 0) THEN
     GO TO   470
   ELSE IF (itord == 0) THEN
     GO TO   400
   ELSE
     GO TO   430
   END IF
   
! THE TOROIDAL RING BASIC INTEGRALS WILL BE COMPUTED IN
! LOCATIONS D(1),...,D(6)
   
   400 CONTINUE
   d(20) = (rp ** jp1)
   
! COMPUTE I(J,1)
   
   d(1) = d(20) * (phib ** jp1) / djp1
   
! COMPUTE I(J,2)
   
   d(2) = (phib ** (jp1+1)) / (djp1 + 1.0E0)
   d(10) = 1.0E0
   DO  i = 1,20
     ip = jp1 + 2 * i + 1
     d(11) = 2 * i + 1
     d(10) = d(10) * d(11) * (d(11)-1.0E0)
     d(12) = (-1.0E0)** i  * phib ** ip / ((djp1 + d(11)) * d(10))
     d(13) =  ABS( d(12) / d(2) )
     d(2) = d(2) + d(12)
     IF (d(13) <= 1.0E-10) GO TO 415
   END DO
   CALL mesage(-30,26,idel)
   415 CONTINUE
   d(2) = d(20) * d(2)
   
! COMPUTE I(J,3)
   
   d(3) = (phib ** jp1) / djp1
   d(10) = 1.0E0
   DO  i = 1,20
     ip = jp1 + 2 * i
     d(11) = 2 * i
     d(10) = d(10) * d(11) * (d(11) - 1.0E0)
     d(12) = (-1.0E0)** i  * phib ** ip / ((djp1 + d(11)) * d(10))
     d(13) =  ABS( d(12) / d(3) )
     d(3) = d(3) + d(12)
     IF (d(13) <= 1.0E-10) GO TO 425
   END DO
   CALL mesage(-30,26,idel)
   425 CONTINUE
   d(3) = d(20) * d(3)
   d(26) = djp1
   
! COMPUTE I(J,4)
   
   CALL romber (phib, d(10), ip, d(4), 1, d(21) )
   IF (ip >= 15) CALL mesage (30,26,idel)
   d(4) = d(20) * d(4)
   
! COMPUTE I(J,5)
   
   CALL romber (phib, d(10), ip, d(5), 2, d(21) )
   IF (ip >= 15) CALL mesage (30,26,idel)
   d(5) = d(20) * d(5)
   
! COMPUTE I(J,6)
   
   CALL romber (phib, d(10), ip, d(6), 3, d(21) )
   IF (ip >= 15) CALL mesage (30,26,idel)
   d(6) = d(20) * d(6)
   
! THE TOROIDAL RING REQUIRED INTEGRALS
   
   delint(k  ) = d(30) * d(1) + d(31) * d(2) + d(32) * d(3)
   delint(k+1) = cosa1 * d(2) + sina1 * d(3)
   delint(k+2) = d(33) * d(4) + d(34) * d(5) + d(35) * d(6)
   delint(k+3) = cosa1 * d(3) - sina1 * d(2)
   delint(k+4) = d(34) * (d(6)-d(4))  + d(36) * d(5)
   delint(k+5) = d(33) * d(6) - d(34) * d(5) + d(35) * d(4)
   GO TO 490
   
! THE CONICAL RING BASIC INTEGRALS WILL BE COMPUTED IN
! LOCATIONS D(1) AND D(2)
   
   430 CONTINUE
   
! COMPUTE I(J,1)
   
   d(1) = (s ** jp1) / djp1
   
   IF (j - 1 < 0) THEN
     GO TO   435
   ELSE IF (j - 1 == 0) THEN
     GO TO   440
   ELSE
     GO TO   445
   END IF
   
! COMPUTE I(0,2)
   
   435 CONTINUE
   d(2) = ALOG( (r1 + s*cosa1) / r1 ) / cosa1
   GO TO 460
   
! COMPUTE I(1,2)
   
   440 CONTINUE
   d(2) = (s - (r1/cosa1) * ALOG( (r1 + s*cosa1) / r1 )) / cosa1
   GO TO 460
   
! COMPUTE I(J,2) WHERE J .GT. 1
   
   445 CONTINUE
   d(2) = 1.0E0 / djp1
   d(10) =-s * cosa1 / r1
   DO  i = 1,1000
     d(11) = jp1 + i
     d(12) = (d(10) ** i) / d(11)
     d(2) = d(2) + d(12)
     IF (d(12) < 1.0E-4 ) GO TO 455
   END DO
   CALL mesage(-30,26,idel)
   455 CONTINUE
   d(2) = ( (s ** jp1) / r1 ) * d(2)
   460 CONTINUE
   
! THE CONICAL RING REQUIRED INTEGRALS
   
   delint(k  ) = r1 * d(1) + cosa1 * ((s**(jp1+1)) / (djp1+1.0E0))
   delint(k+1) = sina1 * d(1)
   delint(k+2) = d(35) * d(2)
   delint(k+3) = cosa1 * d(1)
   delint(k+4) = d(34) * d(2)
   delint(k+5) = d(33) * d(2)
   GO TO 490
   
! THE CYLINDRICAL RING BASIC INTEGRALS WILL BE COMPUTED IN
! LOCATIONS D(1) AND D(2)
   
   470 CONTINUE
   
! COMPUTE I(J,1)
   
   d(1) = (s ** jp1) / djp1
   
! COMPUTE I(J,2)
   
   d(2) = d(1) / r1
   
! THE CYLINDRICAL RING REQUIRED INTEGRALS
   
   delint(k  ) = r1 * d(1) + cosa1 * ((s**(jp1+1)) / (djp1+1.0E0))
   delint(k+1) = sina1 * d(1)
   delint(k+2) = d(35) * d(2)
   delint(k+3) = 0.0E0
   delint(k+4) = 0.0E0
   delint(k+5) = 0.0E0
   
   490 CONTINUE
 END DO
 
 
! LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3 TABLE
 
 matidc = matid
 matflg = 7
 eltemp = tempe
 CALL mat(idel)
 
 
! SET MATERIAL PROPERTIES IN LOCAL VARIABLES
 
 ep = e(1)
 et = e(2)
 vpt= anu(1)
 tz = tzero
 vtp= vpt * et / ep
 del = 1.0E0 - vpt * vtp
 
 
! GENERATE THE ELASTIC CONSTANTS MATRIX(2X2)
 
 ee(1) = ep / del
 ee(2) = et * vpt / del
 ee(3) = ee(2)
 ee(4) = et / del
 
 
! FORM THE STIFFNESS MATRIX IN FIELD COORDINATES
 
! COMPUTE CONSTANTS NEEDED IN DMATRX SUBROUTINE
 
 d(1) = ep / et
 d(7) = 0.0E0
 IF (itord == 0) d(7) = 1.0E0 / rp
 d(2) = d(1) * d(7)
 d(3) = d(2) * d(7)
 d(4) = vpt * d(7)
 d(5) =(ep * tm / (d(1) - vpt**2)) * twopi
 d(6) =(ep * (tf**3) / (12.0E0 * (d(1) - vpt**2))) * twopi
 
! CALL THE AMATRIX SUBROUTINE TO COMPUTE THE STIFFNESS MATRIX (10X10)
 
! NOTE THE DOUBLE SUBSCRIPTING USED IN AMATRIX SUBROUTINE IS
! COMPATIBLE WITH THE CALLING PROGRAM. THE DELINT ARRAY OF INTEGRALS
! IS A (11X6) SINGLY SUBSCRIPTED ARRAY (STORED ROWWISE) IN THE CALLING
! PROGRAM AND IT IS A (6X11) DOUBLY SUBSCRIPTED ARRAY (STORED
! COLUMNWISE) IN AMATRX ROUTINE.
 
 
 CALL amatrx(ak(1), vpt, d(1), d(2), d(3), d(4), d(5), d(6)  &
     ,           delint(1) )
 
 
! FORM THE STRESS MATRIX IN FIELD COORDINATES
 
! COMPUTE THE CONSTANTS NEEDED IN THE SCRLM SUBROUTINE
 
 d(1) = 0.0E0
 IF (itord == 0) d(1) = 1.0E0 / rp
 d(2) = 0.0E0
 d(3) = s / 2.0E0
 d(4) = s
 
! CALL THE SCRLM SUBROUTINE TO COMPUTE THE STRESS MATRIX TRANSPOSED
 
! NOTE THE DOUBLE SUBSCRIPTING USED IN THE SCRLM SUBROUTINE IS
! COMPATIBLE WITH THE CALLING PROGRAM. THE SEL ARRAY WILL RETURN WITH
! THE STRESS MATRIX TRANSPOSED (10X15, STORED ROWWISE) BUT IN THE SCRLM
! SUBROUTINE THE STRESS MATRIX IS COMPUTED AS A DOUBLY SUBSCRIPTED
! 15X10 ARRAY (STORED COLUMNWISE).
 
 
 CALL scrlm (sel(1), d(2), ee(1), tm, 0.0E0, rp, a1, r1, d(1), tf)
 
 
! FORM THE TRANSFORMATION MATRIX(10X12) FROM FIELD COORDINATES TO GRID
! POINT DEGREES OF FREEDOM
 
 DO  i = 1,72
   gambqf(i) = 0.0E0
 END DO
 d(1) = s
 d(2) = s ** 2
 d(3) = s ** 3
 d(4) = s ** 4
 d(5) = s ** 5
 gambqf( 3) = 1.0E0
 gambqf(16) = 1.0E0
 gambqf(30) = 0.5E0
 gambqf(39) = -10.0E0 / d(3)
 gambqf(40) = - 6.0E0 / d(2)
 gambqf(42) = - 1.5E0 / d(1)
 gambqf(45) = -gambqf(39)
 gambqf(46) = - 4.0E0 / d(2)
 gambqf(48) =   0.5E0 / d(1)
 gambqf(51) =  15.0E0 / d(4)
 gambqf(52) =   8.0E0 / d(3)
 gambqf(54) =   1.5E0 / d(2)
 gambqf(57) = -gambqf(51)
 gambqf(58) =   7.0E0 / d(3)
 gambqf(60) = - 1.0E0 / d(2)
 gambqf(63) = - 6.0E0 / d(5)
 gambqf(64) = - 3.0E0 / d(4)
 gambqf(66) = - 0.5E0 / d(3)
 gambqf(69) = -gambqf(63)
 gambqf(70) =  gambqf(64)
 gambqf(72) = -gambqf(66)
 DO  i = 1,48
   gambqm(i) = 0.0E0
 END DO
 gambqm( 1) = 1.0E0
 gambqm(17) = 1.0E0
 gambqm(25) = - 3.0E0 / d(2)
 gambqm(29) = - 2.0E0 / d(1)
 gambqm(31) = -gambqm(25)
 gambqm(35) = - 1.0E0 / d(1)
 gambqm(37) =   2.0E0 / d(3)
 gambqm(41) =   1.0E0 / d(2)
 gambqm(43) = -gambqm(37)
 gambqm(47) =  gambqm(41)
 
 
! TRANSFORM THE STIFFNESS MATRIX TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmats(gambq(1), 10, 12, 1, ak(1), 10, 10, 0, d(1) )
 CALL gmmats(d(1), 12, 10,  0, gambq(1), 10, 12, 0, ak(1) )
 
 
! RE-ARRANGE THE TRANSFORMATION MATRIX (GAMBQ) SUCH THAT THE MEMBRANE
! AND FLEXURE TERMS ARE REVERSED
 
 DO  i = 1,72
   d(i) = gambqf(i)
 END DO
 DO  i = 1,48
   gambq(i) = gambqm(i)
 END DO
 DO  i = 1,72
   gambq(i+48) = d(i)
 END DO
 
 
! TRANSFORM THE STRESS MATRIX TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmats (sel(1), 10, 15, 1, gambq(1), 10, 12, 0, d(1) )
 
 
! FORM THE TRANSFORMATION MATRIX (12X12) FROM ELEMENT TO BASIC
! COORDINATES
 
 DO  i = 1,144
   gamrs(i) = 0.0E0
 END DO
 gamrs( 1) =  cosa1
 gamrs( 3) = -sina1
 gamrs(25) =  sina1
 gamrs(27) =  cosa1
 gamrs(40) = -1.0E0
 gamrs(53) =  1.0E0
 gamrs(66) =  1.0E0
 gamrs(79) =  cosa2
 gamrs(81) = -sina2
 gamrs(103)=  sina2
 gamrs(105)=  cosa2
 gamrs(118)= -1.0E0
 gamrs(131)=  1.0E0
 gamrs(144)=  1.0E0
 
 
! TRANSFORM THE STRESS MATRIX FROM ELEMENT TO BASIC COORDINATES
 
 CALL gmmats (  d(1), 15, 12, 0, gamrs(1), 12, 12, 0, sel(1) )
 
 
! TRANSFORM THE STIFFNESS MATRIX FROM ELEMENT TO BASIC COORDINATES
 
 CALL gmmats(gamrs(1), 12, 12, 1, ak(1), 12, 12, 0, d(1) )
 CALL gmmats(d(1), 12, 12,  0, gamrs(1), 12, 12, 0, ak(1) )
 
 
! LOCATE THE TRANSFORMATION MATRICES FROM BASIC TO LOCAL COORDINATES
! FOR THE TWO GRID POINTS AND EXPAND TO (6X6)
 
 DO  i = 1,144
   gambl(i) = 0.0E0
 END DO
 DO  i = 1,2
   CALL transs (ics(i) , d(1))
   k = 78 * (i - 1)
   DO  j = 1,3
     kk = k + 12* (j-1) + 1
     kl = 3 * (j-1) + 1
     kj = k + 12* (j+2) + j + 3
     gambl(kk  ) = d(kl  )
     gambl(kk+1) = d(kl+1)
     gambl(kk+2) = d(kl+2)
     gambl(kj) = 1.0E0
   END DO
 END DO
 
 
 
! TRANSFORM THE STIFFNESS MATRIX FROM BASIC TO LOCAL COORDINATES
 
 CALL gmmats (gambl(1), 12, 12, 1, ak(1), 12, 12, 0, d(1) )
 CALL gmmats (d(1), 12, 12,  0, gambl(1), 12, 12, 0, ak(1) )
 
 
! TRANSFORM THE STRESS MATRIX FROM BASIC TO LOCAL COORDINATES
 
 CALL gmmats (sel(1), 15, 12, 0, gambl(1), 12, 12, 0, d(1) )
 
 DO  i = 1,180
   sel(i) = d(i)
 END DO
 
 
! FORM THE THERMAL STRESS VECTOR (30X1)
 
! THE MEMBRANE TEMPERATURE TERMS WILL BE STORED IN TS(1),...,TS(15) AND
! THE FLEXURE GRADIENT TEMP. TERMS WILL BE STORED IN TS(16),...,TS(30)
 
 
! COMPUTE CONSTANTS NEEDED IN THE THERMAL STRESS CALCULATIONS
 
 d(1) = 0.0E0
 d(2) = s / 2.0E0
 d(3) = s
 d(4) = ee(1) * alf(1) + ee(2) * alf(2)
 d(5) = ee(3) * alf(1) + ee(4) * alf(2)
 d(6) = (ee(1)-ee(2)) * alf(1) + (ee(3)-ee(4)) * alf(2)
 d(7) = tf ** 3 / 12.0E0
 d(8) = tm / s
 d(9) = d(7) / s
 
! START THE LOOP TO FORM THE THERMAL STRESS VECTORS AT EACH OF THE
! THREE STRESS POINTS
 
 DO  i = 1,3
   CALL solve1(a1, r1, rp, d(i), d(12), d(13), d(14), 0.0E0)
   k = 5 * (i - 1)
   kk = k + 15
   ts(k +1) = tm * d(4)
   ts(k +2) = tm * d(5)
   ts(k +3) = d(7) * d(4)
   ts(k +4) =-d(7) * d(5)
   ts(k +5) = d(7) * d(12) * d(6)
   ts(kk+1) = d(8) * d(i)  * d(4)
   ts(kk+2) = d(8) * d(i)  * d(5)
   ts(kk+3) = d(9) * d(i)  * d(4)
   ts(kk+4) =-d(9) * d(i)  * d(5)
   ts(kk+5) = d(9) * (d(4) + d(i) * d(12) * d(6))
 END DO
 
 
 RETURN
END SUBROUTINE stord1
