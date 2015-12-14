SUBROUTINE ttordr (ti, pg)
     
 
!*****
! THIS ROUTINE COMPUTES THE THERMAL LOAD FOR AN AXI-SYMMETRIC
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
 
 
 
 REAL, INTENT(IN)                         :: ti(2)
 REAL, INTENT(OUT)                        :: pg(1)
 
 DIMENSION          iecpt(18)
 DIMENSION          gambqf(72),    gambqm(48)
 DIMENSION          ee(4), gambq(144), gamrs(144)
 DIMENSION          aki(36),  delint(42)
 DIMENSION          igp(2),   ics(2)
 DIMENSION          gambl(144)
 DIMENSION          d( 36),   r(2),     z(2),     alph(2)
 DIMENSION          fme(40),  ffe(40),  tl(12)
 
 COMMON /condas/ consts(5)
 COMMON   /trimex/ ecpt(18)
 COMMON   /matin/ matidc        ,matflg  &
     ,                  eltemp        ,stress  &
     ,                  sinth         ,costh
 COMMON   /matout/ e(3)          ,anu(3)  &
     ,                  rho           ,g(3)  &
     ,                  alf(3)        ,tzero,    gsube
 
 EQUIVALENCE ( consts(2) , twopi  )
 EQUIVALENCE ( consts(4) , degra  )
 EQUIVALENCE        (iecpt(1) , ecpt(1))
 EQUIVALENCE        (a1, alph(1)), (a2, alph(2))
 EQUIVALENCE        (r1, r(1)),    (r2, r(2))
 EQUIVALENCE        (z1, z(1)),    (z2, z(2))
 EQUIVALENCE        (gambqm(1), gambq(1))
 EQUIVALENCE        (gambqf(1), gambq(49))
 EQUIVALENCE        (delint(1), gambq(1))
 EQUIVALENCE        (fme(1),    gambq(43))
 EQUIVALENCE        (ffe(1),    gambq(83))
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
 
 DO  jp1 = 1,7
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
 
 
! CALL THE FCURL   SUBROUTINE TO FORM THE FOUR (2X10) MATRICES OF
! INTEGRALS (TRANSPOSED)
 
! COMPUTE CONSTANTS NEEDED IN FCURL  SUBROUTINE
 
 d(1) = 0.0E0
 IF (itord == 0) d(1) = 1.0E0 / rp
 
! NOTE THE DOUBLE SUBSCRIPTING USED IN  FCURL  SUBROUTINE IS
! COMPATIBLE WITH THE CALLING PROGRAM. THE DELINT ARRAY OF INTEGRALS
! IS A ( 7X6) SINGLY SUBSCRIPTED ARRAY (STORED ROWWISE) IN THE CALLING
! PROGRAM AND IT IS A (6X 7) DOUBLY SUBSCRIPTED ARRAY (STORED
! COLUMNWISE) IN FCURL  ROUTINE.
 
 
 CALL fcurl (fme(1), fme(21), ffe(1), ffe(21), delint(1), s, d(1))
 
 d(1) = twopi * tm
 d(2) = twopi * (tf **3) / 12.0E0
 DO  i = 1,40
   fme(i) = d(1) * fme(i)
   ffe(i) = d(2) * ffe(i)
 END DO
 
 
! FORM THE THERMAL STRAINS
 
 dtm1 = ti(1) - tz
 dtm2 = ti(2) - ti(1)
 dtf1 = 0.0E0
 dtf2 = 0.0E0
 
! THE TERMS DTF1 AND DTF2 ARE FUNCTIONS OF THE FLEXURAL GRADIENT
! TEMPERATURE BUT SINCE THESE TEMPERATURES ARE NOT AVAILABLE
! THE TERMS WILL BE SET TO ZERO. THEY ARE USUALLY DEFINED AS FOLLOWS,
!     DTF1 = TF(1) - TZ
!     DTF2 = TF(2) - TF(1)
! WHERE TF(1) AND TF(2) ARE THE FLEXURAL GRADIENT TEMPERATURES AT
! GRID POINTS 1 AND 2 RESPECTIVELY.
 
 d(1) = dtm1 * alf(1)
 d(2) = dtm1 * alf(2)
 d(3) = dtm2 * alf(1)
 d(4) = dtm2 * alf(2)
 d(5) = dtf1 * alf(1)
 d(6) = dtf1 * alf(2)
 d(7) = dtf2 * alf(1)
 d(8) = dtf2 * alf(2)
 
 
! FORM THE   THERMAL LOAD   IN FIELD COORDINATES
 
 CALL gmmats (ee(1), 2, 2, 0, d(1), 2, 1, 0, d(11) )
 CALL gmmats (ee(1), 2, 2, 0, d(3), 2, 1, 0, d(13) )
 CALL gmmats (ee(1), 2, 2, 0, d(5), 2, 1, 0, d(15) )
 CALL gmmats (ee(1), 2, 2, 0, d(7), 2, 1, 0, d(17) )
 
 
 CALL gmmats (fme( 1),  2,10, 1, d(11), 2, 1, 0, tl(1) )
 CALL gmmats (fme(21),  2,10,-1, d(13), 2, 1, 0, tl(1) )
 CALL gmmats (ffe( 1),  2,10,-1, d(15), 2, 1, 0, tl(1) )
 CALL gmmats (ffe(21),  2,10,-1, d(17), 2, 1, 0, tl(1) )
 
 
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
 
 
! TRANSFORM THE   THERMAL LOAD   TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmats(gambq(1), 10, 12, 1, tl(1), 10,  1, 0, d(1) )
 
 
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
 
 
 
! TRANSFORM THE   THERMAL LOAD   FROM ELEMENT TO BASIC COORDINATES
 
 CALL gmmats(gamrs(1), 12, 12, 1,  d(1), 12,  1, 0, tl(1) )
 
 
! LOCATE THE TRANSFORMATION MATRICES FROM BASIC TO LOCAL COORDINATES
! FOR THE TWO GRID POINTS AND EXPAND TO (6X6)
 
 DO  i = 1,144
   gambl(i) = 0.0E0
 END DO
 DO  i = 1,2
   CALL gbtran(ics(i),ecpt(4*i+10),d(1))
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
 
 
 
! TRANSFORM THE   THERMAL LOAD   FROM BASIC TO LOCAL COORDINATES
 
 CALL gmmats (gambl(1), 12, 12, 1, tl(1), 12,  1, 0, d(1) )
 
 
 
! ADD THE ELEMENT THERMAL LOAD TO THE STRUCTURE THERMAL LOAD
 
 k = 0
 DO  i = 1,2
   l = igp(i) - 1
   DO  j = 1,6
     k = k + 1
     l = l + 1
     pg(l) = pg(l) +  d(k)
   END DO
 END DO
 
 
 RETURN
END SUBROUTINE ttordr
