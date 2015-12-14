SUBROUTINE stpax1
     
!     THIS ROUTINE IS PHASE I OF STRESS DATA RECOVERY FOR THE AXI-
!     SYMMETRIC WITH A TRAPEZOIDAL CROSS SECTION
 
 
!     ECPT (01) = ELEMENT ID                                I
!     ECPT (02) = SIL A                                     I
!     ECPT (03) = SIL B                                     I
!     ECPT (04) = SIL C                                     I
!     ECPT (05) = SIL D
!     ECPT (06) = MATERIAL ORIENTATION ANGLE(DEGREES)       R
!     ECPT (08) = MATERIAL ID                               I
!     ECPT (09) TO ECPT (22) FOR PHI
!     ECPT (23) = COOR. SYS. FOR GRID POINT A               I
!     ECPT (24) = X-COOR. OF GRID POINT A (IN BASIC COOR)   R
!     ECPT (25) = Z-COOR. OF GRID POINT A (IN BASIC COOR)   R
!     ECPT (26) = 0.0
!     ECPT (27) = COOR. SYS. FOR GRID POINT B
!     ECPT (28) = X-COOR. OF GRID POINT B (IN BASIC COOR)   R
!     ECPT (29) = Z-COOR. OF GRID POINT B (IN BASIC COOR)   R
!     ECPT (30) = 0.0
!     ECPT (31) = COOR. SYS. FOR GRID POINT C               I
!     ECPT (32) = X-COOR. FOR GRID POINT C                  R
!     ECPT (33) = Z-COOR. FOR GRID POINT C                  R
!     ECPT (34) = 0.0
!     ECPT (35) = COOR. SYS. FOR GRID POINT D               I
!     ECPT (36) = X-COOR FOR GRID POINT D                   R
!     ECPT (37) = Z-COOR FOR GRID POINT D                   R
!     ECPT (38) = 0.0
!     ECPT (39) = EL. TEMPERATURE FOR MATERIAL PROP         R
 
!     ANY GROUP OF STATEMENTS PREFACED BY AN IF STATEMENT CONTAINING
!     ...KSYS78 OR LSYS78 ...  INDICATES CODING NECESSARY FOR THIS
!     ELEMENT*S PIEZOELECTRIC CAPABILITY
 
!     KSYS78 = 0   ELASTIC, NON-PIEZOELECTRIC MATERIAL
!     KSYS78 = 1   ELECTRICAL-ELASTIC COUPLED, PIEZOELETRIC MATERIAL
!     KSYS78 = 2   ELASTIC ONLY, PIEZOELECTRIC MATERIAL
!     LSYS78 = .TRUE. IF KSYS78 = 0, OR 2
 
 LOGICAL :: pzmat,lsys78
 INTEGER :: sp(50)
 DIMENSION       iecpt(40),delint(12),teo(45),acurl(208),  &
     ics(4),d1(48),d2(16),acurp1(48),acurp2(16),  &
     gababq(12,12),gbp(4,4),alfb(6),ee(63),wjp(3,4)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 
!     ECPT COMMON BLOCK
 
 COMMON /sdr2x5/ ecpt(39),dum5(61),idel,igp(4),tz,sel(360),ts(06),  &
     ak(144),phi(14),akph2(16),akuph(48),selp1(120), selp2(180),selp3(60)
 COMMON /sdr2x6/ d(144),e1(36),wj(6,12),r(5),z(5)
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e(3),anu(3),rho,g(3),alf(3),tzero,gsube,moskp(9), setmat
 COMMON /matpz / pzout(51)
!     COMMON /MATPZ / CE11,CE12,CE13,CE14,CE15,CE16,CE22,CE23,CE24,CE25,
!                     CE26,CE33,CE34,CE35,CE36,CE44,CE45,CE46,CE55,CE56,
!                     CE66,E11,E12,E13,E14,E15,E16,E21,E22,E23,E24,E25,
!                     E26,E31,E32,E33,E34,E35,E36,EPS11,EPS12,EPS13,
!                     EPS22,EPS23,EPS33,RHO,A1,A2,A12,TREF,GE
 COMMON /system/ ibuf,iout,dum75(75),ksys78
 COMMON /condas/ consts(5)
 EQUIVALENCE     (consts(1),pi),(consts(2),twopi),  &
     (consts(4),degrad),(acurl(1),ak(1)),  &
     (iecpt(1),ecpt(1)),(r(1),r1),(r(2),r2), (r(3),r3),(r(4),r4),(z(1),z1),  &
     (z(2),z2),(z(3),z3),(z(4),z4), (acurp1(1),acurl(145)),(acurp2(1),acurl(193))
 
 lsys78 = .false.
 IF (ksys78 == 0 .OR. ksys78 == 2) lsys78 = .true.
 
!     START EXECUTION
 
!     STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel   = iecpt( 1)
 igp(1) = iecpt( 2)
 igp(2) = iecpt( 3)
 igp(3) = iecpt( 4)
 igp(4) = iecpt( 5)
 matid  = iecpt( 8)
 ics(1) = iecpt(23)
 ics(2) = iecpt(27)
 ics(3) = iecpt(31)
 r(1)   =  ecpt(24)
 d(1)   =  ecpt(26)
 z(1)   =  ecpt(25)
 r(2)   =  ecpt(28)
 z(2)   =  ecpt(29)
 d(2)   =  ecpt(30)
 r(3)   =  ecpt(32)
 z(3)   =  ecpt(33)
 d(3)   =  ecpt(34)
 ics(4) = iecpt(35)
 z(4)   =  ecpt(37)
 d(4)   =  ecpt(38)
 r(4)   =  ecpt(36)
 tempe  =  ecpt(39)
 dgama  =  ecpt( 6)
 
!     TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 DO  i = 1,4
   IF (r(i) <= 0.0) GO TO 910
   IF (d(i) /= 0.0) GO TO 910
 END DO
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = AMIN1(z1,z2,z3,z4)
 z1   = z1 - zmin
 z2   = z2 - zmin
 z3   = z3 - zmin
 z4   = z4 - zmin
 rmin = AMIN1(r1,r2,r3,r4)
 rmax = AMAX1(r1,r2,r3,r4)
 IF (rmax/rmin <= 10.) GO TO 206
 
!     RATIO OF RADII IS TOO LARGE FOR GAUSS QUADRATURE FOR IP=-1
 
 idel1 = idel/1000
 WRITE  (iout,205) ufm,idel1
 205 FORMAT (a23,', TRAPAX ELEMENT',i9,' HAS A MAXIMUM TO MINIMUM ',  &
     'RADIUS RATIO EXCEEDING 10.', /5X,'ACCURACY OF NUMERICAL',  &
     ' INTEGRATION WOULD BE IN DOUBT.')
 GO TO 910
 206 CONTINUE
 
!     FORM THE TRANSFORMMATION MATRIX(12X12) FROM FIELD COOR, TO GRID
!     POINT DEGREES OF FREEDOM
 
 DO  i = 1,144
   gababq( i, 1) = 0.0
 END DO
 gababq( 1, 1) = 1.0
 gababq( 2, 1) = r1
 gababq( 3, 1) = z1
 gababq( 4, 1) = r1*z1
 gababq( 5, 2) = 1.0
 gababq( 6, 2) = r1
 gababq( 7, 2) = z1
 gababq( 8, 2) = gababq(4,1)
 gababq( 9, 3) = 1.0
 gababq(10, 3) = r1
 gababq(11, 3) = z1
 gababq(12, 3) = gababq(4,1)
 gababq( 1, 4) = 1.0
 gababq( 2, 4) = r2
 gababq( 3, 4) = z2
 gababq( 4, 4) = r2*z2
 gababq( 5, 5) = 1.0
 gababq( 6, 5) = r2
 gababq( 7, 5) = z2
 gababq( 8, 5) = gababq(4,4)
 gababq( 9, 6) = 1.0
 gababq(10, 6) = r2
 gababq(11, 6) = z2
 gababq(12, 6) = gababq(4,4)
 gababq( 1, 7) = 1.0
 gababq( 2, 7) = r3
 gababq( 3, 7) = z3
 gababq( 4, 7) = r3*z3
 gababq( 5, 8) = 1.0
 gababq( 6, 8) = r3
 gababq( 7, 8) = z3
 gababq( 8, 8) = gababq(4,7)
 gababq( 9, 9) = 1.0
 gababq(10, 9) = r3
 gababq(11, 9) = z3
 gababq(12, 9) = gababq(4,7)
 gababq( 1,10) = 1.0
 gababq( 2,10) = r4
 gababq( 3,10) = z4
 gababq( 4,10) = r4*z4
 gababq( 5,11) = 1.0
 gababq( 6,11) = r4
 gababq( 7,11) = z4
 gababq( 8,11) = gababq(4,10)
 gababq( 9,12) = 1.0
 gababq(10,12) = r4
 gababq(11,12) = z4
 gababq(12,12) = gababq(4,10)
 
 IF (lsys78) GO TO 305
 gbp(1,1) = 1.0
 gbp(2,1) = r(1)
 gbp(3,1) = z(1)
 gbp(4,1) = r(1)*z(1)
 gbp(1,2) = 1.0
 gbp(2,2) = r(2)
 gbp(3,2) = z(2)
 gbp(4,2) = r(2)*z(2)
 gbp(1,3) = 1.0
 gbp(2,3) = r(3)
 gbp(3,3) = z(3)
 gbp(4,3) = r(3)*z(3)
 gbp(1,4) = 1.0
 gbp(2,4) = r(4)
 gbp(3,4) = z(4)
 gbp(4,4) = r(4)*z(4)
 305 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (12,gababq,12,d(10),0,d(11),ising,sp)
 IF (ising == 2) GO TO 920
 
 IF (ksys78 == 1) CALL invers (4,gbp,4,d(10),0,d(11),ising,sp)
 IF (ising  == 2) GO TO 920
 
!     MODIFY THE TRANSFORMATION MATRIX IF ELEMENT IS A CORE ELEMENT
 
!     CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT
 
!       DELINT(1) = (-1,0)
!       DELINT(02)= (-1,1)
!       DELINT(03)= (-1,2)
!       DELINT(04)= ( 0,0)
!       DELINT(05)= ( 0,1)
!       DELINT(06)= ( 0,2)
!       DELINT(07)= ( 1,0)
!       DELINT(08)= ( 1,1)
!       DELINT(09)= ( 1,2)
!       DELINT(10)= ( 2,0)
!       DELINT(11)= ( 2,1)
!       DELINT(12)= ( 3,0)
 
 i1 = 0
 DO  i = 1,4
   ip = i - 2
   DO  j = 1,3
     iq = j - 1
     i1 = i1 + 1
     IF (i1 /= 12) GO TO 340
     ip = 3
     iq = 0
     340 CONTINUE
     delint(i1) = rzints(ip,iq,r,z,4)
   END DO
 END DO
 
!     LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3
 
 matidc = matid
 matflg = 7
 IF (ksys78 > 0) matflg = 9
 eltemp = tempe
 dgamr  = dgama*degrad
 sinth  = SIN(dgamr)
 costh  = COS(dgamr)
 cosg   = costh
 sing   = sinth
 CALL mat (idel)
 pzmat  = .false.
 IF (setmat == 4. .OR. setmat == 5.) pzmat = .true.
 IF (pzmat) GO TO 410
 ksave  = ksys78
 ksys78 = 0
 lsys78 = .true.
 GO TO 420
 410 rho    = pzout(46)
 alf(1) = pzout(47)
 alf(2) = pzout(48)
 alf(3) = pzout(49)
 tzero  = pzout(50)
 gsube  = pzout(51)
 420 CONTINUE
 IF (setmat == 2.0) GO TO 915
 tz = tzero
 IF (ksys78 > 0) GO TO 500
 
!     SET MATERIAL PROPERTIES IN DOUBLE PRECISION VARIABLES
 
 er  = e(1)
 et  = e(2)
 ez  = e(3)
 vro = anu(1)
 voz = anu(2)
 vzr = anu(3)
 gor = g(1)
 gzo = g(2)
 grz = g(3)
 vor = vro*et/er
 vzo = voz*ez/et
 vrz = vzr*er/ez
 del = 1.0/(1.0 - vro*vor - voz*vzo - vzr*vrz - vro*voz*vzr - vrz*vor*vzo)
 
!     COMPUTE ELASTIC CONSTANTS MATRIX FROM MATERIAL TO ELEMENT AXIS
 
 500 CONTINUE
 DO  i = 1,45
   teo (i) = 0.0
 END DO
 
 IF (ksys78 > 0) GO TO 520
 teo ( 1) = er*(1.0 - voz*vzo)*del
 teo ( 2) = er*(vzr + vzo*vor)*del
 teo ( 3) = ez*(1.0 - vro*vor)*del
 teo ( 4) = er*(vor + vzr*voz)*del
 teo ( 5) = et*(vzo + vro*vzr)*del
 teo ( 6) = et*(1.0 - vrz*vzr)*del
 teo (10) = grz
 teo (15) = gor
 teo (21) = gzo
 GO TO 530
 520 CONTINUE
 
!     PIEZOELECTRIC MATERIAL PROPERTIES STORED IN TEO(22-39)
!     DIELECTRIC MATERIAL PROPERTIES STORED IN TEO(40-45)
!     TEO(22-39) CONTAINS E-TRANSPOSE
 
 teo( 1) = pzout( 1)
 teo( 2) = pzout( 2)
 teo( 3) = pzout( 7)
 teo( 4) = pzout( 3)
 teo( 5) = pzout( 8)
 teo( 6) = pzout(12)
 teo( 7) = pzout( 4)
 teo( 8) = pzout( 9)
 teo( 9) = pzout(13)
 teo(10) = pzout(16)
 teo(11) = pzout( 5)
 teo(12) = pzout(10)
 teo(13) = pzout(14)
 teo(14) = pzout(17)
 teo(15) = pzout(19)
 teo(16) = pzout( 6)
 teo(17) = pzout(11)
 teo(18) = pzout(15)
 teo(19) = pzout(18)
 teo(20) = pzout(20)
 teo(21) = pzout(21)
 
 IF (ksys78 == 2) GO TO 530
 teo(22) = pzout(22)
 teo(23) = pzout(28)
 teo(24) = pzout(34)
 teo(25) = pzout(23)
 teo(26) = pzout(29)
 teo(27) = pzout(35)
 teo(28) = pzout(24)
 teo(29) = pzout(30)
 teo(30) = pzout(36)
 teo(31) = pzout(25)
 teo(32) = pzout(31)
 teo(33) = pzout(37)
 teo(34) = pzout(26)
 teo(35) = pzout(32)
 teo(36) = pzout(38)
 teo(37) = pzout(27)
 teo(38) = pzout(33)
 teo(39) = pzout(39)
 teo(40) =-pzout(40)
 teo(41) =-pzout(41)
 teo(42) =-pzout(42)
 teo(43) =-pzout(43)
 teo(44) =-pzout(44)
 teo(45) =-pzout(45)
 530 CONTINUE
 
 c2   = cosg*cosg
 c4   = c2*c2
 s2   = sing*sing
 s4   = s2*s2
 c2s2 = c2*s2
 c3   = cosg*c2
 s3   = sing*s2
 cs2  = cosg*s2
 sc2  = sing*c2
 cs   = cosg*sing
 
 ee( 1) = teo(1)*c4 + teo(3)*s4 + 2.0*c2s2 * (teo(2) + 2.0 * teo(10))
 ee( 2) = teo(2)*(c4+s4) + c2s2 * (teo(1)+teo(3)-4.0*teo(10))
 ee( 3) = teo(1)*s4 + 2.0*c2s2 * (teo(2)+2.0*teo(10)) +  teo(3)*c4
 ee( 4) = teo(4)*c2 + teo(5)*s2
 ee( 5) = teo(4)*s2 + teo(5)*c2
 ee( 6) = teo(6)
 ee( 7) = cosg*sing*s2 * (teo(2)-teo(3)+2.0*teo(10))  &
     + sing*cosg*c2 * (teo(1)-teo(2)-2.0*teo(10))
 ee( 8) = sing*cosg*c2 * (teo(2)-teo(3)+2.0*teo(10))  &
     + cosg*sing*s2 * (teo(1)-teo(2)-2.0*teo(10))
 ee( 9) = sing*cosg * (teo(4) - teo(5))
 ee(10) = c2s2 * (teo(1)-2.0*teo(2)+teo(3)) + teo(10)*(c2-s2)**2
 ee(12) = 0.0
 ee(13) = 0.0
 ee(15) = teo(15)*c2 + teo(21)*s2
 ee(20) = cosg*sing * (teo(15) - teo(21))
 ee(21) = teo(15)*s2 + teo(21)*c2
 
 IF (lsys78) GO TO 540
 
!     PIEZOELECTRIC MATERIAL PROPERTIES IN ELEMENT COORDINATES
 
 ee(37) = c3*teo(22) - s3*teo(26) + cs2*(teo(25)+2.0*teo(32)) -  &
     sc2*(teo(23)+2.0*teo(31))
 ee(38) = c3*teo(23) + s3*teo(25) + cs2*(teo(26)-2.0*teo(31)) +  &
     sc2*(teo(22)-2.0*teo(32))
 ee(39) = s2*teo(27) + c2*teo(24) - 2.0*cs*teo(33)
 ee(40) = c3*teo(25) - s3*teo(23) + cs2*(teo(22)-2.0*teo(32)) -  &
     sc2*(teo(26)-2.0*teo(31))
 ee(41) = c3*teo(26) + s3*teo(22) + cs2*(teo(23)+2.0*teo(31)) +  &
     sc2*( teo(25)+2.0*teo(32))
 ee(42) = s2*teo(24) + c2*teo(27) + 2.0*cs*teo(33)
 ee(43) = cosg*teo(28) - sing*teo(29)
 ee(44) = cosg*teo(29) + sing*teo(28)
 ee(45) = teo(30)
 ee(46) = c3*teo(31) + s3*teo(32) - cs2*(teo(23)-teo(26)+teo(31)) +  &
     sc2*(-teo(32)-teo(25)+teo(22))
 ee(47) = c3*teo(32) - s3*teo(31) - cs2*(teo(25)-teo(22)+teo(32)) +  &
     sc2*(teo(23)+teo(31)-teo(26))
 ee(48) = (c2-s2)*teo(33) + cs*(teo(24)-teo(27))
 ee(49) = c2*teo(34) + s2*teo(38) - cs*(teo(35)+teo(37))
 ee(50) = c2*teo(35) - s2*teo(37) + cs*(teo(34)-teo(38))
 ee(51) = cosg*teo(36) - sing*teo(39)
 ee(52) = c2*teo(37) - s2*teo(35) - cs*(teo(38)-teo(34))
 ee(53) = c2*teo(38) + s2*teo(34) + cs*(teo(35)+teo(37))
 ee(54) = cosg*teo(39) + sing*teo(36)
 
!     DIELECTRIC MATERIAL PROPERTIES IN ELEMENT COORDINTES
 
 ee(55) = s2*teo(43) - 2.0*cs*teo(41) + c2*teo(40)
 ee(56) = (c2-s2)*teo(41) - cs*(teo(43)-teo(40))
 ee(57) =-sing*teo(44) + cosg*teo(42)
 ee(59) = c2*teo(43) + 2.0*cs*teo(41) + s2*teo(40)
 ee(60) = cosg*teo(44) + sing*teo(42)
 ee(63) = teo(45)
 540 CONTINUE
 
!     COMPUTE HARMONIC COEFFICIENT
 
 iecpt(1) = iecpt(1) - (iecpt(1)/1000)*1000 - 1
 ajho  = iecpt(1)
 ajjho = ajho*ajho
 
!     FORM THE ELEMENT STIFFNESS MATRIX IN FIELD SYSTEM
 
 acurl ( 1) = (ee(6) + ajjho * ee(15)) * delint(1)
 acurl ( 2) = (ee(4) + ee(6) + ajjho * ee(15)) * delint(4)
 acurl ( 3) = (ee(6) + ajjho * ee(15)) * delint(2) +  ee(9) * delint (4)
 acurl ( 4) = (ee(4) + ee(6) + ajjho * ee(15))* delint(5) +  ee(9) * delint(7)
 acurl ( 5) = ajho * (ee(6) + ee(15)) * delint(1)
 acurl ( 6) = ajho * ee(6) * delint(4)
 acurl ( 7) = ajho * (ee(6) + ee(15))* delint(2) - ajho * ee(20) * delint(4)
 acurl ( 8) = ajho * ee(6) * delint(5) - ajho * ee(20) * delint(7)
 acurl ( 9) = ajjho * ee(20) * delint(1)
 acurl (10) = delint(4) * (ee(9) + ajjho*ee(20))
 acurl (11) = delint(4) * ee(5) + ajjho * delint(2) * ee(20)
 acurl (12) = delint(7) * ee(5) + delint(5)*(ee(9)+ajjho*ee(20))
 acurl (14) = (ee(1) + 2.0 * ee(4) + ee(6) + ajjho * ee(15)) * delint(7)
 acurl (15) = (ee(4) + ee(6) + ajjho * ee(15)) * delint(5)  &
     + (ee(7) + ee(9)) * delint(7)
 acurl (16) = (ee(1) + 2.0 * ee(4) + ajjho * ee(15) + ee(6))  &
     * delint(8) + (ee(7) + ee(9)) * delint(10)
 acurl (17) = ajho * (ee(4) + ee(6) + ee(15)) * delint(4)
 acurl (18) = ajho * (ee(4) + ee(6)) * delint(7)
 acurl (19) = ajho * (ee(4) + ee(6) + ee(15)) * delint(5) - ajho  &
     * ee(20) * delint(7)
 acurl (20) = ajho * (ee(4) + ee(6)) * delint(8) - ajho * ee(20) * delint(10)
 acurl (21) = ajjho * ee(20) * delint(4)
 acurl (22) = delint(7) * (ee(7) + ee(9) + ajjho*ee(20))
 acurl (23) = delint(7)*(ee(2)+ee(5))+ajjho*delint(5)*ee(20)
 acurl (24) = delint(10)*(ee(2)+ee(5))+delint(8)*(ee(7)+ee(9)) +  &
     delint(8)*ajjho*ee(20)
 acurl (27) = (ee(6) + ajjho * ee(15)) * delint(3) + 2.0  &
     * ee(9) * delint(5) + ee(10) * delint (7)
 acurl (28) = (ee(4) + ee(6) + ajjho * ee(15)) * delint(6)  &
     + ee(10) * delint(10) + (ee(7) + 2.0 * ee(9)) * delint (8)
 acurl (29) = ajho * (ee(6) + ee(15)) * delint(2) + ajho * ee(9) * delint(4)
 acurl (30) = ajho * ee(6) * delint(5) + ajho * ee(9) * delint(7)
 acurl (31) = ajho * (ee(6) + ee(15)) * delint(3) + ajho * (ee(9)  &
     - ee(20)) * delint(5)
 acurl (32) = ajho * (ee(9) - ee(20))* delint(8) + ajho * ee(6) * delint(6)
 acurl (33) = ajjho * ee(20) * delint(2)
 acurl (34) = delint(7)*ee(10) + delint(5)*(ee(9) + ajjho*ee(20))
 acurl (35) = delint(7)*ee(8) + delint(5)*ee(5) + ajjho*delint(3) * ee(20)
 acurl (36) = delint(10)*ee(8)+delint(8)*(ee(5)+ee(10)) +  &
     delint(6)*(ee(9)+ajjho*ee(20))
 acurl (40) = (ee(1) + 2.0 * ee(4) + ee(6) + ajjho * ee(15))  &
     * delint(9) + (2.0 * ee(7) + 2.0 * ee(9)) * delint(11)+ ee(10) * delint(12)
 acurl (41) = ajho * (ee(4) + ee(6) + ee(15)) * delint(5)  &
     + ajho * ee(9) * delint(7)
 acurl (42) = ajho * (ee(4) + ee(6)) * delint(8) + ajho * ee(9) * delint(10)
 acurl (43) = ajho * (ee(4) + ee(6) + ee(15))* delint(6)  &
     + ajho * (ee(9) - ee(20)) * delint(8)
 acurl (44) = ajho * (ee(4) + ee(6)) * delint(9) + ajho  &
     * (ee(9) - ee(20)) * delint(11)
 acurl (45) = ajjho * ee(20) * delint(5)
 acurl (46) = delint(8)*(ee(7) + ee(9) + ajjho*ee(20)) + delint(10) * ee(10)
 acurl (47) = delint(8)*(ee(2) + ee(5)) + delint(10)*ee(8) +  &
     ajjho*delint(6)*ee(20)
 acurl (48) = delint(11)*(ee(2)+ee(5)+ee(10)) + delint(12)*ee(8) +  &
     delint(9)*(ee(7)+ee(9)+ajjho*ee(20))
 acurl (53) = (ee(15) + ajjho * ee(6)) * delint(1)
 acurl (54) = ajjho * ee(6) * delint(4)
 acurl (55) = (ee(15) + ajjho * ee(6)) * delint(2) - ee(20) * delint(4)
 acurl (56) = ajjho * ee(6) * delint(5) - ee(20) * delint(7)
 acurl (57) = ajho * ee(20) * delint(1)
 acurl (58) = ajho*delint(4)*(ee(9)+ee(20))
 acurl (59) = ajho*(delint(4)*ee(5) + delint(2)*ee(20))
 acurl (60) = ajho*(delint(7)*ee(5)+delint(5)*(ee(9)+ee(20)))
 acurl (66) = ajjho * ee(6) * delint(7)
 acurl (67) = ajjho * ee(6) * delint(5)
 acurl (68) = ajjho * ee(6) * delint(8)
 acurl (69) = 0.0
 acurl (70) = ajho*delint(7)*ee(9)
 acurl (71) = ajho*delint(7)*ee(5)
 acurl (72) = ajho*(delint(10)*ee(5)+delint(8)*ee(9))
 acurl (79) = (ee(15) + ajjho * ee(6)) * delint(3) - 2.0  &
     * ee(20) * delint(5) + ee(21) * delint(7)
 acurl (80) = ajjho * ee(6) * delint(6) - ee(20) * delint(8)  &
     + ee(21) * delint(10)
 acurl (81) = ajho * (ee(20) * delint(2) - ee(21) * delint(4))
 acurl (82) = ajho*(delint(5)*(ee(9)+ee(20))-delint(7)*ee(21))
 acurl (83) = ajho*(delint(5)*(ee(5)-ee(21))+delint(3)*ee(20))
 acurl (84) = ajho*(delint(8)*(ee(5)-ee(21))+delint(6)*(ee(9) + ee(20)))
 acurl (92) = ee(21) * delint(12) + ajjho * ee(6) * delint(9)
 acurl (93) =-ajho * ee(21) * delint(7)
 acurl (94) = ajho*(delint(8)*ee(9)-delint(10)*ee(21))
 acurl (95) = ajho* delint(8) * (ee(5)-ee(21))
 acurl (96) = ajho*(delint(11)*(ee(5)-ee(21))+delint(9)*ee(9))
 acurl (105) = ajjho * ee(21) * delint(1)
 acurl (106) = ajjho*delint(4)*ee(21)
 acurl (107) = ajjho*delint(2)*ee(21)
 acurl (108) = ajjho*delint(5)*ee(21)
 acurl (118) = delint(7)*(ee(10)+ajjho*ee(21))
 acurl (119) = delint(7)*ee(8)+ajjho*delint(5)*ee(21)
 acurl (120) = delint(10)*ee(8)+delint(8)*(ee(10)+ajjho*ee(21))
 acurl (131) = delint(7)*ee(3)+ajjho*delint(3)*ee(21)
 acurl (132) = delint(10)*ee(3)+delint(8)*ee(8)+ajjho*delint(6)* ee(21)
 acurl (144) = delint(12)*ee(3) + 2.0*delint(11)*ee(8) +  &
     delint(9)*(ee(10)+ajjho*ee(21))
 
 IF (lsys78) GO TO 550
 acurl(145) = delint(1)*ajho*(ajho*ee(51)-ee(45))
 acurl(146) = delint(4)*(ee(43)+ajho*(ajho*ee(51)-ee(49)-ee(45)))
 acurl(147) = delint(2)*ajho*(ajho*ee(51)-ee(45))+delint(4)*  &
     (ee(44)-ajho*ee(50))
 acurl(148) = delint(5)*(ee(43)+ajho*(ajho*ee(51)-ee(49)-ee(45)))  &
     + delint(7)*(ee(44)-ajho*ee(50))
 acurl(149) = delint(4)*ajho*(ajho*ee(51)-ee(45)-ee(39))
 acurl(150) = delint(7)*(ee(43)+ee(37)+ajho*(ajho*ee(51)-ee(49)  &
     - ee(45)-ee(39)))
 acurl(151) = delint(5)*ajho*(ajho*ee(51)-ee(45)-ee(39))+delint(7)  &
     * (ee(44)+ee(38)-ajho*ee(50))
 acurl(152) = delint(8)*(ee(43)+ee(37)+ajho*(ajho*ee(51)-ee(49)-  &
     ee(45)-ee(39)))+delint(10)*(ee(44)+ee(38)-ajho* ee(50))
 acurl(153) = delint(2)*ajho*(ajho*ee(51)-ee(45))-delint(4)*ajho * ee(48)
 acurl(154) = delint(5)*(ee(43)+ajho*(ajho*ee(51)-ee(49)-ee(45)))  &
     + delint(7)*(ee(46)-ajho*ee(48))
 acurl(155) = delint(3)*ajho*(ajho*ee(51)-ee(45))+delint(5)*  &
     (ee(44)-ajho*(ee(50)+ee(48)))+delint(7)*ee(47)
 acurl(156) = delint(6)*(ee(43)+ajho*(ajho*ee(51)-ee(49)-ee(45)))  &
     + delint(8)*(ee(46)+ee(44)-ajho*(ee(50)+ee(48)))+ delint(10)*ee(47)
 acurl(157) = delint(5)*ajho*(ajho*ee(51)-ee(45)-ee(39))-delint(7)  &
     * ajho*ee(48)
 acurl(158) = delint(8)*(ee(43)+ee(47)+ajho*(ajho*ee(51)-ee(49)-  &
     ee(45)-ee(39)))-delint(10)*(ee(46)-ajho*ee(48))
 acurl(159) = delint(6)*ajho*(ajho*ee(51)-ee(45)-ee(39))+delint(8)  &
     * (ee(44)+ee(38)-ajho*(ee(50)+ee(48)))+delint(10)* ee(47)
 acurl(160) = delint(9)*(ee(43)+ee(37)+ajho*(ajho*ee(51)-ee(49)-  &
     ee(45)-ee(39)))+delint(11)*(ee(46)+ee(44)+ee(38)-  &
     ajho*(ee(50)+ee(48)))+delint(12)*ee(47)
 acurl(161) = delint(1)*ajho*(ee(51)-ajho*ee(45))
 acurl(162) = delint(4)*(-ee(49)+ajho*(ee(51)+ee(43)-ajho*ee(45)))
 acurl(163) = delint(2)*ajho*(ee(51)-ajho*ee(45))+delint(4)*  &
     (ajho*ee(44)-ee(50))
 acurl(164) = delint(5)*(-ee(49)+ajho*(ee(51)+ee(43)-ajho*ee(51)))  &
     + delint(7)*(ajho*ee(44)-ee(50))
 acurl(165) =-delint(4)*ajjho*ee(45)
 acurl(166) = delint(7)*ajho*(ee(43)-ajho*ee(45))
 acurl(167) = delint(7)*ajho*ee(44)-delint(5)*ajjho*ee(45)
 acurl(168) = delint(8)*ajho*(ee(43)-ajho*ee(45))+delint(10)* ajho*ee(44)
 acurl(169) = delint(2)*ajho*(ee(51)-ajho*ee(45))-delint(4)*ajho* ee(54)
 acurl(170) = delint(5)*(-ee(49)+ajho*(ee(51)+ee(43)-ajho*ee(45)))  &
     + delint(7)*(ee(52)-ajho*ee(54))
 acurl(171) = delint(3)*ajho*(ee(51)-ajho*ee(45))+delint(5)*  &
     (ajho*(ee(44)-ee(54))-ee(50))+delint(7)*ee(53)
 acurl(172) = delint(6)*(-ee(49)+ajho*(ee(51)+ee(43)-ajho*ee(45)))  &
     + delint(8)*(ee(52)-ee(50)+ajho*(ee(44)-ee(54))) + delint(10)*ee(53)
 acurl(173) =-delint(5)*ajjho*ee(45)-delint(7)*ajho*ee(54)
 acurl(174) = delint(8)*ajho*(ee(43)-ajho*ee(45))+delint(10)*  &
     (ee(54)-ajho*ee(54))
 acurl(175) =-delint(6)*ajjho*ee(45)+delint(8)*ajho*(ee(44)-  &
     ee(54))+delint(10)*ee(53)
 acurl(176) = delint(9)*ajho*(ee(43)-ajho*ee(45))+delint(11)*  &
     (ee(52)+ajho*(ee(44)-ee(54)))+delint(12)*ee(53)
 acurl(177) = delint(1)*ajjho*ee(54)
 acurl(178) = delint(4)*ajho*(ajho*ee(54)-ee(52))
 acurl(179) = delint(2)*ajjho*ee(54)-delint(4)*ajho*ee(53)
 acurl(180) = delint(5)*ajho*(ajho*ee(54)-ee(52))-delint(7)*ajho * ee(53)
 acurl(181) = delint(4)*ajho*(ajho*ee(54)-ee(48))
 acurl(182) = delint(7)*(ee(46)+ajho*(ajho*ee(54)-ee(52)-ee(48)))
 acurl(183) = delint(5)*ajho*(ajho*ee(54)-ee(48))+delint(7)*  &
     (ee(47)-ajho*ee(53))
 acurl(184) = delint(8)*(ee(46)+ajho*(ajho*ee(54)-ee(52)-ee(48)))  &
     + delint(10)*(ee(47)-ajho*ee(53))
 acurl(185) = delint(2)*ajjho*ee(54)-delint(4)*ajho*ee(42)
 acurl(186) = delint(5)*ajho*(ajho*ee(54)-ee(52))+delint(7)*(ee(40)  &
     - ajho*ee(42))
 acurl(187) = delint(3)*ajjho*ee(54)-delint(5)*ajho*(ee(53)+ee(42))  &
     + delint(7)*ee(41)
 acurl(188) = delint(6)*ajho*(ajho*ee(54)-ee(52))+delint(8)*  &
     (ee(40)-ajho*(ee(53)+ee(42)))+delint(10)*ee(41)
 acurl(189) =-delint(5)*ajho*ee(48)+delint(4)*ajjho*ee(54)  &
     - delint(7)*ajho*ee(42)
 acurl(190) = delint(8)*(ee(46)-ajho*ee(48))+delint(7)*ajho*  &
     (ajho*ee(54)-ee(52))+delint(10)*(ee(40)-ajho*ee(42))
 acurl(191) =-delint(6)*ajho*ee(48)+delint(5)*ajjho*ee(54)+  &
     delint(8)*(ee(47)-ajho*ee(42))-delint(7)*ajho*ee(53) + delint(10)*ee(41)
 acurl(192) = delint(9)*(ee(46)-ajho*ee(48))+delint(8)*ajho*  &
     (ajho*ee(54)-ee(52))+delint(11)*(ee(47)+ee(40)-  &
     ajho*ee(42))-delint(10)*ajho*ee(53)+delint(12)*ee(41)
 
 acurl(193) = delint(1)*ajjho*ee(63)
 acurl(194) = delint(4)*ajho*(ajho*ee(63)-ee(57))
 acurl(195) = delint(2)*ajjho*ee(63)-delint(4)*ajho*ee(60)
 acurl(196) = delint(5)*ajho*(ajho*ee(63)-ee(57))-delint(7)* ajho*ee(60)
 acurl(197) = delint(4)*ajho*(ajho*ee(63)-ee(57))
 acurl(198) = delint(7)*(ajjho*ee(63)-2.0*ajho*ee(57)+ee(55))
 acurl(199) = delint(5)*ajho*(ajho*ee(63)-ee(57))+delint(7)*(ee(56)  &
     - ajho*ee(60))
 acurl(200) = delint(8)*(ajjho*ee(63)-2.0*ajho*ee(57)+ee(55))  &
     + delint(10)*(ee(56)-ajho*ee(60))
 acurl(201) = delint(2)*ajjho*ee(63)-delint(4)*ajho*ee(60)
 acurl(202) = delint(5)*ajho*(ajho*ee(63)-ee(57))+delint(7)*  &
     (ee(56)-ajho*ee(60))
 acurl(203) = delint(3)*ajjho*ee(63)-delint(5)*2.0*ajho*ee(60)  &
     + delint(7)*ee(59)
 acurl(204) = delint(6)*ajho*(ajho*ee(63)-ee(57))+delint(8)*  &
     (ee(56)-2.0*ajho*ee(60))+delint(10)*ee(59)
 acurl(205) = delint(5)*ajho*(ajho*ee(63)-ee(57))-delint(7)* ajho*ee(60)
 acurl(206) = delint(8)*(ajjho*ee(63)-2.0*ee(57)+ee(55))+delint(10)  &
     * (ee(56)-ajho*ee(60))
 acurl(207) = delint(6)*ajho*(ajho*ee(63)-ee(57))+delint(8)*(ee(56)  &
     - 2.0*ajho*ee(60))+delint(10)*ee(59)
 acurl(208) = delint(9)*(ajjho*ee(63)-2.0*ajho*ee(57)+ee(55))+  &
     2.0*delint(11)*(ee(56)-ajho*ee(60))+delint(12)*ee(59)
 550 CONTINUE
 
!     TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM FIELD SYSTEM
!     TO GRID POINT DEGREES OF FREEDOM
 
!     EXPAND ACURL INTO (12X12)
 
 DO  ib = 2,12
   ic = 13*ib - 25
   i  = ic
   DO  j = ib,12
     ic = ic + 12
     i  = i + 1
     acurl(ic) = acurl(i)
   END DO
 END DO
 
 dgama = pi
 IF (ajho == 0.0) dgama = twopi
 DO  i = 1,144
   acurl(i) = acurl(i)*dgama
 END DO
 
 IF (lsys78) GO TO 638
 DO  i = 145,208
   acurl(i) = acurl(i)*dgama
 END DO
 638 CONTINUE
 
 CALL gmmats (gababq,12,12,1, acurl, 12,12,0, d )
 CALL gmmats (     d,12,12,0, gababq,12,12,0, ak)
 
 IF (lsys78) GO TO 639
 CALL gmmats (gababq,12,12,1, acurp1,12,4,0, d1)
 CALL gmmats (d1,12,4,0, gbp,4,4,0, akuph)
 CALL gmmats (gbp,4,4,1, acurp2,4,4,0, d2)
 CALL gmmats (d2,4,4,0, gbp,4,4,0, akph2)
 639 CONTINUE
 
!     ********** COORDINATE SYSTEM NOT POSSIBLE ***********************
!     *** WITH RINGAX.  THE FOLLOWING CODE WILL IMPLEMENT IT  *********
!     IF FOLLOWING CODE IS IMPLEMENTED MUST BE MODIFIED FOR PIEZO-
!     ELECTRIC
 
!     ZERO OUT THE AKI MATRIX
 
!.    DO 700 I = 1,144
!.    AKI(I) = 0.0
!.700 CONTINUE
!.    DO 800 I = 1,4
!.    CALL TRANSS (ICS(I),D(1))
!.    K = 39 * (I-1) + 1
!.    DO 800 J = 1, 3
!.    KK = K + 12*(J-1)
!.    JJ = 3*(J-1) + 1
!.    AKI (KK  ) = D (JJ  )
!.    AKI (KK+1) = D (JJ+1)
!.    AKI (KK+2) = D (JJ+2)
!.800 CONTINUE
 
!     TRANSFORM THE STIFFNESS MATRIX FROM BASIC TO LOCAL COORDINATES
 
!.    CALL GMMATS (AKI(1),12,12,1, AK(1),12,12,0, D(1))
!.    CALL GMMATS (D(1),12,12,0, AKI(1),12,12,0, AK(1))
 
!     COMPUTE THE FIFTH GRID POINT
 
 r(5) = (r1 + r2 + r3 + r4)/4.0
 z(5) = (z1 + z2 + z3 + z4)/4.0
 
!     FORM WJ MATRIX
 
 DO  iki = 1,5
   DO  i = 1,72
     wj(i, 1) = 0.0
   END DO
   rsum     = r(iki)
   zsum     = z(iki)
   zdr      = zsum/rsum
   wj(1, 2) = 1.0
   wj(1, 4) = zsum
   wj(2,11) = 1.0
   wj(2,12) = rsum
   wj(3, 1) = 1.0/rsum
   wj(3, 2) = 1.0
   wj(3, 3) = zdr
   wj(3, 4) = zsum
   wj(3, 5) = ajho/rsum
   wj(3, 6) = ajho
   wj(3, 7) = ajho*zdr
   wj(3, 8) = ajho*zsum
   wj(4, 3) = 1.0
   wj(4, 4) = rsum
   wj(4,10) = 1.0
   wj(4,12) = zsum
   wj(5, 1) =-ajho/rsum
   wj(5, 2) =-ajho
   wj(5, 3) =-ajho*zdr
   wj(5, 4) =-ajho*zsum
   wj(5, 5) =-1.0/rsum
   wj(5, 7) =-zdr
   wj(6, 7) = 1.0
   wj(6, 8) = rsum
   wj(6, 9) =-ajho/rsum
   wj(6,10) =-ajho
   wj(6,11) =-ajho*zdr
   wj(6,12) =-ajho*zsum
   
   IF (lsys78) GO TO 1060
   
!     FORM WJP MATRIX
   
   DO  i = 1,3
     DO  j = 1,4
       wjp(i,j) = 0.0
     END DO
   END DO
   
   wjp(1,2) = 1.0
   wjp(1,4) = zsum
   wjp(2,3) = 1.0
   wjp(2,4) = rsum
   wjp(3,1) =-ajho/rsum
   wjp(3,2) =-ajho
   wjp(3,3) =-ajho*zdr
   wjp(3,4) =-ajho*zsum
   1060 CONTINUE
   
!     EXPAND EE(21) INTO E1(36)
   
   DO  i = 1,36
     e1( i) = 0.0
   END DO
   e1( 1) = ee( 1)
   e1( 2) = ee( 2)
   e1( 3) = ee( 4)
   e1( 4) = ee( 7)
   e1( 7) = ee( 2)
   e1( 8) = ee( 3)
   e1( 9) = ee( 5)
   e1(10) = ee( 8)
   e1(13) = ee( 4)
   e1(14) = ee( 5)
   e1(15) = ee( 6)
   e1(16) = ee( 9)
   e1(19) = ee( 7)
   e1(20) = ee( 8)
   e1(21) = ee( 9)
   e1(22) = ee(10)
   e1(29) = ee(15)
   e1(36) = ee(21)
   
!     COMPUTE THE STRESS MATRICES
   
   k = 72*(iki-1) + 1
   CALL gmmats (wj,12,6,1, gababq,12,12,0, d(1))
   CALL gmmats (e1(1),6,6,0, d(1),6,12,0, sel(k))
   
   IF (lsys78) GO TO 1070
   kp1 = 24*(iki-1) + 1
   CALL gmmats (wjp,4,3,1, gbp,4,4,0, d2(1))
   CALL gmmats (ee(37),6,3,0, d2(1),3,4,0, selp1(kp1))
   kp2 = 36*(iki-1) + 1
   CALL gmmats (ee(37),6,3,1, d(1),6,12,0, selp2(kp2))
   kp3 = 12*(iki-1) + 1
   CALL gmmats (ee(55),3,3,0, d2(1),3,4,0, selp3(kp3))
   1070 CONTINUE
   
!     ** COORDINATE SYSTEMS NOT POSSIBLE WITH RINGAX *******************
!     ** THE FOLLOWING CODE WILL IMPLEMENT IT **************************
!     ** NOTE THAT WJ IS SEL(K) IN FOLLOWING GMMATS ********************
   
!     ** IF FOLLOWING CODE IS IMPLEMENTED MUST BE MODIFIED FOR PIEZO-
!     ELECTRIC TRANSFORM THE STRESS MATRIX FROM BASIC TO LOCAL
!     COORDINATES
!..   CALL GMMATS (WJ,6,12,0, AKI(1),12,12,0, SEL(K) )
   
   
 END DO
 
!     COMPUTE THE THERMAL STRAIN
 
 alfb(1) = alf(1)
 alfb(2) = alf(3)
 alfb(3) = alf(2)
 alfb(4) = 0.0
 alfb(5) = 0.0
 alfb(6) = 0.0
 
!     COMPUTE THE THERMAL STRESS
 
 ts(1) = ee(1)*alfb(1) + ee(2)*alfb(2) + ee(4)*alfb(3)
 ts(2) = ee(2)*alfb(1) + ee(3)*alfb(2) + ee(5)*alfb(3)
 ts(3) = ee(4)*alfb(1) + ee(5)*alfb(2) + ee(6)*alfb(3)
 ts(4) = ee(7)*alfb(1) + ee(8)*alfb(2) + ee(9)*alfb(3)
 ts(5) = 0.0
 ts(6) = 0.0
 
!     SAVE ECPT(9) TO ECP(22)
 
 DO  iki = 1,14
   phi (iki) = ecpt(8+iki)
 END DO
 GO TO 940
 
!     SET FATAL ERROR FLAG AND ALLOWING ERROR MESSAGES TO ACCUMLATE
 
 910 i = 37
 GO TO 930
 915 i = 126
 GO TO 930
 920 i = 26
 930 CALL mesage (-30,i,idel)
 940 IF (.NOT.pzmat) ksys78 = ksave
 RETURN
END SUBROUTINE stpax1
