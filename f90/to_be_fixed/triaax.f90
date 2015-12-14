SUBROUTINE triaax
     
!     THIS SUBROUTINE COMPUTES THE STIFFNESS AND MASS MATRICES FOR THE
!     ASSYMMETRIC RING WITH A TRIANGULAR CROSS SECTION, TO BE USED BY
!     THE ELEMENT MATRIX GENERATOR.
 
!     SINGLE PRECISION VERSION
 
!     THIS CURRENT VERSION ALLOWS FOR COORDINATE SYSTEM
!     THIS SUBROUTINE USES THE ADDITIONAL ROUTINES DKL, DELTKL
 
!     THE ECPT FOR THE TRIAX ELEMENT IS
 
!     ECPT (01) = ELEMENT ID                              I
!     ECPT (02) = SIL A                                   I
!     ECPT (03) = SIL B                                   I
!     ECPT (04) = SIL C                                   I
!     ECPT (05) = MATERIAL ORIENTATION ANGLE(DEGREES)     R
!     ECPT (07) = MATERIAL ID                             I
!     ECPT (08) TO ECPT(21) = STRESS PHASE ANG.           R
!     ECPT (22) = CORD. SYS. GRID POINT A (NOT USED)      I
!     ECPT (23) = R-CORD OF GRID A                        R
!     ECPT (24) = Z-CORD OF GRID A                        R
!     ECPT (25) = 0.0                                     R
!     ECPT (26) = CORD. SYS. GRID POINT B (NOT USED)      I
!     ECPT (27) = R-CORD OF GRID B                        R
!     ECPT (28) = Z-CORD OF GRID B                        R
!     ECPT (29) = 0.0                                     R
!     ECPT (30) = CORD. SYS. GRID POINT C (NOT USED)      I
!     ECPT (31) = R-CORD OF GRID C                        R
!     ECPT (32) = Z-CORD OF GRID C                        R
!     ECPT (33) = 0.0                                     R
!     ECPT (34) = EL. TEMPERATURE FOR MATERIAL PROP       R
 
!     ANY GROUP OF STATEMENTS PREFACED BY AN IF STATEMENT CONTAINING
!     ...KSYS78 OR LSYS78 ...  INDICATES CODING NECESSARY FOR THIS
!     ELEMENT*S PIEZOELECTRIC CAPABILITY
 
!     KSYS78 = 0   ELASTIC, NON-PIEZOELECTRIC MATERIAL
!     KSYS78 = 1   ELECTRICAL-ELASTIC COUPLED, PIEZOELETRIC MATERIAL
!     KSYS78 = 2   ELASTIC ONLY, PIEZOELECTRIC MATERIAL
!     LSYS78 = .TRUE. IF KSYS78 = 0, OR 2
 
 
 LOGICAL :: nogo,heat,pzmat,lsys78
 INTEGER :: dict(11),elid,estid,isort(3)
 REAL :: gababq(9,9),r(3),z(3),ee(63),teo(45),delint(12),  &
     ecpt(10),delm(12),bmass(9,9),akm(81),akjm(81), amt(9),gababp(3,3)
 DOUBLE PRECISION :: consts
 DIMENSION       iecpt(34),aki(81),akt(16),acurl(117),d(81),ak(81),  &
     akj(144),ics(3),d1(27),d2(9),acurp1(27),acurp2(9), akuph(27),akph2(9),akip(9)
 COMMON /system/ ksystm(77),ksys78
 COMMON /emgprm/ ixtra,dum(14),ismb(3),iprec,nogo,heat,icmbar
 COMMON /emgdic/ dxx,ldict,ngrids,elid,estid
 COMMON /emgest/ idel,igp(3),dgama,dm1,matid,spa(14),ics1,r1,zz1,  &
     zer,ics2,r2,zz2,zer2,ics3,r3,zz3,zer3,tempe
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e(3),anu(3),rho,g(3),alf(3),tzero,gsube,moskp(9), setmat
 COMMON /condad/ consts(5)
 COMMON /matpz / pzout(51)
 
!     COMMON /MATPZ / CE11,CE12,CE13,CE14,CE15,CE16,CE22,CE23,CE24,CE25,
!                     CE26,CE33,CE34,CE35,CE36,CE44,CE45,CE46,CE55,CE56,
!                     CE66,E11,E12,E13,E14,E15,E16,E21,E22,E23,E24,E25,
!                     E26,E31,E32,E33,E34,E35,E36,EPS11,EPS12,EPS13,
!                     EPS22,EPS23,EPS33,RHO,A1,A2,A12,TREF,GE
 
 EQUIVALENCE     (iecpt(1),ecpt(1),idel), (dict(5),dict5),  &
     (z(1),z1), (z(2),z2), (z(3),z3),  &
     (aki(1),gababq(1,1)), (bmass(1,1),akm(1)),  &
     (consts(1),pi), (consts(4),degrad),  &
     (consts(2),twopi), (akip(1),gababp(1,1)),  &
     (acurp1(1),acurl(82)), (acurp2(1),acurl(109))
 DATA            idel2,jax / 0, 4HTRIA/
 
 lsys78 = .false.
 IF (ksys78 == 0 .OR. ksys78 == 2) lsys78 = .true.
 idel1 = idel/1000
 
!     INITALIZE
 
 DO  i = 1,3
   r(i) = ecpt(4*i+19)
   z(i) = ecpt(4*i+20)
   ics(i) = iecpt(4*i+18)
 END DO
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 12
 dict(4) = 15
 ipr     = iprec
 
 IF (r1 <= 0. .OR. r2 <= 0. .OR. r3 <= 0.) GO TO 7770
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = AMIN1(z1,z2,z3)
 z1 = z1 - zmin
 z2 = z2 - zmin
 z3 = z3 - zmin
 
!     FORM TRANSFORMATION MATRIX GABABQ (9X9) FROM FIELD COORDINATES TO
!     GRID POINT DEGREES OF FREEDOM
 
 DO  i = 1,9
   DO  j = 1,9
     gababq(i,j) = 0.
   END DO
 END DO
 
 aa = 1./(r2*z3 + r1*z2 + z1*r3 - z2*r3 - r1*z3 - r2*z1)
 c1 =  aa*(r2*z3 - z2*r3)
 c2 = -aa*(z3 - z2)
 c3 =  aa*(r3 - r2)
 gababq(1,1) = c1
 gababq(1,2) = c2
 gababq(1,3) = c3
 gababq(2,4) = c1
 gababq(2,5) = c2
 gababq(2,6) = c3
 gababq(3,7) = c1
 gababq(3,8) = c2
 gababq(3,9) = c3
 IF (lsys78) GO TO 102
 gababp(1,1) = c1
 gababp(1,2) = c2
 gababp(1,3) = c3
 102 CONTINUE
 c1 = -aa*(r1*z3 - z1*r3)
 c2 =  aa*(z3 - z1)
 c3 = -aa*(r3 - r1)
 gababq(4,1) = c1
 gababq(4,2) = c2
 gababq(4,3) = c3
 gababq(5,4) = c1
 gababq(5,5) = c2
 gababq(5,6) = c3
 gababq(6,7) = c1
 gababq(6,8) = c2
 gababq(6,9) = c3
 IF (lsys78) GO TO 104
 gababp(2,1) = c1
 gababp(2,2) = c2
 gababp(2,3) = c3
 104 CONTINUE
 c1 =  aa*(r1*z2 - z1*r2)
 c2 = -aa*(z2 - z1)
 c3 =  aa*(r2 - r1)
 gababq(7,1) = c1
 gababq(7,2) = c2
 gababq(7,3) = c3
 gababq(8,4) = c1
 gababq(8,5) = c2
 gababq(8,6) = c3
 gababq(9,7) = c1
 gababq(9,8) = c2
 gababq(9,9) = c3
 IF (lsys78) GO TO 110
 gababp(3,1) = c1
 gababp(3,2) = c2
 gababp(3,3) = c3
 110 CONTINUE
 
!     COMPUTE THE INTEGRAL VALUES IN ARRAY DELINT THE ORDER IS INDICATED
!     THE FOLLOWING TABLE
 
!     DELINT(01) = (-1,0)
!     DELINT(02) = (-1,1)
!     DELINT(03) = (-1,2)
!     DELINT(04) = (0, 0)
!     DELINT(05) = (0, 1)
!     DELINT(06) = (1, 0)
 
!     OR FOR THE MASS MATRIX
 
!     DELINT(1)  =  (1,0)
!     DELINT(2)  =  (1,1)
!     DELINT(3)  =  (1,2)
!     DELINT(4)  =  (2,0)
!     DELINT(5)  =  (2,1)
!     DELINT(7)  =  (3,0)
 
 
 IF (ismb(1) == 0) GO TO 180
 ra = (r1 + r2 + r3)/3.0
 za = (z1 + z2 + z3)/3.0
 rh = AMIN1(r1,r2,r3)/10.0
 dr = AMAX1(ABS(r1-r2),ABS(r2-r3),ABS(r3-r1))
 area = (r1*(z2-z3) + r2*(z3-z1) + r3*(z1-z2))/2.0
 
 i1 = 0
 DO  i = 1,2
   ip = i - 2
   DO  j = 1,3
     iq =  j  - 1
     i1 =  i1 + 1
     IF (i1 /= 6) GO TO 120
     ip = 1
     iq = 0
     120 IF (dr > rh) GO TO 130
     delint(i1) = ((ra**ip)*(za**iq))*area
     GO TO 135
     130 delint(i1) = dkls(3,ip,iq,r,z)
     135 delint(i1) = ABS (delint(i1))
   END DO
 END DO
 
!     MASS MATRIX
 
 IF (ismb(2) == 0) GO TO 200
 180 CALL delkls (akj,r,z,0)
 delm (1) = akj(2)
 delm (2) = akj(7)
 delm (3) = akj(8)
 delm (4) = akj(10)
 delm (5) = akj(9)
 delm (7) = akj(12)
 
!     LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3
 
 200 dgamr  = dgama*degrad
 cosg   = COS(dgamr)
 sing   = SIN(dgamr)
 sinth  = sing
 costh  = cosg
 matidc = matid
 matflg = 7
 IF (ksys78 > 0) matflg = 9
 eltemp = tempe
 CALL mat (idel)
 pzmat = .false.
 IF (setmat == 4. .OR. setmat == 5.) pzmat = .true.
 IF (pzmat) GO TO 210
 ksave  = ksys78
 ksys78 = 0
 lsys78 = .true.
 GO TO 220
 210 rho    = pzout(46)
 alf(1) = pzout(47)
 alf(2) = pzout(48)
 alf(3) = pzout(49)
 tzero  = pzout(50)
 gsube  = pzout(51)
 220 CONTINUE
 IF (setmat == 2.) GO TO 7780
 dict5 = g sub e
 IF (ksys78 > 0) GO TO 249
 
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
 del = 1./(1. - vro*vor  - voz*vzo  - vzr*vrz - vro*voz*vzr - vrz*vor*vzo)
 
!     COMPUTE ELASTIC CONSTANTS MATRIX FROM MATERIAL TO ELEMENT AXIS
 
 249 CONTINUE
 DO  i = 1,45
   teo(i) = 0.
 END DO
 
 IF (ksys78 > 0) GO TO 251
 teo(1) = er*(1.  - voz*vzo)*del
 teo(2) = er*(vzr + vzo*vor)*del
 teo(3) = ez*(1.  - vro*vor)*del
 teo(4) = er*(vor + vzr*voz)*del
 teo(5) = et*(vzo + vro*vzr)*del
 teo(6) = et*(1.  - vrz*vzr)*del
 teo(10)= grz
 teo(15)= gor
 teo(21)= gzo
 GO TO 252
 251 CONTINUE
 
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
 IF (ksys78 == 2) GO TO 252
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
 252 CONTINUE
 c2 = cosg*cosg
 c4 = c2  *c2
 s2 = sing*sing
 s4 = s2  *s2
 c2s2 = c2*s2
 c3 = cosg*c2
 s3 = sing*s2
 cs2= cosg*s2
 sc2= sing*c2
 cs = cosg*sing
 
 ee( 1) = teo(1)*c4 + teo(3)*s4 + 2.*c2s2 *(teo(2) + 2.*teo(10))
 ee( 2) = teo(2)*(c4+s4) + c2s2*(teo(1)+teo(3) - 4.*teo(10))
 ee( 3) = teo(4)*c2 + teo(5)*s2
 ee( 4) = cosg*sing*s2*(teo(2) - teo(3) + 2.*teo(10))  &
     + sing*cosg*c2*(teo(1) - teo(2) - 2.*teo(10))
 ee( 7) = ee(2)
 ee( 8) = teo(1)*s4 + 2.*c2s2*(teo(2) + 2.*teo(10))+ teo(3)*c4
 ee( 9) = teo(4)*s2 + teo(5)*c2
 ee(10) = sing*cosg*c2*(teo(2) - teo(3) + 2.*teo(10))  &
     + cosg*sing*s2*(teo(1) - teo(2) - 2.*teo(10))
 ee(13) = ee(3)
 ee(14) = ee(9)
 ee(15) = teo(6)
 ee(16) = sing*cosg*(teo(4)-teo(5))
 ee(19) = ee(4)
 ee(20) = ee(10)
 ee(21) = ee(16)
 ee(22) = c2s2*(teo(1) - 2.*teo(2) + teo(3)) + teo(10)*(c2-s2)**2
 ee(29) = teo(15)*c2 + teo(21)*s2
 ee(30) = sing*cosg*(teo(15)-teo(21))
 ee(35) = ee(30)
 ee(36) = teo(15)*s2 + teo(21)*c2
 
 IF (lsys78) GO TO 254
 
!     PIEZOELECTRIC MATERIAL PROPERTIES IN ELEMENT COORDINATES
 
 ee(37) = c3*teo(22) - s3*teo(26) + cs2*(teo(25)+2.0*teo(32)) -  &
     sc2*(teo(23) + 2.0*teo(31))
 ee(38) = c3*teo(23) + s3*teo(25) + cs2*(teo(26)-2.0*teo(31)) +  &
     sc2*(teo(22) - 2.0*teo(32))
 ee(39) = s2*teo(27) + c2*teo(24) - 2.0*cs*teo(33)
 ee(40) = c3*teo(25) - s3*teo(23) + cs2*(teo(22)-2.0*teo(32)) -  &
     sc2*(teo(26) - 2.0*teo(31))
 ee(41) = c3*teo(26) + s3*teo(22) + cs2*(teo(23)+2.0*teo(31)) +  &
     sc2*(teo(25) + 2.0*teo(32))
 ee(42) = s2*teo(24) + c2*teo(27) + 2.0*cs*teo(33)
 ee(43) = cosg*teo(28) - sing*teo(29)
 ee(44) = cosg*teo(29) + sing*teo(28)
 ee(45) = teo(30)
 ee(46) = c3*teo(31) + s3*teo(32) - cs2*(teo(23)-teo(26)+teo(31)) +  &
     sc2*(-teo(32) - teo(25) + teo(22))
 ee(47) = c3*teo(32) - s3*teo(31) - cs2*(teo(25)-teo(22)+teo(32)) +  &
     sc2*(teo(23) + teo(31) - teo(26))
 ee(48) = (c2-s2)*teo(33) + cs*(teo(24) - teo(27))
 ee(49) = c2*teo(34) + s2*teo(38) - cs*(teo(35) + teo(37))
 ee(50) = c2*teo(35) - s2*teo(37) + cs*(teo(34) - teo(38))
 ee(51) = cosg*teo(36) - sing*teo(39)
 ee(52) = c2*teo(37) - s2*teo(35) - cs*(teo(38) - teo(34))
 ee(53) = c2*teo(38) + s2*teo(34) + cs*(teo(35) + teo(37))
 ee(54) = cosg*teo(39) + sing*teo(36)
 
!     DIELECTRIC MATERIAL PROPERTIES IN ELEMENT COORDINTES
 
 ee(55) = s2*teo(43) - 2.0*cs*teo(41) + c2*teo(40)
 ee(56) = (c2-s2)*teo(41) - cs*(teo(43)-teo(40))
 ee(57) =-sing*teo(44) + cosg*teo(42)
 ee(59) = c2*teo(43) + 2.0*cs*teo(41) + s2*teo(40)
 ee(60) = cosg*teo(44) + sing*teo(42)
 ee(63) = teo(45)
 254   CONTINUE
 
!     COMPUTE HARMONIC COEFFICIENT
 
 mjho = MOD(iecpt(1),1000) - 1
 ajho = mjho
 ajjho= ajho*ajho
 rhod = rho *pi
 IF (ajho   == 0.) rhod = 2.*rhod
 IF (ismb(1) == 0) GO TO 300
 
!     FORM THE ELEMENT STIFFNESS MATRIX IN FIELD SYSTEM
 
 acurl( 1) = (ee(15) + ajjho*ee(29))*delint(1)
 acurl( 2) = (ee(03) + ee(15) + ajjho*ee(29))*delint(4)
 acurl( 3) = (ee(15) + ajjho*ee(29))*delint(2) + ee(16)*delint(4)
 acurl( 4) = (ee(15) + ee(29))*ajho*delint(1)
 acurl(05) = ee(15)*ajho*delint(4)
 acurl(06) = (ee(15) + ee(29))*ajho*delint(2) - ee(30)*ajho*delint(4)
 acurl(07) = ajjho*delint(1)*ee(35)
 acurl(08) = (ee(16) + ajjho*ee(35))*delint(4)
 acurl(09) = ee(9)*delint(4) + ajjho*delint(2)*ee(35)
 acurl(11) = (ee(1) + 2.*ee(3) + ee(15) + ajjho*ee(29))*delint(6)
 acurl(12) = (ee(3) + ee(15) + ajjho*ee(29))*delint(5)  &
     + (ee(4) + ee (16))*delint(6)
 acurl(13) = (ee(3) + ee(15) + ee(29))*ajho*delint(4)
 acurl(14) = (ee(3) + ee(15))*delint(6)*ajho
 acurl(15) = (ee(3) + ee(15) + ee(29))*ajho*delint(5) - ajho*ee(30)*delint(6)
 acurl(16) = ajjho*delint(4)*ee(35)
 acurl(17) = (ee(4) + ee(16) + ajjho*ee(35))*delint(6)
 acurl(18) = (ee(2) + ee(9))*delint(6) + ajjho*delint(5)*ee(35)
 acurl(21) = (ee(15) + ajjho*ee(29))*delint(3) + ee(22)  &
     * delint(6) +  2.*ee(16)*delint(5)
 acurl(22) = (ee(15) + ee(29))*ajho*delint (2) + ajho * delint(4)*ee(16)
 acurl(23) = ee(15)*ajho*delint(5) + ajho*delint(6)*ee(16)
 acurl(24) = (ee(15) + ee(29))*ajho*delint(3) + (ee(16) - ee(30))  &
     * ajho*delint(5)
 acurl(25) = ajjho*delint(2)*ee(35)
 acurl(26) = ee(22)*delint(6) + (ee(21) + ajjho*ee(35))*delint(5)
 acurl(27) = ee(9)*delint(5) + ee(10)*delint(6) + ajjho * delint(3)*ee(35)
 acurl(31) = (ee(29) + ajjho*ee(15))*delint(1)
 acurl(32) = ee(15)*ajjho*delint(4)
 acurl(33) = (ee(29) + ajjho*ee(15))*delint(2) - ee(30)*delint(4)
 acurl(34) = ajho*delint(1)*ee(35)
 acurl(35) = ajho*(ee(16) + ee(35))*delint(4)
 acurl(36) = ee(9)*ajho*delint(4) + ajho*delint(2)*ee(35)
 acurl(41) = ajjho*delint(6)*ee(15)
 acurl(42) = ee(15)*ajjho*delint(5)
 acurl(43) = 0.
 acurl(44) = ajho*delint(6)*ee(16)
 acurl(45) = ee(9)*ajho*delint(6)
 acurl(51) = (ee(29) + ajjho*ee(15))*delint(3) + ee(36)  &
     * delint(6) - 2.*ee(35)*delint(5)
 acurl(52) = ajho*(delint(2)*ee(30) - delint(4)*ee(36))
 acurl(53) = -ee(36)*ajho*delint(6) + ajho*(ee(16) + ee(35)) * delint(5)
 acurl(54) = (ee(9) - ee(36))*ajho*delint(5) + ajho * delint(3)*ee(35)
 acurl(61) = ee(36)*ajjho*delint(1)
 acurl(62) = ee(36)*ajjho*delint(4)
 acurl(63) = (ee(36))*ajjho*delint(2)
 acurl(71) = (ee(22) + ajjho*ee(36))*delint(6)
 acurl(72) = ee(36)*ajjho*delint(5) + ee(20)*delint(6)
 acurl(81) = ee(36)*ajjho*delint(3) + ee(8)*delint(6)
 IF (lsys78) GO TO 256
 acurl(82) =-(ee(45) - ajho*ee(51))*ajho*delint(1)
 acurl(83) = (ee(43) - ajho*ee(45) - ajho*ee(49) + ajjho*ee(51)) * delint(4)
 acurl(84) = (ee(44) - ajho*ee(50))*delint(4) - (ee(45)  &
     - ajho*ee(51))*ajho*delint(2)
 acurl(85) =-(ee(39) + ee(45) - ajho*ee(51))*ajho*delint(4)
 acurl(86) = (ee(37) + ee(43) - (ee(39) + ee(45) + ee(49)  &
     - ajho*ee(51))*ajho)*delint(6)
 acurl(87) = (ee(38) + ee(44) - ajho*ee(50))*delint(6) - (ee(39)  &
     + ee(45) - ajho*ee(51))*ajho*delint(5)
 acurl(88) =-(ee(45) - ajho*ee(51))*ajho*delint(2) - ee(48)*ajho * delint(4)
 acurl(89) = (ee(43) - ajho*ee(45) - ajho*ee(49) + ajjho*ee(51))  &
     * delint(5) + (ee(46) - ee(48)*ajho)*delint(6)
 acurl(90) = (ee(44) - ajho*ee(48) - ajho*ee(50))*delint(5)  &
     + ee(47)*delint(6) - (ee(45)-ajho*ee(51))*ajho*delint(3)
 acurl(91) =-(ee(45)*ajho - ee(51))*ajho*delint(1)
 acurl(92) = (ajho*ee(43) - ajjho*ee(45) - ee(49) + ajho*ee(51)) * delint(4)
 acurl(93) = (ee(44)*ajho - ee(50))*delint(4) - (ee(45)*ajho  &
     - ee(51))*ajho*delint(2)
 acurl(94) =-ee(45)*ajjho*delint(4)
 acurl(95) = (ee(43) - ajho*ee(45))*ajho*delint(6)
 acurl(96) = ee(44)*ajho*delint(6) - ee(45)*ajjho*delint(5)
 acurl(97) =-(ee(45)*ajho - ee(51))*ajho*delint(2) - ee(54)*ajho * delint(4)
 acurl(98) = (ee(43)*ajho - ajjho*ee(45) - ee(49) + ee(51)*ajho)  &
     * delint(5) + (ee(52) - ajho*ee(54))*delint(6)
 acurl(99) = (ee(44)*ajho - ee(50) - ee(54)*ajho)*delint(5)  &
     + ee(53)*delint(6) - (ee(45)*ajho-ee(51))*ajho*delint(3)
 acurl(100)= ee(54)*ajjho*delint(1)
 acurl(101)=-(ee(52) - ee(54)*ajho)*ajho*delint(4)
 acurl(102)=-(ee(53)*delint(4) - ee(54)*ajho*delint(2))*ajho
 acurl(103)=-(ee(48) - ee(54)*ajho)*ajho*delint(4)
 acurl(104)= (ee(46) - ee(48)*ajho - ee(52)*ajho+ee(54)*ajjho) *  delint(6)
 acurl(105)= (ee(47) - ee(53)*ajho)*delint(6) - (ee(48) - ee(54)  &
     * ajho)*ajho*delint(5)
 acurl(106)= ee(54)*ajjho*delint(2) - ee(42)*ajho*delint(4)
 acurl(107)= (ee(40) - ee(42)*ajho)*delint(6) - (ee(52) - ee(54)  &
     * ajho)*ajho*delint(5)
 acurl(108)= ee(41)*delint(6) + (-ee(42) - ee(53))*ajho*delint(5)  &
     + ee(54)*ajjho*delint(3)
 acurl(109)= ee(63)*ajjho*delint(1)
 acurl(110)= (-ee(57) + ee(63)*ajho)*ajho*delint(4)
 acurl(111)=-ee(60)*ajho*delint(4) + ee(63)*ajjho*delint(2)
 acurl(112)= acurl(110)
 acurl(113)= (ee(55) - 2.0*ee(57)*ajho + ee(63)*ajjho)*delint(6)
 acurl(114)= (ee(56) - ee(60)*ajho)*delint(6) + (-ee(57) + ee(63)  &
     * ajho)*ajho*delint(5)
 acurl(115)= acurl(111)
 acurl(116)= acurl(114)
 acurl(117)= ee(59)*delint(6) - 2.0*ee(60)*ajho*delint(5) + ee(63)  &
     * ajjho*delint(3)
 256 CONTINUE
 
!     EXPAND ACURL INTO (9X9)
 
 DO  ib = 2,9
   ic = 10*ib - 19
   i  = ic
   DO  j = ib,9
     ic = ic + 9
     i  = i  + 1
     acurl(ic) = acurl(i)
   END DO
 END DO
 dgama = pi
 IF (ajho == 0.) dgama = twopi
 DO  i = 1,81
   acurl(i) = acurl(i)*dgama
 END DO
 IF (lsys78) GO TO 300
 
 DO  i = 82,117
   acurl(i) = acurl(i)*dgama
 END DO
 
 300 IF (ismb(2) == 0) GO TO 400
 IF (icmbar < 0) GO TO 350
 
!     CONSISTENT MASS IN FIELD COORDINATES
 
 DO  i = 1,9
   DO  j = 1,9
     bmass(i,j) = 0.
   END DO
 END DO
 bmass(1,1) = rhod*delm(1)
 bmass(1,2) = rhod*delm(4)
 bmass(1,3) = rhod*delm(2)
 bmass(2,1) = rhod*delm(4)
 bmass(2,2) = rhod*delm(7)
 bmass(2,3) = rhod*delm(5)
 bmass(3,1) = rhod*delm(2)
 bmass(3,2) = rhod*delm(5)
 bmass(3,3) = rhod*delm(3)
 bmass(4,4) = rhod*delm(1)
 bmass(4,5) = rhod*delm(4)
 bmass(4,6) = rhod*delm(2)
 bmass(5,4) = rhod*delm(4)
 bmass(5,5) = rhod*delm(7)
 bmass(5,6) = rhod*delm(5)
 bmass(6,4) = rhod*delm(2)
 bmass(6,5) = rhod*delm(5)
 bmass(6,6) = rhod*delm(3)
 bmass(7,7) = rhod*delm(1)
 bmass(7,8) = rhod*delm(4)
 bmass(7,9) = rhod*delm(2)
 bmass(8,7) = rhod*delm(4)
 bmass(8,8) = rhod*delm(7)
 bmass(8,9) = rhod*delm(5)
 bmass(9,7) = rhod*delm(2)
 bmass(9,8) = rhod*delm(5)
 bmass(9,9) = rhod*delm(3)
 GO TO 400
 350 area  = (r1*(z2-z3) + r2*(z3-z1) + r3*(z1-z2))/2.
 convm =  rhod*(r1 + r2 + r3)/3.*area
 
!     TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM FIELD SYSTEM
!     TO GRID POINT DEGREES OF FREEDOM
 
 400 IF (ismb(1) == 0) GO TO 410
 CALL gmmats (aki,9,9,1,acurl,9,9,0,d)
 CALL gmmats (d,  9,9,0,aki,  9,9,0,ak)
 IF (lsys78) GO TO 405
 CALL gmmats (aki,9,9,1,acurp1,9,3,0,d1)
 CALL gmmats (d1, 9,3,0,akip,  3,3,0,akuph)
 CALL gmmats (akip,3,3,1,acurp2,3,3,0,d2)
 CALL gmmats (d2,  3,3,0,akip,  3,3,0,akph2)
 405 CONTINUE
 
 IF (ismb(2) == 0 .OR. icmbar < 0) GO TO 450
 410 CALL gmmats (aki,9,9,1,bmass,9,9,0,d)
 CALL gmmats (d,  9,9,0,aki,  9,9,0,akm)
 
 450 DO  i = 1,81
   akj(i)  = 0.
   akjm(i) = 0.
 END DO
 DO  i = 82,117
   akj(i) = 0.0
 END DO
 
 GO TO 480
 
!     COORDINATE SYSTEMS POSSIBLE WITH RINGAX THRU CODE BELOW
! **  IF FOLLOWING CODE IS IMPLEMENTED MUST BE MODIFIED FOR PIEZOELECTRI
 
!     DO 470 I = 1,3
!     IF (ICS(I) .EQ. 0) GO TO 470
!     K = 9*(I-1) + 1
!     CALL TRANSS (ICS(I),D(K))
! 470 CONTINUE
 
!     CREATE AN ARRAY OF SORTED GRID POINTS
 
 480 DO  i = 1,3
   isort(i) = igp(i)
 END DO
 i = -3
 484 j = 0
 DO  k = 1,3
   IF (isort(k) < j) CYCLE
   j = isort(k)
   l = k
 END DO
 isort(l) = i
 i = i + 1
 IF (i < 0) GO TO 484
 DO  i = 1,3
   isort(i) = -isort(i)
 END DO
 
!     TRANSFORM 3 X 3 TO 6 X 6 FOR COORD SYSTEM TRANSFORMATIONS
 
 DO  isil = 1,3
   ipp = isort(isil)
   
   ir1 = 3*(isil-1) + 1
   DO  ii = 1,3
     i   = isort(ii)
     ic1 = 3*(ii-1) + 1
     irc = (ir1 -1)*9 + ic1
     DO  j = 1,3
       j1   = (j-1)*4 + 1
       ircc = irc + (j-1)*9
       IF (ismb(1) == 0) GO TO 495
       akt(j1  ) = ak(ircc  )
       akt(j1+1) = ak(ircc+1)
       akt(j1+2) = ak(ircc+2)
       IF (lsys78) GO TO 492
       m = ircc/3 + 1
       n = (m-1)/9 + 1 + (ii-1)*9 + (j-1)*3
       akt(j1+3) = akuph(m)
       akt(j1+15-j*3) = akuph(n)
       akt(16) = akph2(ir1+ii-1)
       492 CONTINUE
       
       495 IF (ismb(2) == 0 .OR. icmbar < 1) CYCLE
       j1 = (j-1)*3 + 1
       amt(j1  ) =akm(ircc  )
       amt(j1+1) =akm(ircc+1)
       amt(j1+2) =akm(ircc+2)
     END DO
     
     GO TO 540
     
!    ABOVE GO TO MAKES CST CODE BELOW INTO DEAD CODE
!    COORDINATE SYSTEM TRANSFORMATION CODE
! ** IF FOLLOWING CODE IS IMPLEMENTED MUST BE MODIFIED FOR PIEZOELECTRIC
     
!     IF (ICS(IPP) .EQ. 0) GO TO 520
!     IAA = 9*(IPP-1) + 1
!     CALL GMMATS (D(IAA),3,3,1,AKT(1),3,3,0,D(28))
!     CALL GMMATS (D(IAA),3,3,1,AMT(1),3,3,0,D(37))
!     DO 510 J = 1,9
!     AKT(J) = D(J+27)
! 510 AKM(J) = D(J+36)
     
! 520 IF (ICS(I) .EQ. 0) GO TO 540
!     IAI = 9*(I-1) + 1
!     CALL GMMATS (AKT(1),3,3,0,D(IAI),3,3,0,D(28))
!     CALL GMMATS (AMT(1),3,3,0,D(IAI),3,3,0,D(37))
!     DO 530 J = 1,9
!     AKT(J) = D(J+27)
! 530 AMT(J) = D(J+36)
     
!     NOW INSERT  AKT AND AMT INTO THE OVERALL STIFFNESS MATRICES
!     ACCORDING TO INCREASING SIL VALUE
     
     540 DO  ij = 1,3
       DO  jj = 1,3
         ki = (ij-1)*3 + jj
         iout = (ipp-1)*27 + (i-1)*3 + (ij-1)*9 + jj
         akjm(iout)= amt(ki)
       END DO
     END DO
     DO  ij = 1,4
       DO  jj = 1,4
         ki = (ij-1)*4 + jj
         iout = (ipp-1)*48 + (i-1)*4 + (ij-1)*12 + jj
         akj(iout) = akt(ki)
       END DO
     END DO
   END DO
 END DO
 
!     NOW OUTPUT THE MATRIX VIA EMG OUT
 
 dict(2) = 1
 IF (ismb(1) == 0) GO TO 650
 CALL emgout (akj,akj,144,1,dict,1,ipr)
 650 IF (ismb(2) == 0 .AND. .NOT.pzmat) ksys78 = ksave
 IF (ismb(2) == 0) RETURN
 dict(3) = 9
 dict(4) = 7
 IF (icmbar < 0) GO TO 670
 CALL emgout (akjm,akjm,81,1,dict,2,ipr)
 GO TO 700
 
!     GENERATE LUMPED MASS MATRIX HERE
 
 670 DO  i = 1,9
   akjm(i) = convm/3.0
 END DO
 dict(2) = 2
 CALL emgout (akjm,akjm,9,1,dict,2,ipr)
 700 IF (.NOT.pzmat) ksys78 = ksave
 RETURN
 
!     ERROR EXITS
 
 7770 i = 37
 7777 IF (idel1 == idel) GO TO 7778
 idel2  = idel1
 ics(1) = idel1
 ics(2) = jax
 CALL mesage (30,i,ics)
 7778 nogo = .true.
 GO TO 700
 7780 i = 126
 GO TO 7777
END SUBROUTINE triaax
