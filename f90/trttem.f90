SUBROUTINE trttem ( ti, pg)
!****
! THIS ROUTINE COMPUTES THE THERMAL LOAD FOR THE ASSYMMETRIC RING ELE
! WITH A TRIANGULAR CROSS SECTION
!****
! ECPT (01) = ELEMENT ID                         I
! ECPT (02) = SIL A                              I
! ECPT (03) = SIL B                              I
! ECPT (04) = SIL C                              I
! ECPT (05) = MATERIAL ORIENTATION ANGLE(DEGREES)R
! ECPT (07) = MATERIAL ID                        I
! ECPT (08) TO  ECPT (21) = PHI                  R
! ECPT (22) = COOR. SYS. FOR GRID POINT A        I
! ECPT (23) = R-CORD OF GRID A                   R
! ECPT (24) = Z-CORD OF GRID A                   R
! ECPT (25) = 0.0                                R
! ECPT (26) = CORD. SYS. GRID POINT B (NOT USED) I
! ECPT (27) = R-CORD OF GRID B                   R
! ECPT (28) = Z-CORD OF GRID B                   R
! ECPT (29) = 0.0                                R
! ECPT (30) = CORD. SYS. GRID POINT C (NOT USED) I
! ECPT (31) = R-CORD OF GRID C                   R
! ECPT (32) = Z-CORD OF GRID C                   R
! ECPT (33) = 0.0                                R
! ECPT (34) = EL. TEMPERATURE FOR MATERIAL       R
 
 
 REAL, INTENT(IN)                         :: ti(3)
 REAL, INTENT(OUT)                        :: pg(1)
 DIMENSION r(3), z(3), gababq(9,9),delint(12),teo(21)  &
     ,         dtt(9),igp(3),iecpt(34),ics(3), fij(9) ,         d(3)  ,tl(9)
 
!  . ECPT COMMON BLOCK
 COMMON    /trimex/  ecpt(34)
 
!  . MATERIAL INPUT AND OUTPUT...
 COMMON    /matin/ matidc            ,matflg  &
     ,                  eltemp            ,stress  &
     ,                  sinth             ,costh
 
 COMMON    /matout/ e(3)              ,anu(3)  &
     ,                  rho               ,g(3)  &
     ,                  alf(3)            ,tzero      ,gsube  &
     ,                   moskp(9)            ,setmat
 COMMON    /condas/  consts(5)
 EQUIVALENCE (iecpt(1), ecpt(1)),  (z(1), z1),  (z(2), z2)  &
     ,             ( z(3), z3) ,        ( r(1), r1),  ( r(2), r2),  (r(3), r3)
 EQUIVALENCE (consts(1),pi), (consts(4),degrad)
 
! START EXECUTION
 
! STORE ECPT PARAMETERS IN LOCAL VARIABLES
 idel = iecpt(1)
 igp(1) = iecpt(2)
 igp(2) = iecpt(3)
 igp(3) = iecpt(4)
 matid  = iecpt(07)
 ics(1) = iecpt(22)
 ics(2) = iecpt(26)
 ics(3) = iecpt(30)
 r(1)   = ecpt(23)
 r(2)   = ecpt(27)
 r(3)   = ecpt(31)
 z(2)   = ecpt(28)
 d(2)   = ecpt(29)
 z(1)   = ecpt(24)
 d(1)   = ecpt(25)
 z(3)   = ecpt(32)
 d(3)   = ecpt(33)
 dgama  = ecpt(05)
 tempe  = ecpt(34)
 
! COMPUTE THE ELEMENT COORDINATES
 zmin = AMIN1(z1, z2, z3)
 z1 = z1 - zmin
 z2 = z2 - zmin
 z3 = z3 - zmin
 
! FORM THE TRANSFORMATION MATRIX GABABQ (9X9) FROM FIELD COORDINATES TO
! GRID POINT DEGREES OF FREEDOM
 DO  i = 1,9
   DO  j = 1,9
     gababq (i,j) = 0.000
   END DO
 END DO
 aa = r2 * z3 + r1 * z2 + z1 * r3 - z2 * r3 - r1 * z3 - r2 * z1
 aa = 1.0E0 / aa
 c1 = aa * ( r2 * z3 - z2 * r3)
 c2 = - aa * ( z3 - z2 )
 c3 = aa * ( r3 - r2 )
 gababq (1,1) = c1
 gababq (2,4) = c1
 gababq (3,7) = c1
 gababq (1,2) = c2
 gababq (2,5) = c2
 gababq (3,8) = c2
 gababq (1,3) = c3
 gababq(2,6) = c3
 gababq(3,9) = c3
 c1 = -aa * ( r1 * z3 - z1 * r3)
 c2 =  aa * ( z3 - z1 )
 c3 = -aa * ( r3 - r1 )
 gababq(4,1) = c1
 gababq(4,2) = c2
 gababq(4,3) = c3
 gababq(5,4) = c1
 gababq(5,5) = c2
 gababq(5,6) = c3
 gababq(6,7) = c1
 gababq(6,8) = c2
 gababq(6,9) = c3
 c1 = aa * ( r1 * z2 - z1 * r2)
 c2 = -aa * ( z2 - z1)
 c3 = aa * ( r2 - r1 )
 gababq(7,1) = c1
 gababq(7,2) = c2
 gababq(7,3) = c3
 gababq(8,4) = c1
 gababq(8,5) = c2
 gababq(8,6) = c3
 gababq(9,7) = c1
 gababq(9,8) = c2
 gababq(9,9) = c3
 
! LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3
 dgamr = dgama * degrad
 cosg = COS (dgamr)
 sing = SIN (dgamr)
 costh = cosg
 sinth = sing
 matidc = matid
 matflg = 7
 eltemp = tempe
 CALL mat (idel)
 IF (setmat == 2.0) GO TO 910
 
!  . SET MATERIAL PROPERTIES IN LOCAL VARIABLES...
 er = e(1)
 et = e(2)
 ez = e(3)
 vro = anu(1)
 voz = anu(2)
 vzr = anu(3)
 gor = g(1)
 gzo = g(2)
 grz = g(3)
 vor = vro * et / er
 vzo = voz * ez / et
 vrz = vzr * er / ez
 del = 1.0E0 / (1.0E0 - vro * vor - voz * vzo - vzr * vrz  &
     - vro * voz * vzr - vrz * vor * vzo )
 
! COMPUTE ELASTIC CONSTANTS MATRIX FROM MATERIAL TO ELEMENT AXIS
 DO  i = 1,21
   teo (i) = 0.0E0
 END DO
 teo (1) = er * ( 1.0E0 - voz * vzo) * del
 teo (2) = er * ( vzr + vzo * vor ) * del
 teo (3) = ez * ( 1.0E0 - vro * vor ) * del
 teo (4) = er * ( vor + vzr * voz) * del
 teo (5) = et * (vzo + vro * vzr ) * del
 teo (6) = et * ( 1.0E0 - vrz * vzr ) * del
 teo (10) = grz
 teo (15) = gor
 teo (21) = gzo
 c2 = cosg * cosg
 c4 = c2 * c2
 s2 = sing * sing
 s4 = s2 * s2
 c2s2 = c2 * s2
 ee01 = teo(1)*c4 + teo(3)*s4 + (teo(2) + 2.0E0*teo(10)) * 2.0E0*c2s2
 ee02 = teo(2) * (s4+c4) + (teo(1) + teo(3) - 4.0E0*teo(10))*c2s2
 ee03 = teo(4)*c2 + teo(5)*s2
 ee04 = sing*cosg * (teo(1)*c2 - teo(3)*s2 + (teo(2) + 2.0E0  &
     *teo(10)) * (s2-c2))
 ee08 = teo(1)*s4 + teo(3)*c4 + (2.0E0*teo(2) + 4.0E0*teo(10)) * c2s2
 ee09 = teo(4)*s2 + teo(5)*c2
 ee10 = sing*cosg * (teo(1)*s2 - teo(3)*c2 + (teo(2) + 2.0E0  &
     * teo(10)) * (c2 - s2))
 ee15 = teo(6)
 ee16 = sing*cosg * (teo(4) - teo(5))
 
! COMPUTE HARMONIC COEFFICIENT
 ajho = iecpt(1) - (iecpt(1) /1000) * 1000 - 1
 
!  . CALCULATE THE INTEGRAL VALUES IN DELINT...
 
!         DELINT(4) = 0,0
!         DELINT(5) = 0,1
!         DELINT(6) = 1,0
 
 delint(4) = ais (3,0,0,r,z)
 delint(5) = ais (3,0,1,r,z)
 delint(6) = ais (3,1,0,r,z)
 
 t1 = ee01*alf(1) + ee02*alf(3) + ee03*alf(2)
 t2 = ee02*alf(1) + ee08*alf(3) + ee09*alf(2)
 t3 = ee03*alf(1) + ee09*alf(3) + ee15*alf(2)
 t4 = ee04*alf(1) + ee10*alf(3) + ee16*alf(2)
! GENERATE DTT MATRIX
 dtt (1) = delint (4) * t3
 dtt (2) = delint (6) * ( t1 + t3 )
 dtt (3) = delint (5) * t3 + delint(6) * t4
 dtt (4) = ajho * dtt (1)
 dtt (5) = ajho * delint (6)*t3
 dtt (6) = ajho * delint (5) * t3
 dtt (7) = 0.0
 dtt (8) = delint (6) * t4
 dtt (9) = delint (6) * t2
 
! TRANSFORM THE THERMAL LOAD TO GRID POINT DEGREES OF FREEDOM
 CALL gmmats (gababq, 9, 9, 1, dtt, 9, 1, 0, fij)
 t = tzero
 IF (ajho > 0.0) t = 0.0
 t = ((ti(1) + ti(2) + ti(3))/3.0E0 - t) * pi
 IF ( ajho == 0.0 ) t = t * 2.0E0
 DO  i = 1, 9
   tl(i) = t * fij(i)
 END DO
 
!**** THE FOLLOWING CODE REMOVED.  CORD.SYS. NOT POSSIBLE WITH RINGAX **
!.959 FIJ(I) = T*FIJ(I)
!.
!. LOCATE THE TRANSFORMATION MATRICES FOR THE THREE GRID POINTS
!.    DO 750 I=1,81
!.750 AKI(I) = 0.0
!.     DO 800 I = 1,3
!.    CALL GBTRAN(ICS(I),IECPT(4*I+22),DTT(1))  **R,TH,Z NEEDED**
!.    K=30*(I-1) + 1
!.     DO 800 J=1,3
!.    KK = K+9*(J-1)
!.    JJ=3*(J-1)+1
!.    AKI(KK) = DTT(JJ)
!.    AKI(KK+1) = DTT(JJ+1)
!.    AKI(KK+2) = DTT(JJ+2)
!.800 CONTINUE
!.
!. TRANSFORM THE THERMAL LOAD FROM BASIC TO LOCAL COORD...
!.    CALL GMMATS (AKI(1),9,9,1, FIJ(1),9,1,0, TL(1))
 
! ADD THE ELEMENT THERMAL LOAD TO THE STRUCTURE THERMAL LOAD
 k = 0
 DO  i = 1, 3
   l = igp(i) - 1
   DO  j = 1,3
     k = k + 1
     l = l + 1
     pg(l) = pg(l) +  tl(k)
   END DO
 END DO
 GO TO 920
 910 CALL mesage (-30,37,ecpt(1))
 920 RETURN
END SUBROUTINE trttem
