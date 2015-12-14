SUBROUTINE tpztem (ti,pg)
     
!     THIS ROUTINE COMPUTES THE THERMAL LOAD FOR THE AXI-SYMMETRIC
!     TRAPEZOIDAL CROSS SECTION RING
 
!     ECPT COMMON BLOCK IS,
 
!     ECPT( 1) = ELEMENT ID                                I
!     ECPT( 2) = SIL A                                     I
!     ECPT( 3) = SIL B                                     I
!     ECPT( 4) = SIL C                                     I
!     ECPT( 5) = SIL D
!     ECPT( 6) = MATERIAL ORIENTATION ANGLE(DEGREES)       R
!     ECPT( 8) = MATERIAL ID                               I
!     ECPT( 9) TO ECPT (22) FOR PHI
!     ECPT(23) = COOR. SYS. FOR GRID POINT A               I
!     ECPT(24) = X-COOR. OF GRID POINT A (IN BASIC COOR)   R
!     ECPT(25) = Z-COOR. OF GRID POINT A (IN BASIC COOR)   R
!     ECPT(26) = 0.0
!     ECPT(27) = COOR. SYS. FOR GRID POINT B
!     ECPT(28) = X-COOR. OF GRID POINT B (IN BASIC COOR)   R
!     ECPT(29) = Z-COOR. OF GRID POINT B (IN BASIC COOR)   R
!     ECPT(30) = 0.0
!     ECPT(31) = COOR. SYS. FOR GRID POINT C               I
!     ECPT(32) = X-COOR. FOR GRID POINT C                  R
!     ECPT(33) = Z-COOR. FOR GRID POINT C                  R
!     ECPT(34) = 0.0
!     ECPT(35) = COOR. SYS. FOR GRID POINT D               I
!     ECPT(36) = X-COOR FOR GRID POINT D                   R
!     ECPT(37) = Z-COOR FOR GRID POINT D                   R
!     ECPT(38) = 0.0
!     ECPT(39) = EL. TEMPERATURE FOR MATERIAL PROP         R
 
 
 REAL, INTENT(OUT)                        :: ti(4)
 REAL, INTENT(OUT)                        :: pg(1)
 INTEGER :: sp(36)
 DIMENSION  r(4),z(4),gababq(12,12),delint(15),  &
     d(144),teo(21),htn(12,4),igp(4),iecpt(39),ics(4), h(4,4),aki(144),tl(12)
 COMMON /trimex/ ecpt(39)
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e(3),anu(3),rho,g(3),alf(3),tzero,gsube,moskp(9), setmat
 COMMON /condas/ consts(5)
 COMMON /system/ ibuf,iout
 EQUIVALENCE     (iecpt(1),ecpt(1)),(z(1),z1),(z(2),z2),  &
     (z(3),z3),(r(1),r1),(r(2),r2),(r(3),r3),  &
     (r(4),r4),(z(4),z4),(gababq(1,1),aki(1)), (consts(1),pi),(consts(4),degrad)
 DATA    idel2 , jax / 0, 4HAX    /
 
!     START EXECUTION
 
!     STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel   = iecpt( 1)
 igp(1) = iecpt( 2)
 igp(2) = iecpt( 3)
 igp(3) = iecpt( 4)
 igp(4) = iecpt( 5)
 matid  = iecpt( 8)
 ics(1) = iecpt(23)
 ics(3) = iecpt(31)
 ics(2) = iecpt(27)
 r(1)   = ecpt (24)
 r(2)   = ecpt (28)
 r(3)   = ecpt (32)
 ics(4) = iecpt(35)
 z(1)   = ecpt (25)
 d(1)   = ecpt (26)
 z(2)   = ecpt (29)
 d(2)   = ecpt (30)
 z(3)   = ecpt (33)
 d(3)   = ecpt (34)
 z(4)   = ecpt (37)
 d(4)   = ecpt (38)
 r(4)   = ecpt (36)
 tempe  = ecpt (39)
 dgama  = ecpt ( 6)
 idel1  = idel/1000
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = AMIN1(z1,z2,z3,z4)
 z1   = z1 - zmin
 z2   = z2 - zmin
 z3   = z3 - zmin
 z4   = z4 - zmin
 
!     FATAL IF RATIO OF RADII IS TOO LARGE FOR GUASS QUADRATURE
 
 rmin = AMIN1(r1,r2,r3,r4)
 rmax = AMAX1(r1,r2,r3,r4)
 IF (rmin == 0.0) GO TO 206
 IF (rmax/rmin > 10.) GO TO 915
 
 206 CONTINUE
 IF (ABS(z1-z2) > .001) GO TO 910
 IF (ABS(z3-z4) > .001) GO TO 910
 d(5) = (r1+r4)/2.0
 d(6) = (r2+r3)/2.0
 IF (d(5) == 0.0) GO TO 210
 IF (ABS((r1-r4)/d(5)) > .005) GO TO 210
 r1 = d(5)
 r4 = d(5)
 210 CONTINUE
 IF (d(6) == 0.0) GO TO 220
 IF (ABS((r2-r3)/d(6)) > .005) GO TO 220
 r2 = d(6)
 r3 = d(6)
 220 CONTINUE
 
!     FORM THE TRANSFORMMATION MATRIX(12X12) FROM FIELD COOR, TO GRID
!     POINT DEGREES OF FREEDOM
 
 DO  i = 1, 144
   gababq( i, 1) = 0.000
 END DO
 gababq( 1, 1) = 1.000
 gababq( 2, 1) = r1
 gababq( 3, 1) = z1
 gababq( 4, 1) = r1*z1
 gababq( 5, 2) = 1.000
 gababq( 6, 2) = r1
 gababq( 7, 2) = z1
 gababq( 8, 2) = gababq(4,1)
 gababq( 9, 3) = 1.000
 gababq(10, 3) = r1
 gababq(11, 3) = z1
 gababq(12, 3) = gababq(4,1)
 gababq( 1, 4) = 1.000
 gababq( 2, 4) = r2
 gababq( 3, 4) = z2
 gababq( 4, 4) = r2*z2
 gababq( 5, 5) = 1.000
 gababq( 6, 5) = r2
 gababq( 7, 5) = z2
 gababq( 8, 5) = gababq(4,4)
 gababq( 9, 6) = 1.000
 gababq(10, 6) = r2
 gababq(11, 6) = z2
 gababq(12, 6) = gababq(4,4)
 gababq( 1, 7) = 1.000
 gababq( 2, 7) = r3
 gababq( 3, 7) = z3
 gababq( 4, 7) = r3*z3
 gababq( 5, 8) = 1.000
 gababq( 6, 8) = r3
 gababq( 7, 8) = z3
 gababq( 8, 8) = gababq(4,7)
 gababq( 9, 9) = 1.000
 gababq(10, 9) = r3
 gababq(11, 9) = z3
 gababq(12, 9) = gababq(4,7)
 gababq( 1,10) = 1.000
 gababq( 2,10) = r4
 gababq( 3,10) = z4
 gababq( 4,10) = r4*z4
 gababq( 5,11) = 1.000
 gababq( 6,11) = r4
 gababq( 7,11) = z4
 gababq( 8,11) = gababq(4,10)
 gababq( 9,12) = 1.000
 gababq(10,12) = r4
 gababq(11,12) = z4
 gababq(12,12) = gababq(4,10)
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (12,gababq,12,d(10),0,d(11),ising,sp)
 
!     CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT
 
!       DELINT(04) = (0,0)
!       DELINT(05) = (0,1)
!       DELINT(06) = (0,2)
!       DELINT(07) = (1,0)
!       DELINT(08) = (1,1)
!       DELINT(09) = (1,2)
!       DELINT(10) = (2,0)
!       DELINT(11) = (2,1)
!       DELINT(12) = (2,2)
!       DELINT(13) = (3,0)
!       DELINT(14) = (3,1)
!       DELINT(15) = (3,2)
 
 i1 = 3
 DO  i = 1,4
   ip = i - 1
   DO  j = 1,3
     iq = j - 1
     i1 = i1 + 1
     delint(i1) = rzints(ip,iq,r,z,4)
   END DO
 END DO
 
!     LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3
 
 matidc = matid
 matflg = 7
 eltemp = tempe
 dgamr  = dgama*degrad
 sinth  = SIN(dgamr)
 costh  = COS(dgamr)
 sing   = sinth
 cosg   = costh
 CALL mat (idel)
 IF (setmat == 2.0) GO TO 910
 
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
 dela= 1.0/(1.0 - vro*vor - voz*vzo - vzr*vrz - vro*voz*vzr - vrz*vor*vzo)
 
!     COMPUTE ELASTIC CONSTANTS MATRIX FROM MATERIAL TO ELEMENT AXIS
 
 DO  i = 1,21
   teo( i) = 0.0
 END DO
 teo( 1) = er*(1.0 - voz*vzo)*dela
 teo( 2) = er*(vzr + vzo*vor)*dela
 teo( 3) = ez*(1.0 - vro*vor)*dela
 teo( 4) = er*(vor + vzr*voz)*dela
 teo( 5) = et*(vzo + vro*vzr)*dela
 teo( 6) = et*(1.0 - vrz*vzr)*dela
 teo(10) = grz
 teo(15) = gor
 teo(21) = gzo
 sing2   = sing**2
 cosg2   = cosg**2
 sing4   = sing**4
 cosg4   = cosg**4
 ee01    = teo(1)*cosg4 + teo(3)*sing4 + (teo(2) + 2.0*teo(10))*2.0  &
     *sing2*cosg2
 ee02    = teo(2)*(sing4 + cosg4) + (teo(1) + teo(3) - 4.0*teo(10))  &
     *sing2* cosg2
 ee03    = teo(1)*sing4 + teo(3)*cosg4 + (2.0*teo(2) + 4.0*teo(10))  &
     *sing2*cosg2
 ee04    = teo(4)*cosg2 + teo(5)*sing2
 ee05    = teo(4)*sing2 + teo(5)*cosg2
 ee06    = teo(6)
 ee07    = (teo(1)*cosg2 - teo(3)*sing2 + (teo(2) + 2.0*teo(10))  &
     *(sing2 - cosg2))*sing*cosg
 ee08    = (teo(1)*sing2 - teo(3)*cosg2 + (teo(2) + 2.0*teo(10))  &
     *(cosg2 - sing2))*sing*cosg
 ee09    = sing*cosg*(teo(4) - teo(5))
 
!     COMPUTE HARMONIC COEFFICIENT
 
 ajho = iecpt(1) - (iecpt(1)/1000)*1000 - 1
 
!     COMPUTE THERMAL LOAD
 
 a1 = ee01*alf(1) + ee02*alf(3) + ee04*alf(2)
 a2 = ee02*alf(1) + ee03*alf(3) + ee05*alf(2)
 a3 = ee04*alf(1) + ee05*alf(3) + ee06*alf(2)
 a4 = ee07*alf(1) + ee08*alf(3) + ee09*alf(2)
 
!     FORM HTN MATRIX
 
 htn( 1,1) = a3*delint(4)
 htn( 1,2) = a3*delint(7)
 htn( 1,3) = a3*delint(5)
 htn( 1,4) = a3*delint(8)
 htn( 2,1) = (a1+a3)*delint(7)
 htn( 2,2) = (a1+a3)*delint(10)
 htn( 2,3) = (a1+a3)*delint(8)
 htn( 2,4) = (a1+a3)*delint(11)
 htn( 3,1) = a3*delint( 5) + a4*delint( 7)
 htn( 3,2) = a3*delint( 8) + a4*delint(10)
 htn( 3,3) = a3*delint( 6) + a4*delint( 8)
 htn( 3,4) = a3*delint( 9) + a4*delint(11)
 htn( 4,1) = (a1+a3)*delint( 8) + a4*delint(10)
 htn( 4,2) = (a1+a3)*delint(11) + a4*delint(13)
 htn( 4,3) = (a1+a3)*delint( 9) + a4*delint(11)
 htn( 4,4) = (a1+a3)*delint(12) + a4*delint(14)
 htn( 5,1) = ajho*a3*delint(4)
 htn( 5,2) = ajho*a3*delint(7)
 htn( 5,3) = ajho*a3*delint(5)
 htn( 5,4) = ajho*a3*delint(8)
 htn( 6,1) = ajho*a3*delint(7)
 htn( 6,2) = ajho*a3*delint(10)
 htn( 6,3) = ajho*a3*delint(8)
 htn( 6,4) = ajho*a3*delint(11)
 htn( 7,1) = ajho*a3*delint(5)
 htn( 7,2) = ajho*a3*delint(8)
 htn( 7,3) = ajho*a3*delint(6)
 htn( 7,4) = ajho*a3*delint(9)
 htn( 8,1) = ajho*a3*delint(8)
 htn( 8,2) = ajho*a3*delint(11)
 htn( 8,3) = ajho*a3*delint(9)
 htn( 8,4) = ajho*a3*delint(12)
 htn( 9,1) = 0.0
 htn( 9,2) = 0.0
 htn( 9,3) = 0.0
 htn( 9,4) = 0.0
 htn(10,1) = a4*delint(7)
 htn(10,2) = a4*delint(10)
 htn(10,3) = a4*delint(8)
 htn(10,4) = a4*delint(11)
 htn(11,1) = a2*delint(7)
 htn(11,2) = a2*delint(10)
 htn(11,3) = a2*delint(8)
 htn(11,4) = a2*delint(11)
 htn(12,1) = a2*delint(10) + a4*delint( 8)
 htn(12,2) = a2*delint(13) + a4*delint(11)
 htn(12,3) = a2*delint(11) + a4*delint( 9)
 htn(12,4) = a2*delint(14) + a4*delint(12)
 
!     COMPUTE LITTLE H MATRIX (INVERSE OF PARTITION OF GABABQ)
 
 IF (ABS(r2-r1) < 1.0E-16) GO TO 930
 IF (ABS(r3-r4) < 1.0E-16) GO TO 930
 IF (ABS(z4-z1) < 1.0E-16) GO TO 930
 a      = 1.0/((r2-r1)*(r3-r4)*(z4-z1))
 r34a   = a*(r3-r4)
 r21a   = a*(r2-r1)
 h(1,1) = r34a*r2*z4
 h(1,2) =-r1*z4*r34a
 h(1,3) = r4*z1*r21a
 h(1,4) =-r3*z1*r21a
 h(2,1) =-z4*r34a
 h(2,2) = z4*r34a
 h(2,3) =-z1*r21a
 h(2,4) = z1*r21a
 h(3,1) =-r2*a*(r2-r4)
 h(3,2) = r1*r34a
 h(3,3) =-r4*r21a
 h(3,4) = r3*r21a
 h(4,1) = r34a
 h(4,2) =-r34a
 h(4,3) = r21a
 h(4,4) =-r21a
 
!     COMPUTE TI
 
 dgamr = tzero
 IF (ajho > 0.0) dgamr = 0.0
 DO  i = 1,4
   ti(i) = ti(i) - dgamr
 END DO
 
!     COMPUTE THE THEMAL LOAD IN FIELD COORDINATES
 
 CALL gmmats (h,  4, 4,1, ti(1),4,1,0, tl(1))
 CALL gmmats (htn,4,12,1, tl(1),4,1,0, d(1) )
 
!     TRANSFORM THE THERMAL LOAD TO GRID POINT DEGREES OF FREEDOM
!     ***  COORDINATE SYSTEMS NOT POSSIBLE  *******
!     ***  WITH RINGAX.  THE FOLLOWING CODE WILL IMPLEMENT IT. LRK ***
!     ***  THE FOLLOWING GMMATS HAS D(20) INSTEAD OF TL(1)        ****
 
 CALL gmmats (gababq,12,12,1, d(1),12,1,0,tl(1))
 
!     LOCATE THE TRANSFORMATION MATRICES FOR THE THREE GRID POINTS
!.    DO 750 I = 1,144
!.    AKI (I) = 0.0
!.750 CONTINUE
!.    DO 800 I = 1,4
!.    CALL GBTRAN (ICS(I),IECPT(4*I+20),D)   $ THIS IS WRONG ANYWAY
!.    K = 39*(I-1) + 1
!.    DO 800 J = 1,3
!.    KK = K + 12*(J-1)
!.    JJ = 3*(J-1) + 1
!.    AKI(KK  ) = D(JJ  )
!.    AKI(KK+1) = D(JJ+1)
!.    AKI(KK+2) = D(JJ+2)
!.800 CONTINUE
 
!     ADD THE ELEMENT THERMAL LOAD TO THE STRUCTURE THERMAL LOAD
 
!.    CALL GMMATS ( AKI(1), 12, 12, 1, D(20), 12, 1, 0, TL(1) )
 
 dgamr = pi
 IF (ajho == 0.0) dgamr = 2.0*pi
 
 DO  i = 1,12
   tl(i) = dgamr*tl(i)
 END DO
 
 k = 0
 DO  i = 1,4
   l = igp(i) - 1
   DO  j = 1,3
     k = k + 1
     l = l + 1
     pg(l) = pg(l) + tl(k)
   END DO
 END DO
 GO TO 950
 
 910 i = 37
 GO TO 925
 915 i = 218
 GO TO 935
 925 j =-30
 GO TO 945
 930 i = 31
 GO TO 925
 935 j = 30
 IF (idel1 == idel2) GO TO 950
 idel2 = idel1
 sp(2) = jax
 945 sp(1) = idel1
 CALL mesage (j,i,sp)
 950 RETURN
END SUBROUTINE tpztem
