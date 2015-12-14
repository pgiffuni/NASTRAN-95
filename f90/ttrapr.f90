SUBROUTINE ttrapr (ti,pg)
     
!     THIS ROUTINE COMPUTES THE THERMAL LOAD FOR THE TRAPEZOIDAL
!     CROSS SECTION RING
 
!     ECPT FOR THE TRAPEZOIDAL RING
!                                                          TYPE
!     ECPT( 1) ELEMENT IDENTIFICATION                        I
!     ECPT( 2) SCALAR INDEX NO. FOR GRID POINT A             I
!     ECPT( 3) SCALAR INDEX NO. FOR GRID POINT B             I
!     ECPT( 4) SCALAR INDEX NO. FOR GRID POINT C             I
!     ECPT( 5) SCALAR INDEX NO. FOR GRID POINT D             I
!     ECPT( 6) MATERIAL ORIENTATION ANGLE(DEGREES)           R
!     ECPT( 7) MATERIAL IDENTIFICATION                       I
!     ECPT( 8) COOR. SYS. ID. FOR GRID POINT A               I
!     ECPT( 9) X-COOR. OF GRID POINT A (IN BASIC COOR.)      R
!     ECPT(10) Y-COOR. OF GRID POINT A (IN BASIC COOR.)      R
!     ECPT(11) Z-COOR. OF GRID POINT A (IN BASIC COOR.)      R
!     ECPT(12) COOR. SYS. ID. FOR GRID POINT B               I
!     ECPT(13) X-COOR. OF GRID POINT B (IN BASIC COOR.)      R
!     ECPT(14) Y-COOR. OF GRID POINT B (IN BASIC COOR.)      R
!     ECPT(15) Z-COOR. OF GRID POINT B (IN BASIC COOR.)      R
!     ECPT(16) COOR. SYS. ID. FOR GRID POINT C               I
!     ECPT(17) X-COOR. OF GRID POINT C (IN BASIC COOR.)      R
!     ECPT(18) Y-COOR. OF GRID POINT C (IN BASIC COOR.)      R
!     ECPT(19) Z-COOR. OF GRID POINT C (IN BASIC COOR.)      R
!     ECPT(20) COOR. SYS. ID. FOR GRID POINT D               I
!     ECPT(21) X-COOR. OF GRID POINT D (IN BASIC COOR.)      R
!     ECPT(22) Y-COOR. OF GRID POINT D (IN BASIC COOR.)      R
!     ECPT(23) Z-COOR. OF GRID POINT D (IN BASIC COOR.)      R
!     ECPT(24) EL. TEMPERATURE FOR MATERIAL PROPERTIES       R
 
 
 
 REAL, INTENT(OUT)                        :: ti(4)
 REAL, INTENT(OUT)                        :: pg(1)
 DIMENSION  iecpt(24),d(22),gambq(64),r(4),z(4),  &
     teo(16),ee(16),delint(12),gamqs(96),q(32),  &
     gambl(144),alfb(4),igp(4),ics(4),sp(24),hprim(16), tl(12),ts(4),jrz(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /condas/ consts(5)
 COMMON /trimex/ ecpt(24)
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e(3),anu(3),rho,g(3),alf(3),tzero
 COMMON /system/ ibuf,iout
 EQUIVALENCE     (consts(2),twopi),(consts(4),degra), (iecpt(1),ecpt(1)),  &
     (r(1),r1),(r(2),r2),(r(3),r3),(r(4),r4),  &
     (z(1),z1),(z(2),z2),(z(3),z3),(z(4),z4),  &
     (gambl(1),ee(1)),(gambl(17),teo(1)),  &
     (gambl(33),alfb(1)),(gambl(37),ts(1)),  &
     (gambl(41),delint(1)),(gambl(1),gambq(1)),  &
     (gambl(65),q(1)),(gambl(97),hprim(1)), (gambl(113),sp(1)),(gambl(1),gamqs(1))
 
!     STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel   = iecpt( 1)
 igp(1) = iecpt( 2)
 igp(2) = iecpt( 3)
 igp(3) = iecpt( 4)
 igp(4) = iecpt( 5)
 matid  = iecpt( 7)
 ics(1) = iecpt( 8)
 ics(2) = iecpt(12)
 ics(3) = iecpt(16)
 ics(4) = iecpt(20)
 r(1)   = ecpt ( 9)
 d(1)   = ecpt (10)
 z(1)   = ecpt (11)
 r(2)   = ecpt (13)
 d(2)   = ecpt (14)
 z(2)   = ecpt (15)
 r(3)   = ecpt (17)
 d(3)   = ecpt (18)
 z(3)   = ecpt (19)
 r(4)   = ecpt (21)
 d(4)   = ecpt (22)
 z(4)   = ecpt (23)
 tempe  = ecpt (24)
 dgama  = ecpt ( 6)
 
!     TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 DO  i = 1,4
   IF (r(i) < 0.0) CALL mesage (-30,37,idel)
   IF (d(i) /= 0.0) CALL mesage (-30,37,idel)
 END DO
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = AMIN1(z1,z2,z3,z4)
 z1 = z1 - zmin
 z2 = z2 - zmin
 z3 = z3 - zmin
 z4 = z4 - zmin
 rmin = AMIN1(r1,r2,r3,r4)
 rmax = AMAX1(r1,r2,r3,r4)
 IF (rmin == 0.) GO TO 206
 IF (rmax/rmin <= 10.) GO TO 206
 
!     RATIO OF RADII IS TOO LARGE FOR GAUSS QUADRATURE FOR IP=-1
 
 WRITE  (iout,205) ufm,idel
 205 FORMAT (a23,', TRAPRG ELEMENT',i9,' HAS A MAXIMUM TO MINIMUM ',  &
     'RADIUS RATIO EXCEEDING 10.'/5X,'ACCURACY OF NUMERICAL ',  &
     'INTEGRATION WOULD BE IN DOUBT.')
 CALL mesage (-61,0,0)
 206 CONTINUE
 icore = 0
 j = 1
 DO  i = 1,4
   IF (r(i) /= 0.) CYCLE
   icore  = icore + 1
   jrz(j) = i
   j = 2
 END DO
 IF (icore /= 0 .AND. icore /= 2) CALL mesage (-61,0,0)
 
!     CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT WHERE THE ORDER IS
!     INDICATED BY THE FOLLOWING TABLE
 
!       DELINT( 1) - ( 0,0)
!       DELINT( 2) - ( 0,1)
!       DELINT( 3) - ( 0,2)
!       DELINT( 4) - ( 1,0)
!       DELINT( 5) - ( 1,1)
!       DELINT( 6) - ( 1,2)
!       DELINT( 7) - ( 2,0)
!       DELINT( 8) - ( 2,1)
!       DELINT( 9) - ( 2,2)
!       DELINT(10) - ( 3,0)
!       DELINT(11) - ( 3,1)
 
 i1 = 0
 DO  i = 1,4
   ip = i - 1
   DO  j = 1,3
     iq = j - 1
     i1 = i1 + 1
     IF (i1 == 12) CYCLE
     delint(i1) = rzints(ip,iq,r,z,4)
   END DO
 END DO
 
!     LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3 TABLE
 
 matidc = matid
 matflg = 7
 eltemp = tempe
 CALL  mat (idel)
 
!     SET MATERIAL PROPERTIES IN LOCAL VARIABLES
 
 er  = e(1)
 et  = e(2)
 ez  = e(3)
 vrt = anu(1)
 vtz = anu(2)
 vzr = anu(3)
 grz = g(3)
 tz  = tzero
 vtr = vrt*et/er
 vzt = vtz*ez/et
 vrz = vzr*er/ez
 del = 1.0 - vrt*vtr - vtz*vzt - vzr*vrz - vrt*vtz*vzr - vrz*vtr*vzt
 
!     GENERATE ELASTIC CONSTANTS MATRIX (4X4)
 
 ee( 1) = er*(1.0 - vtz*vzt)/del
 ee( 2) = er*(vtr + vzr*vtz)/del
 ee( 3) = er*(vzr + vtr*vzt)/del
 ee( 4) = 0.0
 ee( 5) = ee(2)
 ee( 6) = et*(1.0 - vrz*vzr)/del
 ee( 7) = et*(vzt + vrt*vzr)/del
 ee( 8) = 0.0
 ee( 9) = ee(3)
 ee(10) = ee(7)
 ee(11) = ez*(1.0 - vrt*vtr)/del
 ee(12) = 0.0
 ee(13) = 0.0
 ee(14) = 0.0
 ee(15) = 0.0
 ee(16) = grz
 
!     FORM TRANSFORMATION MATRIX (4X4) FROM MATERIAL AXIS TO ELEMENT
!     GEOMETRIC AXIS
 
 dgamr   = dgama*degra
 cosg    = COS(dgamr)
 sing    = SIN(dgamr)
 teo( 1) = cosg**2
 teo( 2) = 0.0
 teo( 3) = sing**2
 teo( 4) = sing*cosg
 teo( 5) = 0.0
 teo( 6) = 1.0
 teo( 7) = 0.0
 teo( 8) = 0.0
 teo( 9) = teo(3)
 teo(10) = 0.0
 teo(11) = teo(1)
 teo(12) =-teo(4)
 teo(13) =-2.0*teo(4)
 teo(14) = 0.0
 teo(15) =-teo(13)
 teo(16) = teo(1) - teo(3)
 
!     TRANSFORM THE ELASTIC CONSTANTS MATRIX FROM MATERIAL
!     TO ELEMENT GEOMETRIC AXIS
 
 CALL gmmats (teo,4,4,1, ee, 4,4,0, d )
 CALL gmmats (d  ,4,4,0, teo,4,4,0, ee)
 
!     COMPUTE THE THERMAL STRAIN VECTOR
 
 DO  i = 1,3
   alfb(i) = alf(i)
 END DO
 alfb(4) = 0.0
 
 CALL gmmats (ee(1),4,4,0, alfb(1),4,1,0, ts(1))
 
!     FORM THE Q MATRIX (8X4)
 
 d( 1) = ts(1) + ts(2)
 q( 1) = ts(2)*delint(1)
 q( 2) = ts(2)*delint(4)
 q( 3) = ts(2)*delint(2)
 q( 4) = ts(2)*delint(5)
 q( 5) =  d(1)*delint(4)
 q( 6) =  d(1)*delint(7)
 q( 7) =  d(1)*delint(5)
 q( 8) =  d(1)*delint(8)
 q( 9) = ts(2)*delint(2)
 q(10) = ts(2)*delint(5)
 q(11) = ts(2)*delint(3)
 q(12) = ts(2)*delint(6)
 q(13) =  d(1)*delint(5)
 q(14) =  d(1)*delint(8)
 q(15) =  d(1)*delint(6)
 q(16) =  d(1)*delint(9)
 DO  i = 17,24
   q( i) = 0.0
 END DO
 q(25) = ts(3)*delint( 4)
 q(26) = ts(3)*delint( 7)
 q(27) = ts(3)*delint( 5)
 q(28) = ts(3)*delint( 8)
 q(29) = ts(3)*delint( 7)
 q(30) = ts(3)*delint(10)
 q(31) = ts(3)*delint( 8)
 q(32) = ts(3)*delint(11)
 
!     FORM THE TRANSFORMATION MATRIX (8X8) FROM FIELD COORDINATES TO
!     GRID POINT DEGREES OF FREEDOM
 
 DO  i = 1,64
   gambq(i) = 0.0
 END DO
 gambq( 1) = 1.0
 gambq( 2) = r1
 gambq( 3) = z1
 gambq( 4) = r1*z1
 gambq(13) = 1.0
 gambq(14) = r1
 gambq(15) = z1
 gambq(16) = gambq(4)
 gambq(17) = 1.0
 gambq(18) = r2
 gambq(19) = z2
 gambq(20) = r2*z2
 gambq(29) = 1.0
 gambq(30) = r2
 gambq(31) = z2
 gambq(32) = gambq(20)
 gambq(33) = 1.0
 gambq(34) = r3
 gambq(35) = z3
 gambq(36) = r3*z3
 gambq(45) = 1.0
 gambq(46) = r3
 gambq(47) = z3
 gambq(48) = gambq(36)
 gambq(49) = 1.0
 gambq(50) = r4
 gambq(51) = z4
 gambq(52) = r4*z4
 gambq(61) = 1.0
 gambq(62) = r4
 gambq(63) = z4
 gambq(64) = gambq(52)
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (8,gambq(1),8,d(10),0,d(11),ising,sp)
 
 IF (ising == 2) CALL mesage (-30,26,idel)
 
!     FORM THE HPRIM MATRIX (4X4)
 
 k  = 0
 DO  i = 1,4
   kk =  8*(i-1) - 1
   DO  j = 1,4
     k  = k  + 1
     kk = kk + 2
     hprim(k) = gambq(kk)
   END DO
 END DO
 
!     MODIFY THE TRANSFORMATION MATRIX IF ELEMENT IS A CORE ELEMENT
 
 IF (icore == 0) GO TO 665
 jj1 = 2*jrz(1) - 1
 jj2 = 2*jrz(2) - 1
 
 DO  i = 1,8
   j = 8*(i-1)
   gambq(i    ) = 0.0
   gambq(i+ 16) = 0.0
   gambq(j+jj1) = 0.
   gambq(j+jj2) = 0.
 END DO
 665 CONTINUE
 
!     FORM THE TEMPERATURE VECTOR
 
 DO  i = 1,4
   ti(i) = ti(i) - tzero
 END DO
 
!     COMPUTE THE THERMAL LOAD IN FIELD COORDINATES
 
 CALL gmmats (hprim(1),4,4,0, ti(1),4,1,0, tl(1))
 CALL gmmats (q(1),    8,4,0, tl(1),4,1,0,  d(1))
 
!     TRANSFORM THE THERMAL LOAD TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmats (gambq(1),8,8,1, d(1),8,1,0, tl(1))
 
!     GENERATE THE TRANSFORMATION MATRIX FROM TWO TO THREE DEGREES OF
!     FREEDOM PER POINT
 
 DO  i = 1,96
   gamqs( i) = 0.0
 END DO
 gamqs( 1) = 1.0
 gamqs(15) = 1.0
 gamqs(28) = 1.0
 gamqs(42) = 1.0
 gamqs(55) = 1.0
 gamqs(69) = 1.0
 gamqs(82) = 1.0
 gamqs(96) = 1.0
 
!     TRANSFORM THE   THERMAL LOAD   FROM TWO TO THREE DEGREES OF
!     FREEDOM PER POINT
 
 CALL gmmats (gamqs(1),8,12,1, tl(1),8,1,0, d(10))
 
!     LOCATE THE TRANSFORMATION MATRICES FOR THE FOUR  GRID POINTS
 
 DO  i = 1,144
   gambl(i) = 0.0
 END DO
 DO  i = 1,4
   CALL gbtran (ics(i),ecpt(4*i+4),d(1))
   k  = 39*(i-1) + 1
   DO  j = 1,3
     kk = k + 12*(j-1)
     jj = 3*(j-1) + 1
     gambl(kk  ) = d(jj  )
     gambl(kk+1) = d(jj+1)
     gambl(kk+2) = d(jj+2)
   END DO
 END DO
 
!     TRANSFORM THE   THERMAL LOAD   FROM BASIC TO LOCAL COORDINATES
 
 CALL gmmats (gambl(1),12,12,1, d(10),12,1,0, tl(1))
 DO  i = 1,12
   tl(i) = twopi*tl(i)
 END DO
 
!     ADD THE ELEMENT THERMAL LOAD TO THE STRUCTURE THERMAL LOAD
 
 k = 0
 DO  i = 1,4
   l = igp(i) - 1
   DO  j = 1,3
     k = k + 1
     l = l + 1
     pg(l) = pg(l) + tl(k)
   END DO
 END DO
 
 RETURN
END SUBROUTINE ttrapr
