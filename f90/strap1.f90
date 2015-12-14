SUBROUTINE strap1
     
!     THIS ROUTINE IS PHASE I OF STRESS DATA RECOVERY FOR THE
!     TRAPEZOIDAL CROSS SECTION RING
 
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
 
 DIMENSION       iecpt(24),ics(4),gambq(64),dzero(32),sp(24),  &
     alfb(4),teo(16),ee(16),delint(12),gamqs(96),jrz(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /condas/ consts(5)
 COMMON /sdr2x5/ ecpt(24),dum5(76),idel,igp(4),tz,sel(240),ts(4), ak(144)
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e(3),anu(3),rho,g(3),alf(3),tzero
 COMMON /sdr2x6/ d(144),gambl(144),r(5),z(5)
 COMMON /system/ ibuf,iout
 EQUIVALENCE     (consts(2),twopi),(consts(4),degra),  &
     (iecpt(1),ecpt(1)),(r(1),r1),(r(2),r2),(r(3),r3),  &
     (r(4),r4),(z(1),z1),(z(2),z2),(z(3),z3),(z(4),z4),  &
     (gambl(1),sp(1)),(gambl(1),teo(1)), (gambl(17),delint(1))
 
!     STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel   = iecpt(1)
 igp(1) = iecpt(2)
 igp(2) = iecpt(3)
 igp(3) = iecpt(4)
 igp(4) = iecpt(5)
 matid  = iecpt(7)
 ics(1) = iecpt(8)
 ics(2) = iecpt(12)
 ics(3) = iecpt(16)
 ics(4) = iecpt(20)
 r(1)   = ecpt(9)
 d(1)   = ecpt(10)
 z(1)   = ecpt(11)
 r(2)   = ecpt(13)
 d(2)   = ecpt(14)
 z(2)   = ecpt(15)
 r(3)   = ecpt(17)
 d(3)   = ecpt(18)
 z(3)   = ecpt(19)
 r(4)   = ecpt(21)
 d(4)   = ecpt(22)
 z(4)   = ecpt(23)
 tempe  = ecpt(24)
 dgama  = ecpt(6)
 
!     TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 DO  i = 1,4
   IF (r(i) < 0.0) CALL mesage (-30,37,idel)
   IF (d(i) /= 0.0) CALL mesage (-30,37,idel)
 END DO
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = AMIN1(z1,z2,z3,z4)
 z1   = z1 - zmin
 z2   = z2 - zmin
 z3   = z3 - zmin
 z4   = z4 - zmin
 rmin = AMIN1(r1,r2,r3,r4)
 rmax = AMAX1(r1,r2,r3,r4)
 IF (rmin == 0.) GO TO 206
 IF (rmax/rmin <= 10.) GO TO 206
 
!     RATIO OF RADII IS TOO LARGE FOR GAUSS QUADRATURE FOR IP=-1
 
 WRITE  (iout,205) ufm,idel
 205 FORMAT (a23,', TRAPRG ELEMENT',i9,' HAS A MAXIMUM TO MINIMUM ',  &
     'RADIUS RATIO EXCEEDING 10.', /5X,  &
     'ACCURACY OF NUMERICAL INTEGRATION WOULD BE IN DOUBT.')
 CALL mesage (-30,37,idel)
 206 CONTINUE
 icore = 0
 j = 1
 DO  i = 1,4
   IF (r(i) /= 0.) CYCLE
   icore  = icore + 1
   jrz(j) = i
   j = 2
 END DO
 IF (icore /= 0 .AND. icore /= 2) CALL mesage (-30,37,idel)
 
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
 
!     MODIFY THE TRANSFORMATION MATRIX IF ELEMENT IS A CORE ELEMENT
 
 IF (icore == 0) GO TO 305
 jj1 = 2*jrz(1) - 1
 jj2 = 2*jrz(2) - 1
 
 DO  i = 1,8
   j = 8*(i-1)
   gambq(i    ) = 0.0
   gambq(i+ 16) = 0.0
   gambq(j+jj1) = 0.
   gambq(j+jj2) = 0.
 END DO
 305 CONTINUE
 
!     CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT WHERE THE ORDER IS
!     INDICATED BY THE FOLLOWING TABLE
 
!       DELINT( 1) - (-1,0)
!       DELINT( 2) - (-1,1)
!       DELINT( 3) - (-1,2)
!       DELINT( 4) - ( 0,0)
!       DELINT( 5) - ( 0,1)
!       DELINT( 6) - ( 0,2)
!       DELINT( 7) - ( 1,0)
!       DELINT( 8) - ( 1,1)
!       DELINT( 9) - ( 1,2)
!       DELINT(10) - ( 2,0)
!       DELINT(11) - ( 2,1)
!       DELINT(12) - ( 3,0)
 
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
     IF (icore == 0) GO TO 345
     IF (i1    > 3) GO TO 345
     delint(i1) = 0.0
     CYCLE
     345 CONTINUE
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
 
 CALL gmmats (teo,4,4,1, ee ,4,4,0, d )
 CALL gmmats (d  ,4,4,0, teo,4,4,0, ee)
 
!     FORM THE ELEMENT STIFFNESS MATRIX IN FIELD COORDINATES
 
 ee48   = ee(4) + ee(8)
 d ( 1) = ee(1) + 2.0 * ee(2) + ee(6)
 ak( 1) = ee(6) * delint(1)
 ak( 2) = (ee(2) + ee(6)) * delint(4)
 ak( 3) = ee(6) * delint(2) + ee(8) * delint(4)
 ak( 4) = (ee(2) + ee(6)) * delint(5) + ee(8) * delint(7)
 ak( 5) = 0.0
 ak( 6) = ee(8) * delint(4)
 ak( 7) = ee(7) * delint(4)
 ak( 8) = ee(7) * delint(7) + ee(8) * delint(5)
 ak( 9) = ak(2)
 ak(10) = d(1) * delint(7)
 ak(11) = (ee(2) + ee(6)) * delint(5) + ee48 * delint(7)
 ak(12) = d(1) * delint(8) + ee48 * delint(10)
 ak(13) = 0.0
 ak(14) = ee48 * delint(7)
 ak(15) = (ee(3) + ee(7)) * delint(7)
 ak(16) = (ee(3) + ee(7)) * delint(10) + ee48 * delint(8)
 ak(17) = ak( 3)
 ak(18) = ak(11)
 ak(19) = ee(6) * delint(3) + ee(16)* delint(7) + (ee(8) + ee(14)) * delint(5)
 ak(20) = (ee(2) + ee(6)) * delint(6) + ee(16) * delint(10)  &
     + (ee(8) + ee(13) + ee(14)) * delint(8)
 ak(21) = 0.0
 ak(22) = ee(16) * delint(7) + ee(8) * delint(5)
 ak(23) = ee(7) * delint(5) + ee(15) * delint(7)
 ak(24) = (ee(7) + ee(16)) * delint(8)  &
     + ee(8) *delint(6) + ee(15) * delint(10)
 ak(25) = ak(4)
 ak(26) = ak(12)
 ak(27) = ak(20)
 ak(28) = d(1) * delint(9) + ee(16) * delint(12)  &
     + (ee48 + ee(13) + ee(14)) * delint(11)
 ak(29) = 0.0
 ak(30) = ee(16) * delint(10) + ee48 * delint(8)
 ak(31) = (ee(3) + ee(7)) * delint(8) + ee(15) * delint(10)
 ak(32) = (ee(3) + ee(7) + ee(16)) * delint(11)  &
     + ee(15) * delint(12) + ee48 * delint(9)
 ak(33) = 0.0
 ak(34) = 0.0
 ak(35) = 0.0
 ak(36) = 0.0
 ak(37) = 0.0
 ak(38) = 0.0
 ak(39) = 0.0
 ak(40) = 0.0
 ak(41) = ak( 6)
 ak(42) = ak(14)
 ak(43) = ak(22)
 ak(44) = ak(30)
 ak(45) = 0.0
 ak(46) = ee(16)*delint(7)
 ak(47) = ee(15)*delint(7)
 ak(48) = ee(16)*delint(8) + ee(15) * delint(10)
 ak(49) = ak( 7)
 ak(50) = ak(15)
 ak(51) = ak(23)
 ak(52) = ak(31)
 ak(53) = 0.0
 ak(54) = ak(47)
 ak(55) = ee(11)*delint( 7)
 ak(56) = ee(11)*delint(10) + ee(12) * delint(8)
 ak(57) = ak( 8)
 ak(58) = ak(16)
 ak(59) = ak(24)
 ak(60) = ak(32)
 ak(61) = 0.0
 ak(62) = ak(48)
 ak(63) = ak(56)
 ak(64) = ee(11) * delint(12) + ee(16) * delint(9)  &
     + (ee(12) + ee(13)) * delint(11)
 
 DO  i = 1,64
   ak(i) = twopi*ak(i)
 END DO
 
!     TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM FIELD COORDINATES
!     TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmats (gambq,8,8,1, ak   ,8,8,0, d )
 CALL gmmats (d    ,8,8,0, gambq,8,8,0, ak)
 
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
 
!     TRANSFORM THE STIFFNESS MATRIX FROM TWO TO THREE DEGREES OF
!     FREEDOM PER POINT
 
 CALL gmmats (gamqs(1),8,12,1, ak(1)   ,8, 8,0, d(1) )
 CALL gmmats (d(1)    ,12,8,0, gamqs(1),8,12,0, ak(1))
 
!     LOCATE THE TRANSFORMATION MATRICES FOR THE FOUR  GRID POINTS
 
 DO  i = 1,144
   gambl(i) = 0.0
 END DO
 DO  i = 1,4
   CALL transs (ics(i),d(1))
   k = 39*(i-1) + 1
   DO  j = 1,3
     kk = k + 12*(j-1)
     jj = 3 *(j-1) + 1
     gambl(kk  ) = d(jj  )
     gambl(kk+1) = d(jj+1)
     gambl(kk+2) = d(jj+2)
   END DO
 END DO
 
!     TRANSFORM THE STIFFNESS MATRIX FROM BASIC TO LOCAL COORDINATES
 
 CALL gmmats (gambl(1),12,12,1, ak(1)   ,12,12,0, d(1) )
 CALL gmmats (d(1)    ,12,12,0, gambl(1),12,12,0, ak(1))
 
!     COMPUTE THE FIFTH GRID POINT TO BE THE AVERAGE OF THE FOUR
!     CORNER POINTS
 
 r(5) = (r1 + r2 + r3 + r4)/4.0
 z(5) = (z1 + z2 + z3 + z4)/4.0
 
!     INITIALIZE THE CONSTANT PORTION OF THE D SUB 0 MATRIX
 
 DO  i = 1,32
   dzero(i) = 0.0
 END DO
 dzero( 2) = 1.0
 dzero(10) = 1.0
 dzero(23) = 1.0
 dzero(27) = 1.0
 dzero(30) = 1.0
 
!     START THE LOOP TO COMPUTE THE STRESS MATRIX FOR EACH GRID POINT
 
 DO  j = 1,5
   
!     COMPUTE THE VARIABLE PORTION OF THE D SUB 0 MATRIX
   
   dzero( 4) = z(j)
   IF (icore /= 0) GO TO 875
   dzero( 9) = 1.00/r(j)
   dzero(11) = z(j)/r(j)
   875 CONTINUE
   dzero(12) = z(j)
   dzero(24) = r(j)
   dzero(28) = r(j)
   dzero(32) = z(j)
   
!     COMPUTE THE STRESS MATRIX IN FIELD COORDINATES
   
   CALL gmmats (ee(1),4,4,0, dzero(1),4,8,0, d(1))
   
!     TRANSFORM THE STRESS MATRIX TO GRID POINT DEGREES OF FREEDOM
   
   CALL gmmats (d(1),4,8,0, gambq(1),8,8,0, d(37))
   
!     TRANSFORM THE STRESS MATRIX FROM TWO TO THREE DEGREES OF FREEDOM
!     PER POINT
   
   CALL gmmats (d(37),4,8,0, gamqs(1),8,12,0, d(73))
   
!     TRANSFORM THE STRESS MATRIX FROM BASIC TO LOCAL COORDINATES
   
   k = 48*(j-1) + 1
   CALL gmmats (d(73),4,12,0, gambl(1),12,12,0, sel(k))
   
 END DO
 
!     COMPUTE THE THERMAL STRAIN VECTOR
 
 DO  i = 1,3
   alfb(i) = alf(i)
 END DO
 alfb(4) = 0.0
 
!     COMPUTE THE THERMAL STRESS VECTOR
 
 CALL gmmats (ee(1),4,4,0, alfb(1),4,1,0, ts(1))
 RETURN
END SUBROUTINE strap1
