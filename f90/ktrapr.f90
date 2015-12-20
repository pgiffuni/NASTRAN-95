SUBROUTINE ktrapr
     
!     THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR A AXI-SYMMETRIC
!     RING WITH A TRAPEZOIDAL CROSS SECTION
 
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
 
 DOUBLE PRECISION :: constd,   degrad,   d ,       gambq,    r,  &
     z,        ee48,     teo,      ee,       delint,  &
     ak,       aki,      r1,       r2,       r3,  &
     r4,       z1,       z2,       z3,       z4,  &
     zmin,     dgama,    er,       et,       ez,  &
     vrt,      vtr,      vtz,      vzt,      vzr,  &
     vrz,      grz,      del,      cosg,     sing,  &
     dgamr,    akt,      twopi,    dampc,    rmin, rmax,     rzintd
 DIMENSION        jrz(2),   iecpt(24),aki(36),  akt(9)
 COMMON /system/  ibuf,     iout
 COMMON /condad/  constd(5)
 COMMON /msgx  /  nmsg,     mmsg,     msg(4,1)
 COMMON /sma1io/  dum1(10), ifkgg,    igkgg,    if4gg,    dum2(21)
 COMMON /sma1cl/  iopt4,    k4ggsw,   npvt,     dum4(7),  link(10),  &
     idetck,   dodet,    nogo
 COMMON /sma1et/  ecpt(24), dum5(76)
 COMMON /sma1dp/  d(64),    gambq(64), r(4),    z(4),     teo(16),  &
     ee(16),   delint(12),ak(64),  dgama,    zmin,  &
     er,       et,       ez,       vrt,      vtr,  &
     vtz,      vzt,      vzr,      vrz,      grz,  &
     del,      cosg,     sing,     dgamr,    igp(4), ics(4),   sp(24),   tempe
 COMMON /matin /  matidc,   matflg,   eltemp,   stress,   sinth, costh
 COMMON /matout/  e(3),     anu(3),   rho,      g(3),     alf(3),  &
     tzero,    g_sub_e
 EQUIVALENCE      (constd(2),twopi),  (constd(4),degrad),  &
     (iecpt(1) ,ecpt(1)), (r(1),r1),(r(2),r2),(r(3),r3), (r(4),r4),  &
     (z(1),z1),(z(2),z2),(z(3),z3), (z(4),z4),  &
     (aki(1),gambq(1))  ,(akt(1),gambq(37))
 DATA    irg   /  4HTRAP    /
 
 
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
 r(1)   = ecpt( 9)
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
 dgama  = ecpt( 6)
 
!     CHECK INTERNAL GRID POINTS FOR PIVOT POINT
 
 ipp = 0
 DO  i = 1,4
   IF (npvt == igp(i)) ipp = i
 END DO
 IF (ipp == 0) CALL mesage (-30,34,idel)
 
!     TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 DO  i = 1,4
   IF (r(i) < 0.0D0) GO TO 910
   IF (d(i) /= 0.0D0) GO TO 910
 END DO
 
!     COMPUTE THE ELEMENT COORDINATES
 
 zmin = DMIN1(z1,z2,z3,z4)
 z1 = z1 - zmin
 z2 = z2 - zmin
 z3 = z3 - zmin
 z4 = z4 - zmin
 
!     FATAL IF RATIO OF RADII IS TO LARGE FOR GUASS QUADRATURE FOR
!     IP =-1
 
 rmin = DMIN1(r1,r2,r3,r4)
 rmax = DMAX1(r1,r2,r3,r4)
 IF (rmin == 0.d0) GO TO 206
 IF (rmax/rmin > 10.d0) GO TO 930
 
 206 IF (r1 >= r2 .OR. r4 >= r3 .OR. z4 <= z1) GO TO 910
 IF (DABS(z1-z2) > 1.0D-3) GO TO 910
 IF (DABS(z3-z4) > 1.0D-3) GO TO 910
 d(5) = (r1+r4)/2.0D0
 d(6) = (r2+r3)/2.0D0
 IF (d(5) == 0.0D0) GO TO 210
 IF (DABS((r1-r4)/d(5)) > 0.5D-2) GO TO 210
 r1 = d(5)
 r4 = d(5)
 210 CONTINUE
 IF (d(6) == 0.0D0) GO TO 220
 IF (DABS((r2-r3)/d(6)) > 0.5D-2) GO TO 220
 r(2) = d(6)
 r(3) = d(6)
 220 CONTINUE
 
 icore = 0
 j = 1
 DO  i = 1,4
   IF (r(i) /= 0.d0) CYCLE
   icore  = icore + 1
   jrz(j) = i
   j = 2
 END DO
 IF (icore /= 0 .AND. icore /= 2) GO TO 910
 
!     FORM THE TRANSFORMATION MATRIX (8X8) FROM FIELD COORDINATES TO
!     GRID POINT DEGREES OF FREEDOM
 
 DO  i = 1,64
   gambq(i) = 0.0D0
 END DO
 gambq( 1) = 1.0D0
 gambq( 2) = r1
 gambq( 3) = z1
 gambq( 4) = r1*z1
 gambq(13) = 1.0D0
 gambq(14) = r1
 gambq(15) = z1
 gambq(16) = gambq(4)
 gambq(17) = 1.0D0
 gambq(18) = r2
 gambq(19) = z2
 gambq(20) = r2*z2
 gambq(29) = 1.0D0
 gambq(30) = r2
 gambq(31) = z2
 gambq(32) = gambq(20)
 gambq(33) = 1.0D0
 gambq(34) = r3
 gambq(35) = z3
 gambq(36) = r3*z3
 gambq(45) = 1.0D0
 gambq(46) = r3
 gambq(47) = z3
 gambq(48) = gambq(36)
 gambq(49) = 1.0D0
 gambq(50) = r4
 gambq(51) = z4
 gambq(52) = r4*z4
 gambq(61) = 1.0D0
 gambq(62) = r4
 gambq(63) = z4
 gambq(64) = gambq(52)
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL inverd (8,gambq(1),8,d(10),0,d(11),ising,sp)
 IF (ising == 2) GO TO 920
 
!     MODIFY THE TRANSFORMATION MATRIX IF ELEMENT IS A CORE ELEMENT
 
 IF (icore == 0) GO TO 305
 jj1 = 2*jrz(1) - 1
 jj2 = 2*jrz(2) - 1
 
 DO  i = 1,8
   j = 8*(i-1)
   gambq(i   ) = 0.0D0
   gambq(i+16) = 0.0D0
   gambq(j+jj1)= 0.d0
   gambq(j+jj2)= 0.d0
 END DO
 305 CONTINUE
 
!     CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT WHERE THE ORDER IS
!     INDICATED BY THE FOLLOWING TABLE
 
!        DELINT( 1) - (-1,0)
!        DELINT( 2) - (-1,1)
!        DELINT( 3) - (-1,2)
!        DELINT( 4) - ( 0,0)
!        DELINT( 5) - ( 0,1)
!        DELINT( 6) - ( 0,2)
!        DELINT( 7) - ( 1,0)
!        DELINT( 8) - ( 1,1)
!        DELINT( 9) - ( 1,2)
!        DELINT(10) - ( 2,0)
!        DELINT(11) - ( 2,1)
!        DELINT(12) - ( 3,0)
 
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
     IF (i1 > 3) GO TO 345
     delint(i1) = 0.0D0
     CYCLE
     345 CONTINUE
     delint(i1) = rzintd(ip,iq,r,z,4)
   END DO
 END DO
 
!     LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3 TABLE
 
 matidc = matid
 matflg = 7
 eltemp = tempe
 CALL mat (idel)
 
!     SET MATERIAL PROPERTIES IN DOUBLE PRECISION VARIABLES
 
 er  = e(1)
 et  = e(2)
 ez  = e(3)
 vrt = anu(1)
 vtz = anu(2)
 vzr = anu(3)
 grz = g(3)
 vtr = vrt*et/er
 vzt = vtz*ez/et
 vrz = vzr*er/ez
 del = 1.0D0 - vrt*vtr - vtz*vzt - vzr*vrz - vrt*vtz*vzr - vrz*vtr*vzt
 
!     GENERATE ELASTIC CONSTANTS MATRIX (4X4)
 
 ee(1) = er*(1.0D0-vtz*vzt)/del
 ee(2) = er*(vtr + vzr*vtz)/del
 ee(3) = er*(vzr + vtr*vzt)/del
 ee(4) = 0.0D0
 ee(5) = ee(2)
 ee(6) = et*(1.0D0-vrz*vzr)/del
 ee(7) = et*(vzt+vrt*vzr)/ del
 ee(8) = 0.0D0
 ee(9) = ee(3)
 ee(10)= ee(7)
 ee(11)= ez*(1.0D0-vrt*vtr)/del
 ee(12)= 0.0D0
 ee(13)= 0.0D0
 ee(14)= 0.0D0
 ee(15)= 0.0D0
 ee(16)= grz
 
!     FORM TRANSFORMATION MATRIX (4X4) FROM MATERIAL AXIS TO ELEMENT
!     GEOMETRIC AXIS
 
 dgamr   = dgama*degrad
 cosg    = DCOS(dgamr)
 sing    = DSIN(dgamr)
 teo( 1) = cosg**2
 teo( 2) = 0.0D0
 teo( 3) = sing**2
 teo( 4) = sing*cosg
 teo( 5) = 0.0D0
 teo( 6) = 1.0D0
 teo( 7) = 0.0D0
 teo( 8) = 0.0D0
 teo( 9) = teo(3)
 teo(10) = 0.0D0
 teo(11) = teo(1)
 teo(12) =-teo(4)
 teo(13) =-2.0D0*teo(4)
 teo(14) = 0.0D0
 teo(15) =-teo(13)
 teo(16) = teo(1) - teo(3)
 
!     TRANSFORM THE ELASTIC CONSTANTS MATRIX FROM MATERIAL
!     TO ELEMENT GEOMETRIC AXIS
 
 CALL gmmatd (teo,4,4,1, ee, 4,4,0, d )
 CALL gmmatd (d  ,4,4,0, teo,4,4,0, ee)
 
!     FORM THE ELEMENT STIFFNESS MATRIX IN FIELD COORDINATES
 
 ee48   = ee(4) + ee(8)
 d ( 1) = ee(1) + 2.0D0*ee(2) + ee(6)
 ak( 1) = ee(6)*delint(1)
 ak( 2) = (ee(2) + ee(6))*delint(4)
 ak( 3) = ee(6)*delint(2) + ee(8)*delint(4)
 ak( 4) = (ee(2) + ee(6))*delint(5) + ee(8)*delint(7)
 ak( 5) = 0.0D0
 ak( 6) = ee(8)*delint(4)
 ak( 7) = ee(7)*delint(4)
 ak( 8) = ee(7)*delint(7) + ee(8)*delint(5)
 ak( 9) = ak(2)
 ak(10) = d(1)*delint(7)
 ak(11) = (ee(2) + ee(6))*delint(5) + ee48*delint(7)
 ak(12) = d(1)*delint(8) + ee48*delint(10)
 ak(13) = 0.0D0
 ak(14) = ee48*delint(7)
 ak(15) = (ee(3)+ee(7))*delint(7)
 ak(16) = (ee(3)+ee(7))*delint(10) + ee48*delint(8)
 ak(17) = ak( 3)
 ak(18) = ak(11)
 ak(19) = ee(6)*delint(3) + ee(16)*delint(7) + (ee(8) + ee(14))*delint(5)
 ak(20) = (ee(2) + ee(6))*delint(6) + ee(16)*delint(10) + (ee(8) +  &
     ee(13) + ee(14))*delint(8)
 ak(21) = 0.0D0
 ak(22) = ee(16)*delint(7) + ee(8)*delint(5)
 ak(23) = ee( 7)*delint(5) + ee(15)*delint(7)
 ak(24) = (ee(7) + ee(16))*delint(8) + ee(8)*delint(6) + ee(15)*delint(10)
 ak(25) = ak( 4)
 ak(26) = ak(12)
 ak(27) = ak(20)
 ak(28) = d(1)*delint(9) + ee(16)*delint(12) + (ee48 + ee(13) +  &
     ee(14))*delint(11)
 ak(29) = 0.0D0
 ak(30) = ee(16)*delint(10) + ee48*delint(8)
 ak(31) = (ee(3) + ee(7))*delint(8) + ee(15)*delint(10)
 ak(32) = (ee(3) + ee(7) + ee(16))*delint(11) + ee(15)*delint(12) +  &
     ee48*delint(9)
 ak(33) = 0.0D0
 ak(34) = 0.0D0
 ak(35) = 0.0D0
 ak(36) = 0.0D0
 ak(37) = 0.0D0
 ak(38) = 0.0D0
 ak(39) = 0.0D0
 ak(40) = 0.0D0
 ak(41) = ak( 6)
 ak(42) = ak(14)
 ak(43) = ak(22)
 ak(44) = ak(30)
 ak(45) = 0.0D0
 ak(46) = ee(16)*delint(7)
 ak(47) = ee(15)*delint(7)
 ak(48) = ee(16)*delint(8) + ee(15)*delint(10)
 ak(49) = ak( 7)
 ak(50) = ak(15)
 ak(51) = ak(23)
 ak(52) = ak(31)
 ak(53) = 0.0D0
 ak(54) = ak(47)
 ak(55) = ee(11)*delint( 7)
 ak(56) = ee(11)*delint(10) + ee(12)*delint(8)
 ak(57) = ak( 8)
 ak(58) = ak(16)
 ak(59) = ak(24)
 ak(60) = ak(32)
 ak(61) = 0.0D0
 ak(62) = ak(48)
 ak(63) = ak(56)
 ak(64) = ee(11)*delint(12) + ee(16)*delint(9) + (ee(12) + ee(15))*delint(11)
 
 DO  i = 1,64
   ak(i) = twopi*ak(i)
 END DO
 
!     TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM FIELD COORDINATES
!     TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmatd (gambq,8,8,1, ak,8,8,0, d)
 CALL gmmatd (d,8,8,0, gambq,8,8,0, ak)
 
!     ZERO OUT THE (6X6) MATRIX USED AS INPUT TO THE INSERTION ROUTINE
 
 DO  i = 1,36
   aki(i) = 0.0D0
 END DO
 
!     LOCATE THE TRANSFORMATION MATRICES FOR THE FOUR  GRID POINTS
 
 DO  i = 1,4
   IF (ics(i) == 0) CYCLE
   k = 9*(i-1) + 1
   CALL transd (ics(i),d(k))
 END DO
 
!     START THE LOOP FOR INSERTION OF THE FOUR  (6X6) MATRICES
!     INTO THE MASTER STIFFNESS MATRIX
 
 ir1  = 2*ipp - 1
 iapp = 9*(ipp-1) + 1
 DO  i = 1,4
   
!     PLACE THE APPROIATE (2X2) SUBMATRIX OF THE STIFFNESS MATRIX
!     IN A (3X3) MATRIX FOR TRANSFORMATION
   
   ic1    = 2*i - 1
   irc    = (ir1-1)*8 + ic1
   akt(1) = ak(irc)
   akt(2) = 0.0D0
   akt(3) = ak(irc+1)
   akt(4) = 0.0D0
   akt(5) = 0.0D0
   akt(6) = 0.0D0
   akt(7) = ak(irc+8)
   akt(8) = 0.0D0
   akt(9) = ak(irc+9)
   
!     TRANSFORM THE (3X3) STIFFNESS MATRIX
   
   IF (ics(ipp) == 0) GO TO 820
   CALL gmmatd (d(iapp),3,3,1, akt(1),3,3,0, d(37))
   DO  j = 1,9
     akt(j) = d(j+36)
   END DO
   820 CONTINUE
   IF (ics(i) == 0) GO TO 840
   iai = 9*(i-1) + 1
   CALL gmmatd (akt(1),3,3,0, d(iai),3,3,0, d(37))
   DO  j = 1,9
     akt(j) = d(j+36)
   END DO
   840 CONTINUE
   
!     PLACE THE TRANSFORMED (3X3) MATRIX INTO A (6X6) MATRIX FOR
!     THE INSERTION ROUTINE
   
   j = 0
   DO  j1 = 1,18,6
     DO  j2 = 1,3
       j = j + 1
       k = j1 + j2 - 1
       aki(k) = akt(j)
     END DO
   END DO
   
!     CALL THE INSERTION ROUTINE
   
   CALL sma1b (aki(1),igp(i),-1,ifkgg,0.0D0)
   IF (iopt4 == 0 .OR. gsube == 0.0) CYCLE
   k4ggsw = 1
   dampc  = gsube
   CALL sma1b (aki(1),igp(i),-1,if4gg,dampc)
 END DO
 RETURN
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 910 i = 37
 GO TO 950
 920 i = 26
 GO TO 950
 930 i = 221
! ...     221 WILL PRINT USER MESSAGE 2218
 
 950 IF (nmsg == 0) GO TO 970
 IF (nmsg >= mmsg) RETURN
 DO  j = 1,nmsg
   IF (msg(3,j) == idel .AND. msg(2,j) == i) RETURN
 END DO
 970 ics(1) = idel
 ics(2) = irg
 CALL mesage (30,i,ics)
 nogo = 1
 RETURN
 
END SUBROUTINE ktrapr
