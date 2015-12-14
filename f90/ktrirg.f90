SUBROUTINE ktrirg
     
 
!*****
! THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR A AXI-SYMMETRIC RING
! WITH A TRIANGULAR CROSS SECTION
!*****
 
 
!                        ECPT FOR THE TRIANGULAR RING
 
 
!                                                      TYPE
! ECPT( 1) ELEMENT IDENTIFICATION                        I
! ECPT( 2) SCALAR INDEX NO. FOR GRID POINT A             I
! ECPT( 3) SCALAR INDEX NO. FOR GRID POINT B             I
! ECPT( 4) SCALAR INDEX NO. FOR GRID POINT C             I
! ECPT( 5) MATERIAL ORIENTATION ANGLE(DEGREES)           R
! ECPT( 6) MATERIAL IDENTIFICATION                       I
! ECPT( 7) COOR. SYS. ID. FOR GRID POINT A               I
! ECPT( 8) X-COOR. OF GRID POINT A (IN BASIC COOR.)      R
! ECPT( 9) Y-COOR. OF GRID POINT A (IN BASIC COOR.)      R
! ECPT(10) Z-COOR. OF GRID POINT A (IN BASIC COOR.)      R
! ECPT(11) COOR. SYS. ID. FOR GRID POINT B               I
! ECPT(12) X-COOR. OF GRID POINT B (IN BASIC COOR.)      R
! ECPT(13) Y-COOR. OF GRID POINT B (IN BASIC COOR.)      R
! ECPT(14) Z-COOR. OF GRID POINT B (IN BASIC COOR.)      R
! ECPT(15) COOR. SYS. ID. FOR GRID POINT C               I
! ECPT(16) X-COOR. OF GRID POINT C (IN BASIC COOR.)      R
! ECPT(17) Y-COOR. OF GRID POINT C (IN BASIC COOR.)      R
! ECPT(18) Z-COOR. OF GRID POINT C (IN BASIC COOR.)      R
! ECPT(19) EL. TEMPERATURE FOR MATERIAL PROPERTIES       R
 
 
 DOUBLE PRECISION :: constd, degrad
 DOUBLE PRECISION :: d ,  gambq,    r,   z  &
     ,                  teo, ee,  delint,   ak,  aki ,                  akt
 DOUBLE PRECISION :: r1,  r2,  r3,  z1,  z2,  z3,  zmin, dgama  &
     ,                  dr,  rh,  dz,  zh,  ra,  za,  area  &
     ,                  er,  et,  ez,  vrt, vtr, vtz, vzt  &
     ,                  vzr, vrz, grz, del, cosg,sing,dgamr  &
     ,                  twopi,    dki
 DOUBLE PRECISION :: dampc
 
 DIMENSION          iecpt(19)
 DIMENSION          aki(36),  akt(9)
 
 COMMON /condad/ constd(5)
 COMMON   /sma1io/ dum1(10)  &
     ,                  ifkgg ,                  igkgg,    if4gg,    dum2(21)
 COMMON   /sma1cl/ iopt4,    k4ggsw  &
     ,                  npvt ,                   dum4(7)  &
     ,                  link(10)           ,idetck  &
     ,                  dodet              ,nogo
 COMMON   /sma1et/ ecpt(19)  &
     ,                  dum5(81)
 COMMON   /matin/ matidc             ,matflg  &
     ,                  eltemp             ,stress  &
     ,                  sinth              ,costh
 COMMON   /matout/ e(3)               ,anu(3)  &
     ,                  rho                ,g(3)  &
     ,                  alf(3)        ,tzero,    gsube
 COMMON   /sma1dp/ d(36) ,   gambq(36),     r(3) ,    z(3)  &
     ,                  teo(16),  ee(16),   delint(8),     ak(36)  &
     ,                  dgama,    zmin  &
     ,                  dr,  rh,  dz,  zh,  ra,  za,  area  &
     ,                  er,  et,  ez,  vrt, vtr, vtz, vzt  &
     ,                  vzr, vrz, grz, del, cosg,sing,dgamr  &
     ,                  igp(3) , ics(3) , sp(18) ,                  tempe
 
 EQUIVALENCE ( constd(2) , twopi  )
 EQUIVALENCE ( constd(4) , degrad )
 EQUIVALENCE        (iecpt(1) , ecpt(1))
 EQUIVALENCE   (r(1),r1),     (r(2),r2),     (r(3),r3)  &
     ,             (z(1),z1),     (z(2),z2),     (z(3),z3)
 EQUIVALENCE        (aki(1),  gambq(1))
 EQUIVALENCE        (akt(1),  teo(1))
 
! ----------------------------------------------------------------------
 
! STORE ECPT PARAMETERS IN LOCAL VARIABLES
 
 idel  = iecpt(1)
 igp(1)= iecpt(2)
 igp(2)= iecpt(3)
 igp(3)= iecpt(4)
 matid = iecpt(6)
 ics(1)= iecpt(7)
 ics(2)= iecpt(11)
 ics(3)= iecpt(15)
 r(1)  = ecpt(8)
 d(1)  = ecpt(9)
 z(1)  = ecpt(10)
 r(2)  = ecpt(12)
 d(2)  = ecpt(13)
 z(2)  = ecpt(14)
 r(3)  = ecpt(16)
 d(3)  = ecpt(17)
 z(3)  = ecpt(18)
 tempe = ecpt(19)
 dgama = ecpt(5)
 
 
! CHECK INTERNAL GRID POINTS FOR PIVOT POINT
 
 ipp = 0
 DO  i = 1,3
   IF (npvt == igp(i)) ipp = i
 END DO
 IF (ipp == 0) CALL mesage (-30,34,idel)
 
 
! TEST THE VALIDITY OF THE GRID POINT COORDINATES
 
 ieror1 = 0
 DO  i = 1,3
   IF (r(i) > 0.0D0) CYCLE
   IF (ieror1 /= 0) CYCLE
   CALL mesage (30, 211, idel)
   ieror1 = 1
 END DO
 ieror2 = 0
 DO  i = 1, 3
   IF (d(i) == 0.0D0) CYCLE
   IF (ieror2 /= 0) CYCLE
   CALL mesage (30, 212, idel)
   ieror2 = 1
 END DO
 IF (ieror1 == 0.AND.ieror2 == 0) GO TO 220
 nogo = 2
 RETURN
 220 IF ((r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1) < 0.0D0) GO TO 920
 
 
! COMPUTE THE ELEMENT COORDINATES
 
 zmin = DMIN1(z1, z2, z3)
 z1 = z1 - zmin
 z2 = z2 - zmin
 z3 = z3 - zmin
 
 
 
! FORM THE TRANSFORMATION MATRIX (6X6) FROM FIELD COORDINATES TO GRID
! POINT DEGREES OF FREEDOM
 
 DO  i = 1,36
   gambq(i) = 0.0D0
 END DO
 gambq( 1) = 1.0D0
 gambq( 2) = r1
 gambq( 3) = z1
 gambq(10) = 1.0D0
 gambq(11) = r1
 gambq(12) = z1
 gambq(13) = 1.0D0
 gambq(14) = r2
 gambq(15) = z2
 gambq(22) = 1.0D0
 gambq(23) = r2
 gambq(24) = z2
 gambq(25) = 1.0D0
 gambq(26) = r3
 gambq(27) = z3
 gambq(34) = 1.0D0
 gambq(35) = r3
 gambq(36) = z3
 
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL inverd (6, gambq(1),6 , d(10), 0, d(11) , ising , sp)
 
 IF (ising == 2) GO TO 920
 
 
 
! CALCULATE THE INTEGRAL VALUES IN ARRAY DELINT WHERE THE ORDER IS
! INDICATED BY THE FOLLOWING TABLE
 
!              DELINT( 1) - (-1,0)
!              DELINT( 2) - (-1,1)
!              DELINT( 3) - (-1,2)
!              DELINT( 4) - ( 0,0)
!              DELINT( 5) - ( 0,1)
!              DELINT( 6) - ( 1,0)
!              DELINT( 7) - ( 0,2)
!              DELINT( 8) - ( 1,2)
 
 
! TEST FOR RELATIVE SMALL AREA OF INTEGRATION
! AND IF AREA IS SMALL THEN APPROXIMATE INTEGRALS
 
 dr = DMAX1 ( DABS(r1-r2) , DABS(r2-r3) , DABS(r3-r1) )
 rh = DMIN1 ( r1 , r2 , r3 ) / 10.0D0
 dz = DMAX1 ( DABS(z1-z2) , DABS(z2-z3) , DABS(z3-z1) )
 zh = DMIN1 ( z1 , z2 , z3 ) / 10.0D0
 ra = (r1 + r2 + r3) / 3.0D0
 za = (z1 + z2 + z3) / 3.0D0
 area =(r1*(z2-z3) + r2*(z3-z1) + r3*(z1-z2)) / 2.0D0
 kode = 0
 IF (DABS( (r2-r1)/r2 ) < 1.0D-5) kode = 1
 IF ( dr <= rh  .OR.  dz <= zh ) kode = -1
 
 
 310 CONTINUE
 i1 = 0
 DO  i = 1,3
   ip = i - 2
   DO  j = 1,3
     iq = j - 1
     IF (ip == 1 .AND. iq == 1) CYCLE
     i1 = i1 + 1
     IF (kode < 0) THEN
       GO TO   320
     ELSE IF (kode == 0) THEN
       GO TO   330
     ELSE
       GO TO   340
     END IF
     320 delint(i1) =((ra) ** ip)*((za) ** iq) * area
     CYCLE
     330 delint(i1) =   dki (1,3,1,2,1,3,ip,iq,r,z)  &
         +  dki (3,2,1,2,3,2,ip,iq,r,z)
     CYCLE
     340 CONTINUE
     delint(i1) =   dki (1,3,3,2,1,3,ip,iq,r,z)
   END DO
 END DO
 d(1)      = delint(6)
 delint(6) = delint(7)
 delint(7) = d(1)
 
 
! TEST FOR EXCESSIVE ROUND-OFF ERROR IN INTEGRAL CALCULATIONS
! AND IF IT EXIST APPROXIMATE INTEGRALS
 
 IF (kode < 0) GO TO 500
 DO  i = 1,8
   IF (delint(i) < 0.0D0) GO TO 475
 END DO
 IF (delint(8) <= delint(7)) GO TO 475
 IF (delint(3) >= delint(8)) GO TO 475
 IF (delint(3) > delint(7)) GO TO 475
 GO TO 500
 475 CONTINUE
 kode = -1
 GO TO 310
 500 CONTINUE
 
 
 
! LOCATE THE MATERIAL PROPERTIES IN THE MAT1 OR MAT3 TABLE
 
 matidc = matid
 matflg = 7
 eltemp = tempe
 CALL  mat (idel)
 
 
! SET MATERIAL PROPERTIES IN DOUBLE PRECISION VARIABLES
 
 er = e(1)
 et = e(2)
 ez = e(3)
 vrt = anu(1)
 vtz = anu(2)
 vzr = anu(3)
 grz = g(3)
 vtr = vrt * et / er
 vzt = vtz * ez / et
 vrz = vzr * er / ez
 del = 1.0D0 - vrt*vtr - vtz*vzt - vzr*vrz - vrt*vtz*vzr - vrz*vtr*vzt
 
 
! GENERATE ELASTIC CONSTANTS MATRIX (4X4)
 
 ee(1) = er * (1.0D0 - vtz*vzt) / del
 ee(2) = er * (vtr + vzr*vtz) / del
 ee(3) = er * (vzr + vtr*vzt) / del
 ee(4) = 0.0D0
 ee(5) = ee(2)
 ee(6) = et * (1.0D0 - vrz*vzr) / del
 ee(7) = et * (vzt + vrt*vzr) / del
 ee(8) = 0.0D0
 ee(9) = ee(3)
 ee(10)= ee(7)
 ee(11)= ez * (1.0D0 - vrt*vtr) / del
 ee(12)= 0.0D0
 ee(13)= 0.0D0
 ee(14)= 0.0D0
 ee(15)= 0.0D0
 ee(16)= grz
 
 
! FORM TRANSFORMATION MATRIX (4X4) FROM MATERIAL AXIS TO ELEMENT
! GEOMETRIC AXIS
 
 dgamr = dgama * degrad
 cosg = DCOS(dgamr)
 sing = DSIN(dgamr)
 teo( 1) = cosg ** 2
 teo( 2) = 0.0D0
 teo( 3) = sing ** 2
 teo( 4) = sing * cosg
 teo( 5) = 0.0D0
 teo( 6) = 1.0D0
 teo( 7) = 0.0D0
 teo( 8) = 0.0D0
 teo( 9) = teo(3)
 teo(10) = 0.0D0
 teo(11) = teo(1)
 teo(12) = -teo(4)
 teo(13) = -2.0D0 * teo(4)
 teo(14) = 0.0D0
 teo(15) = -teo(13)
 teo(16) = teo(1) - teo(3)
 
 
! TRANSFORM THE ELASTIC CONSTANTS MATRIX FROM MATERIAL
! TO ELEMENT GEOMETRIC AXIS
 
 CALL gmmatd (teo , 4, 4, 1, ee , 4, 4, 0, d )
 CALL gmmatd (d   , 4, 4, 0, teo, 4, 4, 0, ee)
 
 
 
! FORM THE ELEMENT STIFFNESS MATRIX IN FIELD COORDINATES
 
 ak( 1) = ee(6) * delint(1)
 ak( 2) = (ee(2) + ee(6)) * delint(4)
 ak( 3) = ee(6) * delint(2) + ee(8) * delint(4)
 ak( 4) = 0.0D0
 ak( 5) = ee(8) * delint(4)
 ak( 6) = ee(7) * delint(4)
 ak( 7) = ak(2)
 ak( 8) = (ee(1) + 2.0D0*ee(2) + ee(6)) * delint(6)
 ak( 9) = (ee(2) + ee(6)) * delint(5) + (ee(4) + ee(8)) *delint(6)
 ak(10) = 0.0D0
 ak(11) = (ee(4) + ee(8)) * delint(6)
 ak(12) = (ee(3) + ee(7)) * delint(6)
 ak(13) = ak(3)
 ak(14) = ak(9)
 ak(15) = ee(6) * delint(3) + 2.0D0*ee(8) * delint(5) + ee(16) * delint(6)
 ak(16) = 0.0D0
 ak(17) = ee(8) * delint(5) + ee(16) * delint(6)
 ak(18) = ee(7) * delint(5) + ee(12) * delint(6)
 ak(19) = 0.0D0
 ak(20) = 0.0D0
 ak(21) = 0.0D0
 ak(22) = 0.0D0
 ak(23) = 0.0D0
 ak(24) = 0.0D0
 ak(25) = ak(5)
 ak(26) = ak(11)
 ak(27) = ak(17)
 ak(28) = 0.0D0
 ak(29) = ee(16) * delint(6)
 ak(30) = ee(12) * delint(6)
 ak(31) = ak(6)
 ak(32) = ak(12)
 ak(33) = ak(18)
 ak(34) = 0.0D0
 ak(35) = ak(30)
 ak(36) = ee(11) * delint(6)
 
 DO  i = 1,36
   ak(i) = twopi * ak(i)
 END DO
 
! TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM FIELD COORDINATES
! TO GRID POINT DEGREES OF FREEDOM
 
 CALL gmmatd (gambq , 6, 6, 1, ak , 6, 6, 0, d )
 CALL gmmatd (d  , 6, 6, 0, gambq , 6, 6, 0, ak)
 
 
 
! ZERO OUT THE (6X6) MATRIX USED AS INPUT TO THE INSERTION ROUTINE
 
 DO  i = 1,36
   aki(i) = 0.0D0
 END DO
 
 
! LOCATE THE TRANSFORMATION MATRICES FOR THE THREE GRID POINTS
 
 DO  i = 1,3
   IF (ics(i) == 0) CYCLE
   k = 9 * (i-1) + 1
   CALL transd (ics(i) , d(k))
 END DO
 
 
 
! START THE LOOP FOR INSERTION OF THE THREE (6X6) MATRICES
! INTO THE MASTER STIFFNESS MATRIX
 
 ir1  = 2 * ipp - 1
 iapp = 9 * (ipp-1) + 1
 DO  i = 1,3
   
! PLACE THE APPROIATE (2X2) SUBMATRIX OF THE STIFFNESS MATRIX
! IN A (3X3) MATRIX FOR TRANSFORMATION
   
   ic1 = 2 * i - 1
   irc = (ir1 - 1) * 6 + ic1
   akt(1) = ak(irc)
   akt(2) = 0.0D0
   akt(3) = ak(irc+1)
   akt(4) = 0.0D0
   akt(5) = 0.0D0
   akt(6) = 0.0D0
   akt(7) = ak(irc+6)
   akt(8) = 0.0D0
   akt(9) = ak(irc+7)
   
! TRANSFORM THE (3X3) STIFFNESS MATRIX
   
   IF (ics(ipp) == 0) GO TO 820
   CALL gmmatd (d(iapp) , 3, 3, 1, akt(1) , 3, 3, 0, d(28) )
   DO  j = 1,9
     akt(j) = d(j+27)
   END DO
   820 CONTINUE
   IF (ics(i) == 0) GO TO 840
   iai = 9 * (i - 1) + 1
   CALL gmmatd (akt(1) , 3, 3, 0, d(iai) , 3, 3, 0, d(28) )
   DO  j = 1,9
     akt(j) = d(j+27)
   END DO
   840 CONTINUE
   
! PLACE THE TRANSFORMED (3X3) MATRIX INTO A (6X6) MATRIX FOR
! THE INSERTION ROUTINE
   
   j = 0
   DO  j1 = 1,18,6
     DO  j2 = 1,3
       j = j + 1
       k = j1 + j2 - 1
       aki(k) = akt(j)
     END DO
   END DO
   
! CALL THE INSERTION ROUTINE
   
   CALL sma1b (aki(1) , igp(i), -1, ifkgg, 0.0D0)
   IF (iopt4 == 0 .OR. gsube == 0.0) CYCLE
   k4ggsw = 1
   dampc = gsube
   CALL sma1b (aki(1) , igp(i) , -1,if4gg , dampc )
 END DO
 RETURN
 
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
 915 nogo=1
 RETURN
 920 CALL mesage(30,26,idel)
 GO TO 915
 
END SUBROUTINE ktrirg
