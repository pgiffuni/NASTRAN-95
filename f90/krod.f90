SUBROUTINE krod
!*****
! THIS ROUTINE COMPUTES THE TWO 6 X 6 MATRICES  K(NPVT,NPVT) AND
! K(NPVT,J) FOR A ROD HAVING END POINTS NUMBERED NPVT AND J.
!*****
 
!                        E C P T  F O R  T H E  R O D
 
!                                                                CARD
!                                                 TYPE   TABLE   TYPE
! ECPT( 1)ELEMENT ID.                               I     ECT    CROD
! ECPT( 2)SCALAR INDEX NUMBER FOR GRID POINT A      I     ECT    CROD
! ECPT( 3)SCALAR INDEX NUMBER FOR GRID POINT B      I     ECT    CROD
! ECPT( 4)MATERIAL ID.                              I     EPT    PROD
! ECPT( 5)AREA  (A)                                 R     EPT    PROD
! ECPT( 6)POLAR MOMENT OF INERTIA (J)               R     EPT    PROD
! ECPT( 7) TORSIONAL STRESS COEFF (C)                R    EPT    PROD
! ECPT( 8) NON-STRUCTRAL MASS (MU)                   R    EPT    PROD
! ECPT( 9) COOR. SYS. ID. NO. FOR GRID POINT A       I   BGPDT   GRID
! ECPT(10) X-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(11) Y-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(12) Z-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(13) COOR. SYS. ID. NO. FOR GRID POINT B       I   BGPDT
! ECPT(14) X-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(15) Y-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(16) Z-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(17) ELEMENT TEMPERATURE
 
 LOGICAL :: heat
 
 DOUBLE PRECISION :: x                  ,y  &
     ,                  z                  ,xl  &
     ,                  xn                 ,dscl  &
     ,                  dscr               ,dampc  &
     ,                  d                  ,ke  &
     ,                  ti                 ,dumdp
 
 DIMENSION iecpt(4)
 COMMON /BLANK/icom
 COMMON   /system/ isys
 
! SMA1 I/O PARAMETERS
 
 COMMON   /sma1io/ ifcstm             ,ifmpt  &
     ,                  ifdit              ,idum1  &
     ,                  ifecpt             ,igecpt  &
     ,                  ifgpct             ,iggpct  &
     ,                  ifgei              ,iggei  &
     ,                  ifkgg              ,igkgg  &
     ,                  if4gg              ,ig4gg  &
     ,                  ifgpst             ,iggpst  &
     ,                  inrw               ,outrw  &
     ,                  clsnrw             ,clsrw  &
     ,                  neor               ,eor  &
     ,                  mcbkgg(7)          ,mcb4gg(7)
 
! SMA1 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON   /sma1bk/ icstm              ,ncstm  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6k              ,n6x6k  &
     ,                  i6x64              ,n6x64
 
! SMA1 PROGRAM CONTROL PARAMETERS
 
 COMMON   /sma1cl/ iopt4              ,k4ggsw  &
     ,                  npvt               ,left  &
     ,                  frowic             ,lrowic  &
     ,                  nrowsc             ,tnrows  &
     ,                  jmax               ,nlinks  &
     ,                  link(10)           ,idetck  &
     ,                  dodet              ,nogo
 COMMON/sma1ht/     heat
 
! ECPT COMMON BLOCK
 
 COMMON   /sma1et/ ecpt(17)           ,dumet(83)
 
! INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON   /matin/ matidc             ,matflg  &
     ,                  eltemp             ,stress  &
     ,                  sinth              ,costh
 COMMON   /matout/ e                  ,g  &
     ,                  nu                 ,rho  &
     ,                  alpha              ,tsubo  &
     ,                  gsube              ,sigt  &
     ,                  sigc               ,sigs
 COMMON/hmtout/     fk
 
! LOCAL DOUBLE PRECISION VARIABLES
 
 COMMON   /sma1dp/ x                  ,y  &
     ,                  z                  ,xl  &
     ,                  xn(3)              ,dscl  &
     ,                  dscr               ,dampc  &
     ,                  d(18)              ,ke(36)  &
     ,                  ti(9)              ,dumdp(227)
 
 
! NOTE THAT EQUIVALENCE IS NECESSARY SINCE ECPT IS A MIXED --- INTEGERS
! AND REAL --- ARRAY
 
 EQUIVALENCE (iecpt(1),ecpt(1))
!*****
!  BRANCH ON HEAT FORMULATION.
!*****
 IF( heat ) GO TO 200
 IF (iecpt(2) == npvt) GO TO 10
 IF (iecpt(3) /= npvt) CALL mesage (-30,34,iecpt(1))
 itemp = iecpt(2)
 iecpt(2) = iecpt(3)
 iecpt(3) = itemp
 ka  = 13
 kb  =  9
 GO TO 20
 10 ka  =  9
 kb  =  13
 
! AT THIS POINT KA POINTS TO THE COOR. SYS. ID. OF THE PIVOT GRID POINT.
! SIMILARLY FOR KB AND THE NON-PIVOT GRID POINT.
! NOW COMPUTE THE LENGTH OF THE ROD.
 
! WE STORE THE COORDINATES IN THE D ARRAY SO THAT ALL ARITHMETIC WILL BE
! DOUBLE PRECISION
 
 20 d(1) = ecpt(ka+1)
 d(2) = ecpt(ka+2)
 d(3) = ecpt(ka+3)
 d(4) = ecpt(kb+1)
 d(5) = ecpt(kb+2)
 d(6) = ecpt(kb+3)
 x    = d(1) - d(4)
 y    = d(2) - d(5)
 z    = d(3) - d(6)
 xl = DSQRT (x**2 + y**2 + z**2)
 IF (xl /= 0.0D0) GO TO 30
 CALL mesage(30,26,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULA
 
 nogo=1
 RETURN
 30 CONTINUE
 
! CALCULATE A NORMALIZED DIRECTION VECTOR IN BASIC COORDINATES.
 
 xn(1) = x / xl
 xn(2) = y / xl
 xn(3) = z / xl
 
! LOCATE E = YOUNG-S MODULUS, G = SHEAR MODULUS AND DAMPC = DAMPING
! CONSTANT IN THE MAT1 TABLE AND COMPUTE DSCL = A * E / XL AND
! DSCR = J * G / XL.  A IS ECPT(5) AND J IS ECPT(6)
 
 matidc = iecpt(4)
 matflg = 1
 eltemp = ecpt(17)
 CALL mat (iecpt(1))
 
! WE STORE ECPT(5), ECPT(6), E AND G IN DOUBLE PRECISION LOCATIONS SO
! THAT ALL ARITHMETIC WILL BE DOUBLE PRECISION
 
 d(1) = ecpt(5)
 d(2) = e
 d(3) = ecpt(6)
 d(4) = g
 dscl = d(1) * d(2) / xl
 dscr = d(3) * d(4) / xl
 dampc  = g_sub_e
 
! SET UP THE -N- MATRIX AND STORE AT D(1)
 
 d(1) = xn(1) * xn(1)
 d(2) = xn(1) * xn(2)
 d(3) = xn(1) * xn(3)
 d(4) = d(2)
 d(5) = xn(2) * xn(2)
 d(6) = xn(2) * xn(3)
 d(7) = d(3)
 d(8) = d(6)
 d(9) = xn(3) * xn(3)
 
! ZERO OUT THE 6X6 WHICH WILL BE USED FOR STORAGE OF KGG(NPVT,NONPVT),
! NONPVT = NPVT,J
! KGG(NPVT,NONPVT), NONPVT = NPVT,J
 
 DO  i = 1,36
   ke(i) = 0.0D0
 END DO
 nonpvt = 2
 k2 = 1
 
! IF PIVOT GRID POINT IS IN BASIC COORDINATES, GO TO 70
 
 IF (iecpt(ka) == 0) GO TO 70
 CALL transd (ecpt(ka),ti(1))
 CALL gmmatd (ti(1),3,3,1, d(1),3,3,0, d(10))
 CALL gmmatd (d(10),3,3,0, ti(1),3,3,0, d(1))
 
! AT THIS POINT D(1) CONTAINS THE MATRIX PRODUCT TAT * N * TA
! AND D(10) CONTAINS THE MATRIX PRODUCT TAT * N.
 
 ASSIGN 100 TO iretrn
 GO TO 80
 70 ASSIGN 90 TO iretrn
 
! FILL THE KE MATRIX
 
 80 ke( 1) = dscl * d(k2  )
 ke( 2) = dscl * d(k2+1)
 ke( 3) = dscl * d(k2+2)
 ke( 7) = dscl * d(k2+3)
 ke( 8) = dscl * d(k2+4)
 ke( 9) = dscl * d(k2+5)
 ke(13) = dscl * d(k2+6)
 ke(14) = dscl * d(k2+7)
 ke(15) = dscl * d(k2+8)
 ke(22) = dscr * d(k2  )
 ke(23) = dscr * d(k2+1)
 ke(24) = dscr * d(k2+2)
 ke(28) = dscr * d(k2+3)
 ke(29) = dscr * d(k2+4)
 ke(30) = dscr * d(k2+5)
 ke(34) = dscr * d(k2+6)
 ke(35) = dscr * d(k2+7)
 ke(36) = dscr * d(k2+8)
 CALL sma1b (ke,ecpt(nonpvt),-1,ifkgg,0.0D0)
 IF (iopt4 == 0  .OR.  g_sub_e == 0.0) GO TO 85
 k4ggsw = 1
 CALL sma1b (ke,ecpt(nonpvt),-1,if4gg,dampc)
 
!  RETURN  FROM  FILL  CODE W/ IRETRN = 90 IMPLIES G.P. A WAS IN BASIC
!    .      .     .      .      .     =100 IMPLIES G.P. A WAS NOT BASIC
!    .      .     .      .      .     =140 IMPLIES THE K(NPVT,NONPVT)
!                                      HAS BEEN COMPUTED AND INSERTED
!                                      AND HENCE WE ARE FINISHED.
 
 85 GO TO iretrn , (90,100,140)
 90 k1 = 1
 k2 = 10
 GO TO 110
 100 k1 = 10
 k2 = 1
 110 nonpvt = 3
 
! IF NON-PIVOT GRID POINT IS IN BASIC COORDINATES, GO TO 120
 
 IF (iecpt(kb) == 0) GO TO 120
 CALL transd (ecpt(kb),ti(1))
 
! RECALL THAT D(K1) CONTAINS TAT * N.
 
 CALL gmmatd (d(k1),3,3,0, ti(1),3,3,0, d(k2))
 
! AT THIS POINT D(K2) CONTAINS TAT * N * TB.
 
 GO TO 130
 120 k2 = k1
 130 ASSIGN 140 TO iretrn
 
! SET CONSTANTS NEGATIVE TO PROPERLY COMPUTE K(NPVT,NONPVT)
 
 dscr = -dscr
 dscl = -dscl
 GO TO 80
 
! A TRANSFER TO STATEMENT NO. 140 IMPLIES KGG AND/OR K4GG CALCULATIONS
! HAVE BEEN COMPLETED.
 
 140 RETURN
!*****
!  HEAT FORMULATION.  FIRST COMPUTE LENGTH OF ELEMENT.
!*****
 200 x = ecpt(14) - ecpt(10)
 y = ecpt(15) - ecpt(11)
 z = ecpt(16) - ecpt(12)
 xl= DSQRT(x**2 + y**2 + z**2)
 IF( xl  > 0.0) THEN
   GO TO   400
 END IF
 300 CALL mesage( -30, 26, iecpt(1) )
 
!     GET MATERIAL PROPERTY -K- FROM HMAT ROUTINE
 
 400 matflg = 1
 matidc = iecpt(4)
 eltemp = ecpt(17)
 CALL hmat( iecpt )
 
 xl = DBLE(fk) * DBLE(ecpt(5)) / xl
 
 IF( npvt == iecpt(3) ) xl = -xl
 DO  i = 1,2
   CALL sma1b( xl, iecpt(i+1), npvt, ifkgg, 0.0D0 )
   xl = -xl
 END DO
 RETURN
END SUBROUTINE krod
