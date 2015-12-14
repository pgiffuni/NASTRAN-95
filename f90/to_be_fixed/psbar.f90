SUBROUTINE psbar
!*****
! THIS ROUTINE COMPUTES THE TWO 6 X 6 MATRICES K(NPVT,NPVT) AND
! K(NPVT,J) FOR A BAR ELEMENT HAVING END POINTS NUMBERED NPVT AND J.
!*****
 
!                     E C P T  F O R  T H E  B A R
 
! ECPT( 1)  -  IELID          ELEMENT ID. NUMBER
! ECPT( 2)  -  ISILNO(2)      * SCALAR INDEX NOS. OF THE GRID POINTS
! ECPT( 3)  -    ...          *
! ECPT( 4)  -  SMALLV(3)      $ REFERENCE VECTOR
! ECPT( 5)  -    ...          $
! ECPT( 6)  -    ...          $
! ECPT( 7)  -  ICSSV          COOR. SYS. ID FOR SMALLV VECTOR
! ECPT( 8)  -  IPINFL(2)      * PIN FLAGS
! ECPT( 9)  -    ...          *
! ECPT(10)  -  ZA(3)          $ OFFSET VECTOR FOR POINT A
! ECPT(11)  -    ...          $
! ECPT(12)  -    ...          $
! ECPT(13)  -  ZB(3)          * OFFSET VECTOR FOR POINT B
! ECPT(14)  -    ...          *
! ECPT(15)  -    ...          *
! ECPT(16)  -  IMATID         MATERIAL ID.
! ECPT(17)  -  A              CROSS-SECTIONAL AREA
! ECPT(18)  -  I1             $ AREA MOMENTS OF INERTIA
! ECPT(19)  -  I2             $
! ECPT(20)  -  FJ             POLAR MOMENT OF INERTIA
! ECPT(21)  -  NSM            NON-STRUCTURAL MASS
! ECPT(22)  -  FE             FORCE ELEMENT DESCRIPTIONS (FORCE METHOD)
! ECPT(23)  -  C1             * STRESS RECOVERY COEFFICIENTS
! ECPT(24)  -  C2             *
! ECPT(25)  -  D1             *
! ECPT(26)  -  D2             *
! ECPT(27)  -  F1             *
! ECPT(28)  -  F2             *
! ECPT(29)  -  G1             *
! ECPT(30)  -  G2             *
! ECPT(31)  -  K1             $ AREA FACTORS FOR SHEAR
! ECPT(32)  -  K2             $
! ECPT(33)  -  I12            AREA MOMENT OF INERTIA
! ECPT(34)  -  MCSIDA         COOR. SYS. ID. FOR GRID POINT A
! ECPT(35)  -  GPA(3)         * BASIC COORDINATES FOR GRID POINT A
! ECPT(36)  -    ...          *
! ECPT(37)  -    ...          *
! ECPT(38)  -  MCSIDB         COOR. SYS. ID. FOR GRID POINT B
! ECPT(39)  -  GPB(3)         $ BASIC COORDINATES FOR GRID POINT B
! ECPT(40)  -    ...          $
! ECPT(41)  -    ...          $
! ECPT(42)  -  ELTEMP         AVG. ELEMENT TEMPERATURE
! ECPT(43)  -  EPSIN1         PREVIOUS STRAIN VALUE ONCE REMOVED
! ECPT(44)  -  EPSIN2         PREVIOUS STRAIN VALUE
! ECPT(45)  -  ESTAR          PREVIOUSLY COMPUTED MODULUS OF ELASTICITY
! ECPT(46)  -  V1STAR         * ELEMENT FORCES, INITIALLY ZERO
! ECPT(47)  -   V2STAR        *
! ECPT(48)  -   TSTAR         *
! ECPT(49)  -   M1ASTR        *
! ECPT(50)  -   M2ASTR        *
! ECPT(51)  -   UAIN(6)       $ INCREMENTAL DISPLACEMENT VECTOR AT PT. A
! ECPT(52)  -   ...           $
! ECPT(53)  -   ...           $
! ECPT(54)  -   ...           $
! ECPT(55)  -   ...           $
! ECPT(56)  -   ...           $
! ECPT(57)  -   UBIN(6)       * INCREMENTAL DISPLACEMENT VECTOR AT PT. B
! ECPT(58)  -   ...           *
! ECPT(59)  -   ...           *
! ECPT(60)  -   ...           *
! ECPT(61)  -   ...           *
! ECPT(62)  -   ...           *
 
 REAL :: l                  ,lsq  &
     ,                  lcube              ,i1  &
     ,                  i2                 ,k1  &
     ,                  k2                 ,ke  &
     ,                  kep                ,i12  &
     ,                  nsm                ,lr1  &
     ,                  lr2                ,lb  &
     ,                  l2b3               ,l2b6
 REAL :: m1astr             ,m2astr  &
     ,                  m1a                ,m2a  &
     ,                  m1b                ,m2b  &
     ,                  k1a                ,k2a  &
     ,                  k1b                ,k2b
 LOGICAL :: abasic             ,bbasic  &
     ,                  basic              ,aofset  &
     ,                  bofset             ,offset
 DIMENSION veci(3)            ,vecj(3)  &
     ,                  veck(3)            ,ecpt(100)  &
     ,                  iecpt(100)         ,ipin(10)  &
     ,                  ta(18)             ,tb(9)  &
     ,                  smalv0(6)          ,u(24)  &
     ,                  d(9)               ,sa(72)  &
     ,                  sb(36)             ,fa(6) ,                  fb(6)
 
 COMMON   /pla32e/ ielid              ,isilno(2)  &
     ,                  smallv(3)          ,icssv  &
     ,                  ipinfl(2)          ,za(3)  &
     ,                  zb(3)              ,imatid  &
     ,                  a                  ,i1  &
     ,                  i2                 ,fj  &
     ,                  nsm                ,fe  &
     ,                  c1                 ,c2  &
     ,                  d1                 ,d2  &
     ,                  f1                 ,f2  &
     ,                  g1                 ,g2  &
     ,                  k1                 ,k2 ,                  i12  &
     ,                  mcsida             ,gpa(3)  &
     ,                  mcsidb             ,gpb(3)  &
     ,                  eltemp             ,epsin1  &
     ,                  epsin2             ,estar  &
     ,                  v1star             ,v2star
 COMMON   /pla32e/ tstar              ,m1astr  &
     ,                  m2astr             ,uain(6) ,                  ubin(6)
 COMMON   /pla32s/ ke(144)            ,kep(144)  &
     ,                  dela(6)            ,delb(6)
 COMMON   /pla32c/ gamma              ,gammas
 COMMON   /sout/ iselid             ,sig1a  &
     ,                  sig2a              ,sig3a  &
     ,                  sig4a              ,sigax  &
     ,                  sigamx             ,sigamn  &
     ,                  msten              ,sig1b  &
     ,                  sig2b              ,sig3b  &
     ,                  sig4b              ,sigbmx  &
     ,                  sigbmn             ,mscom ,                  dum14(14)
 COMMON   /matin/ matidc             ,matflg  &
     ,                  temdum             ,plaarg ,                  dum2(2)
 COMMON   /matout/ e sub 0            ,g sub 0  &
     ,                  nu                 ,rho  &
     ,                  alpha              ,t_sub_0  &
     ,                  gsube              ,sigmat  &
     ,                  sigmac             ,sigmas
 
 EQUIVALENCE (ielid,ecpt(1),iecpt(1))  &
     ,                  (ta(10),tb(1)) ,                  (ecpt(71),d(1))  &
     ,                  (e sub 0,plaans) ,                  (sa(37),sb(1))  &
     ,                  (msten,smten) ,                  (mscom,smcom)
 
!-----------------------------------------------------------------------
 
! SET UP POINTERS TO COOR. SYS. IDS., OFFSET VECTORS, AND PIN FLAGS.
! ICSIDA AND ICSIDB ARE COOR. SYS. IDS.
 
 jcsida = 34
 jcsidb = 38
 jofsta = 10
 jofstb = 13
 jpina  =  8
 jpinb  =  9
 icsida = iecpt(34)
 icsidb = iecpt(38)
 
! NORMALIZE THE REFERENCE VECTOR WHICH LIES IN THE FIRST PRINCIPAL AXIS
! PLANE  (FMMS - 36 P. 4)
 
 fl = 0.0
 DO  i = 1,3
   fl = fl + smallv(i)**2
 END DO
 fl = SQRT(fl)
 DO  i = 1,3
   smallv(i) = smallv(i) / fl
 END DO
 
! DETERMINE IF POINT A AND B ARE IN BASIC COORDINATES OR NOT.
 
 abasic = .true.
 bbasic = .true.
 IF (icsida /= 0) abasic = .false.
 IF (icsidb /= 0) bbasic = .false.
 
! COMPUTE THE TRANSFORMATION MATRICES TA AND TB IF NECESSARY
 
 IF (.NOT. abasic) CALL transs (ecpt(jcsida),ta)
 IF (.NOT. bbasic) CALL transs (ecpt(jcsida),tb)
 
! DETERMINE IF WE HAVE NON-ZERO OFFSET VECTORS.
 
 aofset = .true.
 j = jofsta - 1
 DO  i = 1,3
   j = j + 1
   IF (ecpt(j) /= 0.0) GO TO 80
 END DO
 aofset = .false.
 80 bofset = .true.
 j = jofstb - 1
 DO  i = 1,3
   j = j + 1
   IF (ecpt(j) /= 0.0) GO TO 100
 END DO
 bofset = .false.
 
! FORM THE CENTER AXIS OF THE BEAM WITHOUT OFFSETS.
 
 100 veci(1) = ecpt(jcsida+1) - ecpt(jcsidb+1)
 veci(2) = ecpt(jcsida+2) - ecpt(jcsidb+2)
 veci(3) = ecpt(jcsida+3) - ecpt(jcsidb+3)
 
! TRANSFORM THE OFFSET VECTORS IF NECESSARY
 
 IF ( .NOT. aofset  .AND.  .NOT. bofset )  GO TO 150
 
! TRANSFORM THE OFFSET VECTOR FOR POINT A IF NECESSARY.
 
 idela = 1
 j = jofsta - 1
 DO  i = 1,3
   j = j + 1
   dela(i) = ecpt(j)
 END DO
 IF (abasic) GO TO 120
 idela = 4
 CALL gmmats (ta,3,3,0, dela(1),3,1,0, dela(4) )
 
! TRANSFORM THE OFFSET VECTOR FOR POINT B IF NECESSARY
 
 120 idelb = 1
 j = jofstb - 1
 DO  i = 1,3
   j = j + 1
   delb(i) = ecpt(j)
 END DO
 IF (bbasic) GO TO 140
 idelb = 4
 CALL gmmats (tb,3,3,0, delb(1),3,1,0, delb(4) )
 
! SINCE THERE WAS AT LEAST ONE NON-ZERO OFFSET VECTOR RECOMPUTE VECI
 
 140 veci(1) = veci(1) + dela(idela  ) - delb(idelb  )
 veci(2) = veci(2) + dela(idela+1) - delb(idelb+1)
 veci(3) = veci(3) + dela(idela+2) - delb(idelb+2)
 
! COMPUTE THE LENGTH OF THE BIG V (VECI) VECTOR AND NORMALIZE
 
 150 veci(1) = -veci(1)
 veci(2) = -veci(2)
 veci(3) = -veci(3)
 fl = SQRT (veci(1)**2  +  veci(2)**2  +  veci(3)**2)
 DO  i = 1,3
   veci(i) = veci(i) / fl
 END DO
 
! COMPUTE THE SMALL V SUB 0 VECTOR, SMALV0.  ****CHECK THIS LOGIC****
 
 DO  i = 1,3
   smalv0(i) = smallv(i)
 END DO
 isv = 1
 IF (icssv == 0) GO TO 180
 isv = 4
 CALL gmmats (ta,3,3,0, smalv0(1),3,1,0, smalv0(4) )
 
! COMPUTE THE K VECTOR, VECK = VECI  X  SMALV0, AND NORMALIZE
 
 180 veck(1) =  veci(2) * smalv0(isv+2)  -  veci(3) * smalv0(isv+1)
 veck(2) =  veci(3) * smalv0(isv  )  -  veci(1) * smalv0(isv+2)
 veck(3) =  veci(1) * smalv0(isv+1)  -  veci(2) * smalv0(isv)
 fll     =  SQRT ( veck(1)**2  +  veck(2)**2  +  veck(3)**2  )
 veck(1) =  veck(1) / fll
 veck(2) =  veck(2) / fll
 veck(3) =  veck(3) / fll
 
! COMPUTE THE J VECTOR, VECJ = VECK  X  VECI, AND NORMALIZE
 
 vecj(1) =  veck(2) * veci(3)  -  veck(3) * veci(2)
 vecj(2) =  veck(3) * veci(1)  -  veck(1) * veci(3)
 vecj(3) =  veck(1) * veci(2)  -  veck(2) * veci(1)
 fll     =  SQRT ( vecj(1)**2  +  vecj(2)**2  +  vecj(3)**2 )
 vecj(1) =  vecj(1) / fll
 vecj(2) =  vecj(2) / fll
 vecj(3) =  vecj(3) / fll
 
! SET UP INTERMEDIATE VARIABLES FOR ELEMENT STIFFNESS MATRIX CALCULATION
 
 l = fl
 lsq = l**2
 lcube = lsq * l
 
! STORE INCREMENTAL DISPLACEMENT VECTORS IN DOUBLE PRECISION LOCATIONS
 
 DO  i = 1,6
   u(i)    = uain(i)
   u(i+12) = ubin(i)
 END DO
!*****
! COMPUTE ON FIRST PASS C  * E  * U   AND C  * E  * U  ON SECOND PASS
!                        B    B    B       A    A    A
!*****
 ipass = 1
 basic = bbasic
 offset = bofset
 jofset = jofstb
 jcsid = 10
 INDEX = 13
 
! IF THERE ARE OFFSETS FOR THIS POINT, CONSTRUCT THE 3 X 3 MATRIX D.
 
 184 IF (.NOT. offset) GO TO 188
 d(1) = 0.0
 d(2) =  ecpt(jofset+2)
 d(3) = -ecpt(jofset+1)
 d(4) = -d(2)
 d(5) = 0.0
 d(6) =  ecpt(jofset)
 d(7) = -d(3)
 d(8) = -d(6)
 d(9) = 0.0
 
! COMPUTE THE 3 VECTOR  D * U , WHERE U  IS THE VECTOR OF THE 3
!                            R         R
! ROTATIONAL DISPLACEMENTS
 
 CALL gmmats (d,3,3,0, u(INDEX+3),3,1,0, u(INDEX+6))
 
! ADD OFFSET CONTRIBUTION TO THE TRANSLATION COMPONENTS OF THE DISPLACE-
! MENT VECTOR
 
 j = INDEX
 DO  i = 1,3
   u(j) = u(j) + u(j+6)
   j = j + 1
 END DO
 
! TRANSFORM TRANSLATIONAL COMPONENTS TO BASIC COORDINATES IF NECESSARY
 
 188 IF (basic) GO TO 190
 CALL gmmats (ta(jcsid),3,3,0, u(INDEX),3,1,0, u(INDEX+3) )
 
! STORE TRANSFORMED VECTOR BACK INTO ITS ORIGINAL D.P. LOCATION
 
 u(INDEX  ) = u(INDEX+3)
 u(INDEX+1) = u(INDEX+4)
 u(INDEX+2) = u(INDEX+5)
 190 IF (ipass == 2) GO TO 192
 ipass = 2
 basic  = abasic
 offset = aofset
 jofset = jofsta
 jcsid  = 1
 INDEX  = 1
 GO TO 184
 
! FORM THE DIFFERENCE OF THE TRANSLATIONAL COMPONENTS OF THE TRANSFORMED
! DISPLACEMENT VECTORS
 
 192 DO  i = 1,3
   u(i+12) = u(i+12) - u(i)
 END DO
 
! FORM DOT PRODUCT
 
 CALL gmmats (veci,3,1,1, u(13),3,1,0, d(1) )
 
! CALCULATE THE INCREMENTAL ELEMENT STRAIN
 
 deps1 = d(1) / l
 
! PERFORM EXTENSIONAL STRAIN CALCULATIONS
 
 deps2  = epsin2 - epsin1
 eps1   = epsin2 + deps1
 eps2 = epsin2 + (deps1 + gammas**2 * deps2)* (gamma + 1.0E0)  &
     /(gammas + 1.0E0) + gammas * (deps1 - gammas*deps2) * (gamma+1.0E0)**2  &
     /(gammas + 1.0E0)
 
! CALL MAT ROUTINE TO GET SIGMA1 AND SIGMA2 AS FUNCTIONS OF EPS1,EPS2
 
 matidc = imatid
 matflg = 1
 CALL mat (iecpt(1))
 e sub 0 l = e sub 0
 g sub 0 l = g sub 0
 matflg = 6
 plaarg = eps1
 CALL mat (iecpt(1))
 sigma1 = plaans
 plaarg = eps2
 CALL mat (iecpt(1))
 sigma2 = plaans
 
! NOTE THAT E1 IS USED IN THIS ROUTINE ONLY TO UPDATE THE EST (ECPT)
! ENTRY
 
 IF (eps1 == eps2) GO TO 200
 e1 = (sigma2 - sigma1) / (eps2 - eps1)
 GO TO 202
 200 e1 = estar
 
! BEGIN ELEMENT STRESS MATRIX CALCULATIONS.
 
 202 e = estar
 g = estar * g sub 0 l / e sub 0 l
 ei1  = e * i1
 ei2  = e * i2
 IF (k1 == 0.0  .OR.  i12 /= 0.0) GO TO 210
 gak1 = g * a * k1
 r1 = (12.0 * ei1 * gak1) / (gak1 * lcube + 12.0 * l * ei1)
 GO TO 220
 210 r1 =  12.0 * ei1 / lcube
 220 IF (k2 == 0.0  .OR.  i12 /= 0.0) GO TO 230
 gak2 = g * a * k2
 r2 = (12.0 * ei2 * gak2) / (gak2 * lcube + 12.0 * l * ei2)
 GO TO 240
 230 r2 =  12.0 * ei2 / lcube
 
! COMPUTE THE -SMALL- K-S, SK1, SK2, SK3 AND SK4
 
 240 sk1 = 0.25 * r1 * lsq  +  ei1 / l
 sk2 = 0.25 * r2 * lsq  +  ei2 / l
 sk3 = 0.25 * r1 * lsq  -  ei1 / l
 sk4 = 0.25 * r2 * lsq  -  ei2 / l
 
! COMPUTE THE TERMS THAT WILL BE NEEDED FOR THE 12 X 12 MATRIX KE
 
 ael = a * e / l
 lr1 = l * r1 / 2.0
 lr2 = l * r2 / 2.0
 gjl = g * fj / l
 
! CONSTRUCT THE 12 X 12 MATRIX KE
 
 DO  i = 1,144
   ke(i) = 0.0
 END DO
 ke(  1) =  ael
 ke(  7) = -ael
 ke( 14) =  r1
 ke( 18) =  lr1
 ke( 20) = -r1
 ke( 24) =  lr1
 ke( 27) =  r2
 ke( 29) = -lr2
 ke( 33) = -r2
 ke( 35) = -lr2
 ke( 40) =  gjl
 ke( 46) = -gjl
 ke( 51) = -lr2
 ke( 53) =  sk2
 ke( 57) =  lr2
 ke( 59) =  sk4
 ke( 62) =  lr1
 ke( 66) =  sk1
 ke( 68) = -lr1
 ke( 72) =  sk3
 ke( 73) = -ael
 ke( 79) =  ael
 ke( 86) = -r1
 ke( 90) = -lr1
 ke( 92) =  r1
 ke( 96) = -lr1
 ke( 99) = -r2
 ke(101) =  lr2
 ke(105) =  r2
 ke(107) =  lr2
 ke(112) = -gjl
 ke(118) =  gjl
 ke(123) = -lr2
 ke(125) =  sk4
 ke(129) =  lr2
 ke(131) =  sk2
 ke(134) =  lr1
 ke(138) =  sk3
 ke(140) = -lr1
 ke(144) =  sk1
 IF (i12 == 0.0) GO TO 255
 beta = 12.0 * e * i12 / lcube
 lb   = l * beta / 2.0
 l2b3 = lsq * beta / 3.0
 l2b6 = lsq * beta / 6.0
 ke( 15) =  beta
 ke( 17) = -lb
 ke( 21) = -beta
 ke( 23) = -lb
 ke( 26) =  beta
 ke( 30) =  lb
 ke( 32) = -beta
 ke( 36) =  lb
 ke( 50) = -lb
 ke( 54) = -l2b3
 ke( 56) =  lb
 ke( 60) = -l2b6
 ke( 63) =  lb
 ke( 65) = -l2b3
 ke( 69) = -lb
 ke( 71) = -l2b6
 ke( 87) = -beta
 ke( 89) =  lb
 ke( 93) =  beta
 ke( 95) =  lb
 ke( 98) = -beta
 ke(102) = -lb
 ke(104) =  beta
 ke(108) = -lb
 ke(122) = -lb
 ke(126) = -l2b6
 ke(128) =  lb
 ke(132) = -l2b3
 ke(135) =  lb
 ke(137) = -l2b6
 ke(141) = -lb
 ke(143) = -l2b3
 
! DETERMINE IF THERE ARE NON-ZERO PIN FLAGS.
 
 255 ka = iecpt(jpina)
 kb = iecpt(jpinb)
 IF (ka == 0  .AND.  kb == 0) GO TO 325
 
! SET UP THE IPIN ARRAY
 
 DO  i = 1,5
   ipin(i  ) = MOD(ka,10)
   ipin(i+5) = MOD(kb,10) + 6
   IF (ipin(i+5) == 6) ipin(i+5) = 0
   ka = ka / 10
   kb = kb / 10
 END DO
 
! ALTER KE MATRIX DUE TO PIN FLAGS.
 
 DO  i = 1,10
   IF (ipin(i) == 0) CYCLE
   ii = 13 * ipin(i)  -  12
   IF (ke(ii) /= 0.0) GO TO 280
   il = ipin(i)
   ii = ii - il
   DO  j = 1,12
     ii = ii + 1
     ke(ii) = 0.0
     ke(il) = 0.0
     il = il + 12
   END DO
   CYCLE
   280 DO  j = 1,12
     ji  = 12 * (j-1) + ipin(i)
     ij = 12 * (ipin(i) - 1) + j
     DO  ll = 1,12
       jll = 12 * (j-1) + ll
       ill = 12 * (ipin(i) - 1) + ll
       kep(jll) = ke(jll) - (ke(ill)/ke(ii)) * ke(ji)
     END DO
     kep(ij) = 0.0
     kep(ji) = 0.0
   END DO
   DO  k = 1,144
     ke(k) = kep(k)
   END DO
 END DO
 
!        E
! STORE K   IN KEP(1),...,KEP(36) AND
!        AA
 
!        E
! STORE K   IN KEP(37),...,KEP(72)
!        AB
 
 325 j = 0
 DO  i = 1,72,12
   low = i
   lim = low + 5
   DO  k = low,lim
     j = j + 1
     kep(j) = ke(k)
     kep(j+36) = ke(k+6)
   END DO
 END DO
 
!                                                        T
! STORE VECI, VECJ, VECK IN KE(1),...,KE(9) FORMING THE A  MATRIX.
 
 ke(1) = veci(1)
 ke(2) = veci(2)
 ke(3) = veci(3)
 ke(4) = vecj(1)
 ke(5) = vecj(2)
 ke(6) = vecj(3)
 ke(7) = veck(1)
 ke(8) = veck(2)
 ke(9) = veck(3)
 
! SET POINTERS SO THAT WE WILL BE WORKING WITH POINT A.
 
 basic = abasic
 jcsid  = jcsida
 offset = aofset
 jofset = jofsta
 iwbeg  = 0
 ikel   = 1
 iab = 1
 INDEX = isilno(1)
 
! ZERO OUT THE ARRAY WHERE THE 3 X 3 MATRIX AND THE W  AND W  6 X 6
! MATRICES WILL RESIDE.                              A      B
 
 DO  i = 28,108
   ke(i) = 0.0
 END DO
 
! SET UP THE -G- MATRIX.  IG POINTS TO THE BEGINNING OF THE G MATRIX.
! G = AT X TI
 
 360 ig = 1
 IF (basic) GO TO 370
 CALL transs (ecpt(jcsid),ke(10))
 CALL gmmats (ke(1),3,3,0, ke(10),3,3,0, ke(19) )
 ig = 19
 
! IF THERE IS A NON-ZERO OFFSET FOR THE POINT, SET UP THE D 3 X3 MATRIX.
 
 370 IF ( .NOT. offset ) GO TO 380
 ke(10) = 0.0
 ke(11) =  ecpt(jofset+2)
 ke(12) = -ecpt(jofset+1)
 ke(13) = -ke(11)
 ke(14) = 0.0
 ke(15) =  ecpt(jofset)
 ke(16) = -ke(12)
 ke(17) = -ke(15)
 ke(18) = 0.0
 
! FORM THE 3 X 3 PRODUCT H = G X D, I.E., KE(28) = KE(IG) X KE(10)
 
 CALL gmmats (ke(ig),3,3,0, ke(10),3,3,0, ke(28) )
 
 
! FORM THE W  MATRIX OR THE W  MATRIX IN KE(37) OR KE(73) DEPENDING
!           A                B
! UPON WHICH POINT - A OR B - IS UNDER CONSIDERATION.  G WILL BE STORED
! IN THE UPPER LEFT AND LOWER RIGHT CORNERS.  H, IF NON-ZERO, WILL BE
! STORED IN THE UPPER RIGHT CORNER.
 
 
 380 ke(iwbeg+37) = ke(ig  )
 ke(iwbeg+38) = ke(ig+1)
 ke(iwbeg+39) = ke(ig+2)
 ke(iwbeg+43) = ke(ig+3)
 ke(iwbeg+44) = ke(ig+4)
 ke(iwbeg+45) = ke(ig+5)
 ke(iwbeg+49) = ke(ig+6)
 ke(iwbeg+50) = ke(ig+7)
 ke(iwbeg+51) = ke(ig+8)
 ke(iwbeg+58) = ke(ig  )
 ke(iwbeg+59) = ke(ig+1)
 ke(iwbeg+60) = ke(ig+2)
 ke(iwbeg+64) = ke(ig+3)
 ke(iwbeg+65) = ke(ig+4)
 ke(iwbeg+66) = ke(ig+5)
 ke(iwbeg+70) = ke(ig+6)
 ke(iwbeg+71) = ke(ig+7)
 ke(iwbeg+72) = ke(ig+8)
 IF ( .NOT. offset ) GO TO 390
 ke(iwbeg+40) = ke(28)
 ke(iwbeg+41) = ke(29)
 ke(iwbeg+42) = ke(30)
 ke(iwbeg+46) = ke(31)
 ke(iwbeg+47) = ke(32)
 ke(iwbeg+48) = ke(33)
 ke(iwbeg+52) = ke(34)
 ke(iwbeg+53) = ke(35)
 ke(iwbeg+54) = ke(36)
 
!                          E                      E
! FORM THE PRODUCT  S  =  K    X  W   OR  S   = K    X  W  , DEPENDING
!                    A     AA      A       B     AB      B
! UPON WHICH POINT WE ARE WORKING WITH.
 
 390 CALL gmmats (kep(ikel),6,6,0, ke(iwbeg+37),6,6,0, sa(iab) )
 
! IF THE POINT UNDER CONSIDERATION IS POINT B WE ARE FINISHED.  IF NOT,
! SET UP POINTS AND INDICATORS FOR WORKING WITH POINT B.
 
 IF (iwbeg == 36) GO TO 500
 basic = bbasic
 jcsid = jcsidb
 offset = bofset
 jofset = jofstb
 iwbeg  = 36
 ikel   = 37
 iab = 37
 INDEX  = isilno(2)
 DO  i = 28,36
   ke(i) = 0.0
 END DO
 GO TO 360
 
! COMPUTE FORCES AND MOMENTS FROM  S   AND   S   AND DISPLACEMENT
!                                   A         B
! VECTORS
 
 500 CALL gmmats (sa,6,6,0, uain,6,1,0, fa)
 CALL gmmats (sb,6,6,0, ubin,6,1,0, fb)
 fx = a * sigma1
 v1  = -fa(2) - fb(2) + v1star
 v2  = -fa(3) - fb(3) + v2star
 t   = -fa(4) - fb(4) + tstar
 m2a =  fa(5) + fb(5) + m2astr
 m1a = -fa(6) - fb(6) + m1astr
 m1b = m1a - v1*l
 m2b = m2a - v2*l
!*****
! COMPUTE ELEMENT STRESSES AT 4 POINTS
!*****
 
! COMPUTE K1A AND K2A
 
 IF (i12 /= 0.0) GO TO 530
 IF (i1 /= 0.0) GO TO 520
 k1a = 0.0
 GO TO 540
 520 k1a = -m1a / i1
 GO TO 540
 530 k1a = (m2a * i12  -  m1a * i2) / (i1 * i2  -  i12**2)
 k2a = (m1a * i12  -  m2a * i1) / (i1 * i2  -  i12**2)
 GO TO 560
 540 IF (i2 /= 0.0) GO TO 550
 k2a = 0.0
 GO TO 560
 550 k2a = -m2a / i2
 
! COMPUTE SIG1A, SIG2A, SIG3A AND SIG4A
 
 560 sig1a = k1a * c1  +  k2a * c2
 sig2a = k1a * d1  +  k2a * d2
 sig3a = k1a * f1  +  k2a * f2
 sig4a = k1a * g1  +  k2a * g2
 
! COMPUTE K1B AND K2B
 
 IF (i12 /= 0.0) GO TO 580
 IF (i1 /= 0.0) GO TO 570
 k1b = 0.0
 GO TO 590
 570 k1b = -m1b / i1
 GO TO 590
 580 k1b = (m2b * i12  -  m1b * i2) / (i1 * i2  -  i12**2)
 k2b = (m1b * i12  -  m2b * i1) / (i1 * i2  -  i12**2)
 GO TO 610
 590 IF (i2 /= 0.0) GO TO 600
 k2b = 0.0
 GO TO 610
 600 k2b = -m2b / i2
 
! COMPUTE SIG1B, SIG2B, SIG3B AND SIG4B
 
 610 sig1b = k1b * c1  +  k2b * c2
 sig2b = k1b * d1  +  k2b * d2
 sig3b = k1b * f1  +  k2b * f2
 sig4b = k1b * g1  +  k2b * g2
 
! COMPUTE AXIAL STRESS
 
 sigax = 0.0
 IF (a /= 0.0) sigax = fx / a
 
! COMPUTE MAXIMA AND MINIMA
 
 sigamx = sigax + AMAX1(sig1a,sig2a,sig3a,sig4a)
 sigbmx = sigax + AMAX1(sig1b,sig2b,sig3b,sig4b)
 sigamn = sigax + AMIN1(sig1a,sig2a,sig3a,sig4a)
 sigbmn = sigax + AMIN1(sig1b,sig2b,sig3b,sig4b)
 
! COMPUTE MARGIN OF SAFETY IN TENSION
 
 IF(sigmat <= 0.0)GO TO 620
 IF(AMAX1(sigamx,sigbmx) <= 0.0) GO TO 620
 q=sigmat/AMAX1(sigamx,sigbmx)
 smten=q-1.0
 GO TO 630
 620 msten=1
 
!      COMPUTE MARGIN OF SAFETY IN COMPRESSION
 
 630 sigmac=-ABS(sigmac)
 IF(AMIN1(sigamn,sigbmn) >= 0.0) GO TO 640
 w=sigmac/AMIN1(sigamn,sigbmn)
 smcom=w-1.0
 GO TO 650
 640 mscom=1
 650 iselid = ielid
 
! UPDATE EST (ECPT) ENTRIES
 
 epsin1 = epsin2
 epsin2 = eps1
 estar  = e1
 v1star = v1
 v2star = v2
 tstar  = t
 m1astr = m1a
 m2astr = m2a
 RETURN
END SUBROUTINE psbar
