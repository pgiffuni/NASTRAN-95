SUBROUTINE etrbks (iopt)
     
!     THIS SUBROUTINE CALCULATES THE STIFFNESS MATRIX FOR THE BASIC
!     BENDING TRIANGLE. IT IS USED BY SUBROUTINES TRBSCS, QDPLTS,
!     TRPLTS, QUAD1S,TRIA1S,TRIA2S
!     SINGLE PRECISION VERSION
!     IOPT MAY BE VARIED AS FOLLOWS TO PRODUCE APPROPRIATE RESULTS
!     ******************************************************************
 
 
 
!     IOPT = 0   IMPLIES DO COMPLETE BASIC BENDING TRIANGLE.
!     IOPT = 1   IMPLIES COMPUTE ONLY THE NINE (3X3)MATRICES
!                WHICH FORM THE 9X9 K SUPER U - MATRIX.
!     IOPT = 2   SAME AS IOPT = 1,BUT SAVE H-INVERSE AND S...
 
!     ******************************************************************
!     ECPT LIST FOR BASIC BENDING TRIANGLE             NAME IN
!                                                      THIS
!     ECPT                                             ROUTINE   TYPE
!     ******************************************************************
!     ECPT( 1) = ELEMENT ID                            NECPT(1)  INTEGER
!     ECPT( 2) = GRID POINT A                          NGRID(1)  INTEGER
!     ECPT( 3) = GRID POINT B                          NGRID(2)  INTEGER
!     ECPT( 4) = GRID POINT C                          NGRID(3)  INTEGER
!     ECPT( 5) = THETA = ANGLE OF MATERIAL             ANGLE     REAL
!     ECPT( 6) = MATERIAL ID 1                         MATID1    INTEGER
!     ECPT( 7) = I = MOMENT OF INERTIA                 EYE       REAL
!     ECPT( 8) = MATERIAL ID 2                         MATID2    INTEGER
!     ECPT( 9) = T2                                    T2        REAL
!     ECPT(10) = NON-STRUCTURAL-MASS                   FMU       REAL
!     ECPT(11) = Z1                                    Z11       REAL
!     ECPT(12) = Z2                                    Z22       REAL
!     ECPT(13) = COORD. SYSTEM ID 1                    NECPT(13) INTEGER
!     ECPT(14) = X1                                    X1        REAL
!     ECPT(15) = Y1                                    Y1        REAL
!     ECPT(16) = Z1                                    Z1        REAL
!     ECPT(17) = COORD. SYSTEM ID 2                    NECPT(17) INTEGER
!     ECPT(18) = X2                                    X2        REAL
!     ECPT(19) = Y2                                    Y2        REAL
!     ECPT(20) = Z2                                    Z2        REAL
!     ECPT(21) = COORD. SYSTEM ID 3                    NECPT(21) INTEGER
!     ECPT(22) = X3                                    X3        REAL
!     ECPT(23) = Y3                                    Y3        REAL
!     ECPT(24) = Z3                                    Z3        REAL
!     ECPT(25) = ELEMENT TEMPERATURE                   ELTEMP    REAL
!     ******************************************************************
 
 
 INTEGER, INTENT(IN)                      :: iopt
 LOGICAL :: nogo
 INTEGER :: necpt(26),no(2)
 REAL :: j2x2,k
 DIMENSION       d(9),g2x2(4),j2x2(4),s(18),ecpt(25),g(9),  &
     tjte(18),tite(18),ti(9)
 
 COMMON /condas/ consts(5)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     t_sub_0, g sub e,sigten,sigcom,sigshe, g2x211,g2x212,g2x222
 
 COMMON /emgprm/ dum(15),smb(3),iprec,nogo
 COMMON /emgest/ ielid         ,ngrid(3) ,angle         ,matid1  &
     ,eye           ,matid2 ,t2            ,fmu  &
     ,z11           ,z22 ,dummy1        ,x1  &
     ,y1            ,z1 ,dummy2        ,x2  &
     ,y2            ,z2 ,dummy3        ,x3  &
     ,y3            ,z3 ,dumb(76)
 COMMON /emgtrx/ a(225),prod9(9),temp9(9),xsubb,xsubc,ysubc,dict5,  &
     e(18),k(324),aout(324)
 
 EQUIVALENCE    (consts(4),degra)  ,(g2x2(1),a(88))
 EQUIVALENCE    (ielid,ecpt(1),necpt(1)), (d(1),g(1),a(79)),  &
     (tjte(1),a(100))   ,(tite(1),s(1),a(82)), (j2x2(1),a( 92))   ,(ti(1),a(118))
 
 DATA no /81,190/
 
 ntype = 0
 IF( iopt > 0 ) ntype = 1
 IF (ntype == 1) GO TO 100
 eltemp = ecpt(25)
!     SET UP  I, J, K VECTORS STORING AS FOLLOWS AND ALSO CALCULATE
!     X-SUB-B, X-SUB-C, AND Y-SUB-C.
 
!     E(11), E(14), E(17) WILL BE THE I-VECTOR.
!     E(12), E(15), E(18) WILL BE THE J-VECTOR.
!     E( 1), E( 4), E( 7) WILL BE THE K-VECTOR.
 
!     FIND I-VECTOR = RSUBB - RUBA (NON-NORMALIZED)
 e(11) =  x2-x1
 e(14) =  y2-y1
 e(17) =  z2-z1
 
!     FIND LENGTH = X-SUB-B COOR. IN ELEMENT SYSTEM
 xsubb =  SQRT(e(11)**2 + e(14)**2 + e(17)**2)
 IF  (xsubb <= 1.e-6)  GO TO 7770
 
 e(11) =  e(11)/xsubb
 e(14) =  e(14)/xsubb
 e(17) =  e(17)/xsubb
 
!     TAKE RSUBC - RSUBA AND STORE TEMPORARILY IN E(2), E(5), E(8)
 
 e(2) =  x3-x1
 e(5) =  y3-y1
 e(8) =  z3-z1
 
!     X-SUB-C = I . (RSUBC - RSUBA), THUS
 
 xsubc = e(11) * e(2)  +  e(14) * e(5)  +  e(17) * e(8)
 
!     CROSSING I-VECTOR TO (RSUBC - RSUBA) GIVES THE K-VECTOR
!     (NON-NORMALIZED)
 
 e(1) = e(14) * e( 8)  -  e( 5) * e(17)
 e(4) = e( 2) * e(17)  -  e(11) * e( 8)
 e(7) = e(11) * e( 5)  -  e( 2) * e(14)
 
!     FIND LENGTH = Y-SUB-C COOR. IN ELEMENT SYSTEM
 ysubc = SQRT(e(1)**2 + e(4)**2 + e(7)**2)
 IF (ysubc <= 1.e-6) GO TO 7780
 
!     NORMALIZE K-VECTOR WITH Y-SUB-C
 
 e(1) = e(1)/ysubc
 e(4) = e(4) / ysubc
 e(7) = e(7) / ysubc
 
!     NOW HAVING I AND K VECTORS GET -- J = K CROSS I
 
 e(12) = e( 4) * e(17)  -  e(14) * e( 7)
 e(15) = e(11) * e( 7)  -  e( 1) * e(17)
 e(18) = e( 1) * e(14)  -  e(11) * e( 4)
 
!     NORMALIZE J-VECTOR FOR COMPUTER EXACTNESS JUST TO MAKE SURE
 
 temp = SQRT(e(12)**2 + e(15)**2   + e(18)**2)
 e(12) = e(12) / temp
 e(15) = e(15) / temp
 e(18) = e(18) / temp
 e(2) = 0.
 e(3) = 0.
 e(5) = 0.
 e(6) = 0.
 e(8) = 0.
 e(9) = 0.
 e(10) = 0.
 e(13) = 0.
 e(16) = 0.
 
!     CONVERT ANGLE FROM DEGREES TO RADIANS STORING IN THETA.
 
 theta = angle * degra
 sinth = SIN( theta )
 costh = COS( theta )
 IF(ABS(sinth) < 1.0E-06) sinth = 0.0E0
 
!     ******************************************************************
 
!     SETTING UP G MATRIX
 100 inflag = 2
 matid = matid1
 CALL mat( ecpt(1) )
 
!     FILL G-MATRIX WITH OUTPUT FROM MAT ROUTINE
 
 g(1) = g11
 g(2) = g12
 g(3) = g13
 g(4) = g12
 g(5) = g22
 g(6) = g23
 g(7) = g13
 g(8) = g23
 g(9) = g33
 
!     ******************************************************************
 
!  50 COMPUTATION OF D = I.G-MATRIX (EYE IS INPUT FROM THE ECPT)
 
 DO  i=1,9
   d(i) = g(i)*eye
 END DO
 
!     ******************************************************************
 
 area = xsubb * ysubc/2.
 xbar = (xsubb+xsubc)/3.
 ybar = ysubc/3.
 
 xcsq = xsubc ** 2
 ycsq = ysubc ** 2
 xbsq = xsubb ** 2
 xcyc = xsubc * ysubc
 px2 =(xbsq+xsubb*xsubc + xcsq)/6.
 py2  =  ycsq/6.
 pxy2 =   ysubc*(xsubb + 2.*xsubc)/12.
 xbar3 =  3.* xbar
 ybar3 =  3.* ybar
 ybar2 =  2.* ybar
 
!     ******************************************************************
!                 X
!     FILL THE  (K ) MATRIX STORING IN  A(1). . .A(36)
 
 a( 1) = d( 1)
 a( 2) = d( 3)
 a( 3) = d( 2)
 a( 4) = d( 1) * xbar3
 a( 5) = d( 2) * xbar   +  ybar2 * d(3)
 a( 6) = d( 2) * ybar3
 a( 7) = a( 2)
 a( 8) = d( 9)
 a( 9) = d( 6)
 a(10) = d( 3) * xbar3
 a(11) = d( 6) * xbar   +  ybar2 * d(9)
 a(12) = d( 6) * ybar3
 a(13) = a( 3)
 a(14) = a( 9)
 a(15) = d( 5)
 a(16) = d( 2) * xbar3
 a(17) = d( 5) * xbar   +  ybar2 * d(6)
 a(18) = d( 5) * ybar3
 a(19) = a( 4)
 a(20) = a(10)
 a(21) = a(16)
 a(22) = d(1)*9.*px2
 a(23) = d(2)*3.*px2  +  6.*pxy2*d(3)
 a(24) = d(2)*9.*pxy2
 a(25) = a( 5)
 a(26) = a(11)
 a(27) = a(17)
 a(28) = a(23)
 a(29) =  d(5)*px2 + 4.*pxy2*d(6)   +  4.*py2*d(9)
 a(30) =  d(5)*3.*pxy2  +  6.*py2*d(6)
 a(31) = a( 6)
 a(32) = a(12)
 a(33) = a(18)
 a(34) = a(24)
 a(35) = a(30)
 a(36) =  d(5)*9.*py2
 temp =  4.*area
 DO   i=1,36
   a(i) = a(i)*temp
 END DO
 
!     ******************************************************************
 
!     F1LL  (HBAR) MATRIX STORING AT A(37). . .A(72)
 
 DO  i =37,72
   a(i)=0.
 END DO
 
 a(37) = xbsq
 a(40) = xbsq * xsubb
 a(44) = xsubb
 a(49) = -2.*xsubb
 a(52) = -3.*xbsq
 a(55) = xcsq
 a(56) = xcyc
 a(57) = ycsq
 a(58) = xcsq * xsubc
 a(59) = ycsq * xsubc
 a(60) = ycsq * ysubc
 a(62) = xsubc
 a(63) = ysubc*2.
 a(65) =  xcyc*2.
 a(66) =  ycsq*3.
 a(67) =  -2.*xsubc
 a(68) =-ysubc
 a(70) =  -3.*xcsq
 a(71) =-ycsq
 
!     ******************************************************************
 
 IF (t2 == 0.) GO TO 220
 
!     ALL OF OPERATIONS THRU STMT 220
!     ARE NECESSARY IF T2 IS NON-ZERO.
 
!     ******************************************************************
 
 
!     GET THE G2X2 MATRIX
 
 matid = matid2
 inflag = 3
 CALL mat( ecpt(1) )
 IF(g2x211 == 0.0E0 .AND. g2x212 == 0.0E0 .AND. g2x222 == 0.0E0) GO TO 220
 g2x2(1) = g2x211 * t2
 g2x2(2) = g2x212 * t2
 g2x2(3) = g2x212 * t2
 g2x2(4) = g2x222 * t2
 
 determ = g2x2(1) * g2x2(4)  -  g2x2(3) * g2x2(2)
 j2x2(1) = g2x2(4) / determ
 j2x2(2) =-g2x2(2) / determ
 j2x2(3) =-g2x2(3) / determ
 j2x2(4) = g2x2(1) / determ
 
!     ******************************************************************
 
!     (H  ) IS PARTITIONED INTO A LEFT AND RIGHT PORTION AND ONLY THE
!       YQ  RIGHT PORTION IS COMPUTED AND USED AS A  (2X3). THE LEFT
!           2X3 PORTION IS NULL.  THE RIGHT PORTION WILL BE STORED AT
!           A(73)...A(78) UNTIL NOT NEEDED ANY FURTHER.
 
 
 
 temp = 2.*d(2) + 4.*d(9)
 a(73) =  -6.* (j2x2(1)*d(1) + j2x2(2)*d(3))
 a(74) =  -j2x2(1)*temp - 6.*j2x2(2)*d(6)
 a(75) =  -6.*(j2x2(1)*d(6)  + j2x2(2)*d(5))
 a(76) =  -6.*(j2x2(2)*d(1) + j2x2(4)*d(3))
 a(77) =  -j2x2(2)*temp - 6.*j2x2(4)*d(6)
 a(78) =  -6.*(j2x2(2)*d(6) + j2x2(4)*d(5))
 
!     THE ABOVE 6 ELEMENTS NOW REPRESENT THE (H  ) MATRIX (2X3)
!                                              YQ
 
!     NOW FORMING  PRODUCT (G2X2)(H  ) AND STORING AS AN INTERMEDIATE
!     STEP.                        YQ
 
 
 CALL gmmats(g2x2(1),2,2,0, a(73),2,3,0,  a(79))
 
!                                                               Y
!     WITH LAST PRODUCT  FORM  LOWER RIGHT 3 X 3 PARTITION OF (K )
 
!              Y                   T
!     THUS   (K ) PARTITION = (H  ) (LAST PRODUCT)   STORE AT A(85)
!                               YQ
 
 CALL gmmats(a(73),2,3,1,  a(79),2,3,0,  a(85))
 
!                                                     X
!     NOW ADD THE 9 ELEMENTS OF THIS 3X3 PORTION TO (K )
!     PER STEP 5 PAGE -16- MS-17                            Y
!     MULTIPLY IN AREA AT SAME TIME WHICH WAS LEFT OUT OF (K ) ABOVE.
 
 DO  i=1,3
   a(i + 21) = a(i + 21) + a(i + 84) * area
   a(i + 27) = a(i + 27) + a(i + 87) * area
   a(i + 33) = a(i + 33) + a(i + 90)*area
 END DO
 
!     ADD TO 6 OF THE (HBAR) ELEMENTS THE RESULT OF (H  )(H  )
!                                                    UY   YQ
!     THE PRODUCT IS FORMED DIRECTLY IN THE ADDITION PROCESS BELOW.
!     NO (H  ) MATRIX IS ACTUALLY COMPUTED DIRECTLY.
!          UY
 
!     THE FOLLOWING IS THEN PER STEPS 6 AND 7 PAGE -16- MS-17.
 
 DO  i=1,3
   a(i + 39) = a(i + 39) + xsubb * a(i + 72)
   a(i+57) =  a(i+57) + xsubc*a(i+72) + ysubc*a(i+75)
 END DO
 
!     THIS ENDS ADDED COMPUTATION FOR CASE OF T2 NOT ZERO
 
!     ******************************************************************
 
 220 CONTINUE
 
!     AT THIS POINT INVERT  (H) WHICH IS STORED AT A(37). . .A(72)
!     STORE INVERSE BACK IN A(37) . . . A(72)
!     NO NEED TO COMPUTE THE INVERSE SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL invers(6,a(37),6,a(73),0,determ,ising,a(79))
 
!     CHECK TO SEE IF H WAS SINGULAR
 IF (ising == 2) GO TO 7790
 
 
!     ******************************************************************
!              Q   -1
!     FORM   (K )(H  )  AND STORE AT  A(73). . .A(108)
!                 X                     Q
!     NOTE THAT (K ) AT THIS POINT IS (K )...
 
 CALL gmmats(a(1),6,6,0, a(37),6,6,0, a(73))
 
!                    -1 T
!     FORM(K  ) = (H  ) (LAST PRODUCT) STORE AT A(109). . .A(144)
!            II
 
 CALL gmmats(a(37),6,6,1, a(73),6,6,0, a(109))
 
!     ******************************************************************
 
!     FILL S-MATRIX EQUIVALENCED TO A(82)  (S IS  6X3 )
 
 IF (iopt /= 2) GO TO 260
 
!     SAVE H-INVERSE TO BE USED BY TRIANGULAR PLATE ROUTINE.
 
 DO  i = 37,72
   a(i+108) = a(i)
 END DO
 
 260 s(1) = 1.0
 s(2) = 0.
 s( 3) =-xsubb
 s(4) = 0.
 s(5) = 1.
 s(6) = 0.
 s(7) = 0.
 s(8) = 0.
 s(9) = 1.
 s(10)= 1.
 s(11) = ysubc
 s(12) =-xsubc
 s(13)= 0.
 s(14)= 1.
 s(15)= 0.
 s(16)= 0.
 s(17)= 0.
 s(18)= 1.
 
!     ******************************************************************
!                     T
!     FORM   K   = K   = -K   S  STORING AT A(46)   (K   IS 6X3)
!             IA    AI     II                         IA
 
 CALL gmmats(a(109),6,6,0,  s(1),6,3,0,   a(46))
 
!     THIS PRODUCT IS MULTIPLIED BY SCALER -1 BELOW.
 
!                    T
!          (K  ) = (S )(-K  )
!            AA           IA
 
!     NOTE K    HAS NOT BEEN MULTIPLIED ABOVE BY -1, THUS IGNORE MINUS
!           IA                                                   HERE.
 
 CALL gmmats(s(1),6,3,1,    a(46),6,3,0,  a(1))
 
!     NOW MULTIPLY  K   BY SCALER (-1)
!                    IA
 
 DO  i =46,63
   a(i) = -a(i)
 END DO
 
!     AT THIS POINT, STORED BY ROWS ARE
 
!                  K     (6X6) AT A(109). . .A(144)
!                   II
 
!                  K     (6,3) AT  A(46). . .A(63)
!                   IA
 
!                  K     (3X3) AT A(  1). . .A(  9)
!                   AA
 
!     ARRANGE NINE 3X3 MATRICES OF K SUPER U
 DO  i = 28,36
   a(i) = a(i+18)
 END DO
 a(10) = a(46)
 a(11) = a(49)
 a(12) = a(52)
 a(13) = a(47)
 a(14) = a(50)
 a(15) = a(53)
 a(16) = a(48)
 a(17) = a(51)
 a(18) = a(54)
 a(19) = a(55)
 a(20) = a(58)
 a(21) = a(61)
 a(22) = a(56)
 a(23) = a(59)
 a(24) = a(62)
 a(25) = a(57)
 a(26) = a(60)
 a(27) = a(63)
 a(37) = a(109)
 a(38) = a(110)
 a(39) = a(111)
 a(40) = a(115)
 a(41) = a(116)
 a(42) = a(117)
 a(43) = a(121)
 a(44) = a(122)
 a(45) = a(123)
 a(46) = a(112)
 a(47) = a(113)
 a(48) = a(114)
 a(49) = a(118)
 a(50) = a(119)
 a(51) = a(120)
 a(52) = a(124)
 a(53) = a(125)
 a(54) = a(126)
 a(64) = a(127)
 a(65) = a(128)
 a(66) = a(129)
 a(67) = a(133)
 a(68) = a(134)
 a(69) = a(135)
 a(70) = a(139)
 a(71) = a(140)
 a(72) = a(141)
 a(73) = a(130)
 a(74) = a(131)
 a(75) = a(132)
 a(76) = a(136)
 a(77) = a(137)
 a(78) = a(138)
 a(79) = a(142)
 a(80) = a(143)
 a(81) = a(144)
 dict5 = gsube
 IF (ntype /= 1)  GO TO 350
 loopnd = no(iopt)
 DO  i=1,loopnd
   aout(i) = a(i)
 END DO
 RETURN
 350 CONTINUE
 
!     ******************************************************************
 
 
 DO  npivot = 1,3
   
   
!     AT THIS POINT START ASSEMBLY OF 3   6X6 MATRICES FOR I = PIVOT,
!     AND J =1,2,3  IN THE FOLLOWING EQUATION.
   
!                   T         U    T
!        (K  ) = (T  ) (E) (K  ) (E ) (T )
!          IJ      I         IJ         J
   
!     ******************************************************************
   
!     FIRST GET THE PRODUCT APPLICABLE TO ALL 3 K  .
!                                                IJ
!                   T
!              = (T  ) (E)    A 6X3 MATRIX.
!                  I
   
   
!     CHECK TO SEE IF TI-MATRIX IS NEEDED
!     IF THE CSID IS ZERO FOR THE PIVOT POINT SKIP TRANSFORMATION.
   
   IF (necpt(4*npivot+9) == 0) GO TO 420
   
!     GET  TI AND MULTIPLY WITH E TO FILL TITE (THE COMMON PRODUCT)
   
   CALL transs(necpt(4*npivot+9),ti)
   
!     TI IS EQUIVALENCED TO A(118) AND IS 3X3.
   
!     FORM TITE (UPPER AND LOWER) OK OK OK....
   
   CALL gmmats(ti(1),3,3,1, e(1),3,3,0,  tite(1))
   CALL gmmats(ti(1),3,3,1, e(10),3,3,0, tite(10))
   
   GO TO 460
   
! 250 COMING HERE IMPLIES TI NOT USED.
!     JUST SET TITE = E MATRIX
   420 DO  i=1,18
     tite(i) = e(i)
   END DO
   
!     ******************************************************************
!                                                  T
! 280 AT THIS POINT COMMON PRODUCT IS COMPLETE =(T  )(E)  STORED IN TITE
!                                                 I
   
!     THE PIVOT I IS NPIVOT
   460   npt1 = 189
   
!     THE ABOVE SETS A POINTER, NPT1, TO POINT TO 18 FREE DOUBLE PREC.
!     CORE LOCATIONS IN THE A-ARRAY FOR STORAGE OF THE FOLLOWING
!     SUB-PRODUCT.
!                     U   T
!                  (K  )(E )(T )
!                    IJ       J
   
!     ******************************************************************
   
!     LOOP THRU FOR THE 3 - 6X6 K   ARRAYS.
!                                IJ
   DO  j=1,3
!                          T
!     TAKE SUB PRODUCT = (E )(T )..     STORE IN TJTE MATRIX
!                              J
     
!     NOTE.. THE TRANSPOSE OF THE ABOVE IS BEING FOUND AND USED,
!                           T
!                      = (T  )(E),  AND STORED IN TJTE-MATRIX
!                          J        EQUIVALENCED TO A(100)
     
     
!     CHECK TO SEE IF TRANSFORMATION IS NEEDED.
!     IF NOT SKIP TO  480
     
     IF (necpt(4*j +9)  == 0)  GO TO 480
     
     CALL transs (necpt (4*j+9),ti)
     
     CALL gmmats (ti(1),3,3,1, e(1),3,3,0,  tjte(1))
     CALL gmmats (ti(1),3,3,1,  e(10),3,3,0, tjte(10))
     GO TO 500
     
! 480 COMING HERE IF TRANSFORMATION NOT USED
     
! 480 SET TJTE = E
     480 DO  i =1,18
       tjte(i) = e(i)
     END DO
     
!     ******************************************************************
!           T       T
! 880   ( (E )(T ) )  IS COMPLETE AND STORED BY ROWS IN TJTE-MATRIX.
!               J
!                     U   T
!     NOW FORM,    (K  )(E )(T ), STORING AT A(NPT1)
!                    IJ       J
     
!                                  U
!     TO COMPUTE ABOVE USE 3X3   K
!                                 (NPIVOT,J)
!     COMPUTE POINTER TO THIS 3X3.
     
     500 npt2 = 27*npivot+9*j  - 35
     
     CALL gmmats (a(npt2),3,3,0, tjte,6,3,1,   a(npt1))
     
!     ******************************************************************
     
! 950 AT THIS POINT,
!                       U   T
!                    (K  )(E )(T )  IS STORED AT A(NPT1), (3X6).
!                      IJ       J
     
!     AND,              T
!                    (T  )(E)      IS STORED AT TITE(1) = A(82)  (6X3)
!                      I
!     ******************************************************************
     
!     FORMING FINAL PRODUCT, AND STORING AT A(100) THE 6X6.
     
     CALL gmmats (tite(1),6,3,0, a(npt1),3,6,0,  a(100))
     
!     ******************************************************************
     
     
!    NOW STORE  THE 6X6 MATRIX IN AOUT
     
     iout = (npivot-1)*27 + (j-1)*9 + 1
     i = 113
     DO  ii=1,3
       DO  jj=1,3
         ia = i + (ii-1)*6 + jj
         aout(iout) = a(ia)
         iout = iout+1
       END DO
     END DO
     
   END DO
   
 END DO
 
 RETURN
 
!     ERROR RETURNS
 
 7770 CALL mesage (30,31,ecpt(1))
 7777 nogo = .true.
 RETURN
 
 7780 CALL mesage (30,32,ielid)
 GO TO 7777
 
 7790 CALL mesage (30,33,ielid)
 GO TO 7777
!     ******************************************************************
END SUBROUTINE etrbks
