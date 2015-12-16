SUBROUTINE mtrbsc
     
!OMMENT.  ALL WRITE STATEMENTS WHICH HAVE BEEN COMMENTED OUT, HAVE BEEN
!         LEFT IN THE PROGRAMMING FOR ANY FUTURE DEBUGGING USE.
 
 
!      ************* BASIC BENDING TRIANGLE   ELEMENT ROUTINE **********
 
!     CALLS FROM THIS ROUTINE ARE MADE TO. . .
 
!          MAT    - MATERIAL DATA ROUTINE
!          SMA2B  - INSERTION ROUTINE
!          TRANSD - DOUBLE PRECISION TRANSFORMATION SUPPLIER
!          INVERD - DOUBLE PRECISION INVERSE ROUTINE
!          GMMATD - DOUBLE PRECISION MATRIX MULTIPLY AND TRANSPOSE
!          MESAGE - ERROR MESSAGE WRITER
 
 
!     ******************************************************************
 
 DOUBLE PRECISION :: a        ,e        ,xsubb    ,temp  &
     ,xsubc    ,d        ,ysubc    ,xcyc ,xcsq     ,determ   ,ycsq     ,xbsq  &
     ,g2x2     ,j2x2     ,hyq      ,aij ,bij      ,siij     ,sizero   ,mbaraa  &
     ,mar      ,mrr      ,s ,prod9    ,temp9    ,g  &
     ,yprodj   ,xprodi ,fj       ,fj2      ,fi       ,fij
 
 
 DIMENSION d(9)      ,g(9)      ,g2x2(4)   ,j2x2(4)   , s(18)  &
     ,ecpt(1)   ,hyq(6)    ,siij(7,7) ,mbaraa(9) , mar(18) ,mrr(36)
!     DIMENSION MNAME(9)
!     DIMENSION NASTER(130)
!     DATA (MNAME(I), I = 1,9) /6H1(MAA),6H (MAB),6H (MAC),6H (MBA),
!    $6H (MBB),6H (MBC),6H (MCA),6H (MCB),6H (MCC) /
!     DATA NASTER /130*1H*/
 
 COMMON /sma2io/ dum1(10), ifmgg, dum2(25)
 COMMON /sma2cl/  dum3(2), npvt ,                  dumcl(7)  &
     ,                  link(10)           ,nogo
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     t_sub_0, g_sub_e, sigten, sigcom, sigshe, g2x211, g2x212, g2x222, SPACE(2)
 
!     ECPT BLOCK
 COMMON /sma2et/ necpt(1)      ,ngrid(3)  &
     ,angle         ,matid1 ,eye           ,matid2  &
     ,t2            ,fmu ,z11           ,z22  &
     ,dummy1        ,x1 ,y1            ,z1  &
     ,dummy2        ,x2 ,y2            ,z2  &
     ,dummy3        ,x3 ,y3            ,z3            ,dumb(76)
 
 COMMON /sma2dp/    a(225)        ,prod9(9) ,temp9(9)      ,xsubb  &
     ,xsubc         ,ysubc ,e(9)          ,temp  &
     ,xcsq          ,xbsq ,ycsq          ,xcyc  &
     ,aij           ,determ ,bij           ,sizero  &
     ,fj            ,fj2 ,fi            ,fij  &
     ,yprodj        ,xprodi   &
     ,ising         ,dummy(59)
 
 EQUIVALENCE (d(1),g(1),siij(1,1),a(1))  ,(ecpt(1),necpt(1))  &
     ,(g2x2(1),a(10))             ,(j2x2(1),a(14))  &
     ,(hyq(1),a(50))              ,(mbaraa(1),a(136))  &
     ,(mar(1),a(145))             ,(mrr(1),a(163)) ,(s(1),a(82))
 
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
 
!     SETTING UP G MATRIX
 
 inflag = 2
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
 
!     COMPUTATION OF D = I.G-MATRIX (EYE IS INPUT FROM THE ECPT)
 
 DO  i = 1,9
   d(i) = g(i) * DBLE(eye)
 END DO
 
!     F1LL  (HBAR) MATRIX STORING AT A(100). . .A(135)
 xcsq = xsubc ** 2
 ycsq = ysubc ** 2
 xbsq = xsubb ** 2
 xcyc = xsubc * ysubc
 
 DO  i = 100,135
   a(i) = 0.0D0
 END DO
 
 a(100) = xbsq
 a(103) = xbsq * xsubb
 a(107) = xsubb
 a(112) = -2.0D0 * xsubb
 a(115) = -3.0D0 * xbsq
 a(118) = xcsq
 a(119) = xcyc
 a(120) = ycsq
 a(121) = xcsq * xsubc
 a(122) = ycsq * xsubc
 a(123) = ycsq * ysubc
 a(125) = xsubc
 a(126) = ysubc * 2.0D0
 a(128) = xcyc  * 2.0D0
 a(129) = ycsq  * 3.0D0
 a(130) =-2.0D0 * xsubc
 a(131) =-ysubc
 a(133) =-3.0D0 * xcsq
 a(134) =-ycsq
 
 
!     ******************************************************************
 
 IF( t2 == 0.0E0 ) GO TO 110
 
!     ALL OF THE FOLLOWING OPERATIONS THROUGH STATEMENT LABEL 110
!     ARE NECESSARY IF T2 IS NON-ZERO.
 
 
!     GET THE G2X2 MATRIX
 
 matid = matid2
 inflag = 3
 CALL mat( ecpt(1) )
 IF(g2x211 == 0.0E0 .AND. g2x212 == 0.0E0 .AND. g2x222 == 0.0E0) GO TO 110
 g2x2(1) = g2x211 * t2
 g2x2(2) = g2x212 * t2
 g2x2(3) = g2x2(2)
 g2x2(4) = g2x222 * t2
 
 determ = g2x2(1) * g2x2(4)  -  g2x2(3) * g2x2(2)
 j2x2(1) = g2x2(4) / determ
 j2x2(2) =-g2x2(2) / determ
 j2x2(3) = j2x2(2)
 j2x2(4) = g2x2(1) / determ
 
 
!     (H  ) IS PARTITIONED INTO A LEFT AND RIGHT PORTION AND ONLY THE
!       YQ  RIGHT PORTION IS COMPUTED AND USED AS A  (2X3). THE LEFT
!           2X3 PORTION IS NULL.  THE RIGHT PORTION WILL BE STORED AT
!           A(50)...A(55) UNTIL NOT NEEDED ANY FURTHER.
 
 
 
 temp  = 2.0D0 * d(2) + 4.0D0 * d(9)
 hyq(1) = -6.0D0 * (j2x2(1) * d(1) + j2x2(2) * d(3))
 hyq(2) = -j2x2(1) * temp - 6.0D0 * j2x2(2) * d(6)
 hyq(3) = -6.0D0 * (j2x2(1) * d(6) + j2x2(2) * d(5))
 hyq(4) = -6.0D0 * (j2x2(2) * d(1) + j2x2(4) * d(3))
 hyq(5) = -j2x2(2) * temp - 6.0D0 * j2x2(4) * d(6)
 hyq(6) = -6.0D0 * (j2x2(2) * d(6) + j2x2(4) * d(5))
 
!     ADD TO 6 OF THE (HBAR) ELEMENTS THE RESULT OF (H  )(H  )
!                                                    UY   YQ
!     THE PRODUCT IS FORMED DIRECTLY IN THE ADDITION PROCESS BELOW.
!     NO (H  ) MATRIX IS ACTUALLY COMPUTED DIRECTLY.
!          UY
 
!     THE FOLLOWING IS THEN STEP 6 PAGE 8, FMMS-66
 
 DO  i = 1,3
   a(i + 102) = a(i + 102) + xsubb * hyq(i)
   a(i + 120) = a(i + 120) + xsubc * hyq(i)    + ysubc * hyq(i + 3)
 END DO
 
!     THIS ENDS ADDED COMPUTATION FOR CASE OF T2 NOT ZERO
 
!     ******************************************************************
 
 110 CONTINUE
 
!     AT THIS POINT INVERT  (H) WHICH IS STORED AT A(100). . .A(135)
!     STORE INVERSE BACK IN A(100). . A(135)
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL inverd(6,a(100),6,a(136),0,determ,ising,a(142))
 
!     CHECK TO SEE IF H WAS SINGULAR
 IF( ising /= 2 ) GO TO 120
 
 
!     ISING = 2 IMPLIES SINGULAR MATRIX THUS ERROR CONDITION.
 CALL mesage(30,33,ecpt(1))
 
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
 nogo=1
 RETURN
 
 
!HURN OUT INTEGRAL VALUES I   USED IN REFERENCED M MATRICES
!                          IJ                           SEE P.9, FMMS-66
 
!     THE CALCULATION FOR  (I  ) ARE AS FOLLOWS
!                            IJ
!                                                      ***
!         A1  = XSUBB * YSUBC**(J+1) / ((J+1)*(J+2))      *
!           0J                                            *
!                                                         *
!         B   = XSUBC * YSUBC**(J+1) / (J+2)              *
!           0J                                            ** J=0,6
!                                                         *
!         A   = A1   + B                                  *
!           0J    0J    0J                                *
!                                                         *
!         I   = MU * A1                                   *
!           0J         0J                              ***
 
!                                                            ***
!         A1  = I * XSUBB * A      /(I+J+2)                     *
!           IJ               I-1,J                              *
!                                                               *
!         B   = XSUBC**(I+1) * YSUBC**(J+1) /((I+1)*(I+J+2))    *  I=1,6
!           IJ                                                  ** J=0,6
!                                                               *
!         A   = A1   + B                                        *
!           IJ    IJ    IJ                                      *
!                                                               *
!         I     MU * A1                                         *
!           IJ=        IJ                                       *
!                                                            ***
!      NOTE.. LOOPS FOR PROGRAM BEGIN AT 1 INSTEAD OF 0
!                                      I.E.  I = 1,7
!                                            J = 1,7
 
 120 DO  j=1,7
   yprodj = ysubc **j
   fj  = j
   fj2 = j+1
   aij       = xsubb * yprodj /(fj * fj2)
   bij       = xsubc * yprodj / fj2
   siij(1,j) = fmu   * aij
   aij       = aij   + bij
   IF(j == 7) CYCLE
   k = 8 - j
   DO  i = 2,k
     xprodi = xsubc **i
     fi     = i
     fij    = i + j
     aij       = (fi-1.0D0) * xsubb * aij / fij
     bij       = xprodi * yprodj /(fi * fij)
     siij(i,j) = fmu * aij
     aij       = aij + bij
   END DO
   
 END DO
 sizero = siij(1,1) / 3.0D0
 
!HUNK IN NUMBERS FOR (M-BAR-AA)    3X3 MATRIX AS PER MS-48, PAGES 6-10
 
!                    (M  )         3X6 MATRIX
!                      AR
 
!                    (M  )         6X6 MATRIX
!                      RR
 
!     (M-BAR-AA) MATRIX
 
 mbaraa(1) =  siij(1,1)
 mbaraa(2) =  siij(1,2)
 mbaraa(3) = -siij(2,1)
 mbaraa(4) =  siij(1,2)
 mbaraa(5) =  siij(1,3)
 mbaraa(6) = -siij(2,2)
 mbaraa(7) = -siij(2,1)
 mbaraa(8) = -siij(2,2)
 mbaraa(9) =  siij(3,1)
 
!     (M  ) MATRIX
!       AR
 mar( 1) = siij(3,1)
 mar( 2) = siij(2,2)
 mar( 3) = siij(1,3)
 mar( 4) = siij(4,1)
 mar( 5) = siij(2,3)
 mar( 6) = siij(1,4)
 mar( 7) = siij(3,2)
 mar( 8) = siij(2,3)
 mar( 9) = siij(1,4)
 mar(10) = siij(4,2)
 mar(11) = siij(2,4)
 mar(12) = siij(1,5)
 mar(13) =-siij(4,1)
 mar(14) =-siij(3,2)
 mar(15) =-siij(2,3)
 mar(16) =-siij(5,1)
 mar(17) =-siij(3,3)
 mar(18) =-siij(2,4)
 
!     (M  ) MATRIX  A 6X6 SYMMETRIC MATRIX
!       RR
 mrr( 1) = siij(5,1)
 mrr( 2) = siij(4,2)
 mrr( 3) = siij(3,3)
 mrr( 4) = siij(6,1)
 mrr( 5) = siij(4,3)
 mrr( 6) = siij(3,4)
 mrr( 7) = mrr(2)
 mrr( 8) = siij(3,3)
 mrr( 9) = siij(2,4)
 mrr(10) = siij(5,2)
 mrr(11) = siij(3,4)
 mrr(12) = siij(2,5)
 mrr(13) = mrr(3)
 mrr(14) = mrr(9)
 mrr(15) = siij(1,5)
 mrr(16) = siij(4,3)
 mrr(17) = siij(2,5)
 mrr(18) = siij(1,6)
 mrr(19) = mrr( 4)
 mrr(20) = mrr(10)
 mrr(21) = mrr(16)
 mrr(22) = siij(7,1)
 mrr(23) = siij(5,3)
 mrr(24) = siij(4,4)
 mrr(25) = mrr( 5)
 mrr(26) = mrr(11)
 mrr(27) = mrr(17)
 mrr(28) = mrr(23)
 mrr(29) = siij(3,5)
 mrr(30) = siij(2,6)
 mrr(31) = mrr( 6)
 mrr(32) = mrr(12)
 mrr(33) = mrr(18)
 mrr(34) = mrr(24)
 mrr(35) = mrr(30)
 mrr(36) = siij(1,7)
 
 IF(t2 == 0.0) GO TO 146
 IF(g2x211 == 0.0E0 .AND. g2x212 == 0.0E0 .AND. g2x222 == 0.0E0) GO TO 146
 
 mar( 4) = mar( 4)   + hyq(1) * siij(2,1) + hyq(4) * siij(1,2)
 mar( 5) = mar( 5)   + hyq(2) * siij(2,1) + hyq(5) * siij(1,2)
 mar( 6) = mar( 6)   + hyq(3) * siij(2,1) + hyq(6) * siij(1,2)
 mar(10) = mar(10)   + hyq(1) * siij(2,2) + hyq(4) * siij(1,3)
 mar(11) = mar(11)   + hyq(2) * siij(2,2) + hyq(5) * siij(1,3)
 mar(12) = mar(12)   + hyq(3) * siij(2,2) + hyq(6) * siij(1,3)
 mar(16) = mar(16)   - hyq(1) * siij(3,1) - hyq(4) * siij(2,2)
 mar(17) = mar(17)   - hyq(2) * siij(3,1) - hyq(5) * siij(2,2)
 mar(18) = mar(18)   - hyq(3) * siij(3,1) - hyq(6) * siij(2,2)
 mrr( 4) = mrr( 4)   + hyq(1) * siij(4,1) + hyq(4) * siij(3,2)
 mrr( 5) = mrr( 5)   + hyq(2) * siij(4,1) + hyq(5) * siij(3,2)
 mrr( 6) = mrr( 6)   + hyq(3) * siij(4,1) + hyq(6) * siij(3,2)
 mrr(10) = mrr(10)   + hyq(1) * siij(3,2) + hyq(4) * siij(2,3)
 mrr(11) = mrr(11)   + hyq(2) * siij(3,2) + hyq(5) * siij(2,3)
 mrr(12) = mrr(12)   + hyq(3) * siij(3,2) + hyq(6) * siij(2,3)
 mrr(16) = mrr(16)   + hyq(1) * siij(2,3) + hyq(4) * siij(1,4)
 mrr(17) = mrr(17)   + hyq(2) * siij(2,3) + hyq(5) * siij(1,4)
 mrr(18) = mrr(18)   + hyq(3) * siij(2,3) + hyq(6) * siij(1,4)
 mrr(19) = mrr( 4)
 mrr(20) = mrr(10)
 mrr(21) = mrr(16)
 mrr(22) = mrr(22)   + hyq(1) * (hyq(1) * siij(3,1) + 2.0D0  &
     * (siij(5,1) + hyq(4) * siij(2,2))) + hyq(4) * (2.0D0 * siij(4,2)  &
     + hyq(4) * siij(1,3))
 mrr(23) = mrr(23)   + hyq(2) * siij(5,1) + hyq(5) * siij(4,2)  &
     + hyq(1) * (siij(3,3) + hyq(2) * siij(3,1) + hyq(5) * siij(2,2))  &
     + hyq(4) * (siij(2,4) + hyq(2) * siij(2,2) + hyq(5) * siij(1,3))
 mrr(24) = mrr(24)   + hyq(3) * siij(5,1) + hyq(6) * siij(4,2)  &
     + hyq(1) * (siij(2,4) + hyq(3) * siij(3,1) + hyq(6) * siij(2,2))  &
     + hyq(4) * (siij(1,5) + hyq(3) * siij(2,2) + hyq(6) * siij(1,3))
 mrr(25) = mrr( 5)
 mrr(26) = mrr(11)
 mrr(27) = mrr(17)
 mrr(28) = mrr(23)
 mrr(29) = mrr(29)   + hyq(2) * (hyq(2) * siij(3,1) + 2.0D0  &
     * (siij(3,3) + hyq(5) * siij(2,2))) + hyq(5) * (2.0D0 * siij(2,4)  &
     + hyq(5) * siij(1,3))
 mrr(30) = mrr(30)   + hyq(3) * siij(3,3) + hyq(6) * siij(2,4)  &
     + hyq(2) * (siij(2,4) + hyq(3) * siij(3,1) + hyq(6) * siij(2,2))  &
     + hyq(5) * (siij(1,5) + hyq(3) * siij(2,2) + hyq(6) * siij(1,3))
 mrr(31) = mrr( 6)
 mrr(32) = mrr(12)
 mrr(33) = mrr(18)
 mrr(34) = mrr(24)
 mrr(35) = mrr(30)
 mrr(36) = mrr(36)   + hyq(3) * (hyq(3) * siij(3,1) + 2.0D0  &
     * (siij(2,4) + hyq(6) * siij(2,2))) + hyq(6) * (2.0D0 * siij(1,5)  &
     + hyq(6) * siij(1,3))
 
 146 CONTINUE
 
 
!     FILL S-MATRIX EQUIVALENCED TO A(82)  (S IS  6X3 )
 
 s( 1) = 1.0D0
 s( 2) = 0.0D0
 s( 3) =-xsubb
 s( 4) = 0.0D0
 s( 5) = 1.0D0
 s( 6) = 0.0D0
 s( 7) = 0.0D0
 s( 8) = 0.0D0
 s( 9) = 1.0D0
 s(10) = 1.0D0
 s(11) = ysubc
 s(12) =-xsubc
 s(13) = 0.0D0
 s(14) = 1.0D0
 s(15) = 0.0D0
 s(16) = 0.0D0
 s(17) = 0.0D0
 s(18) = 1.0D0
 
!AN NOW COMPUTE 9 (3X3) MASS MATRICES (FMMS-66, PAGES 10-11)
 
 
!                -1 T           -1
!     ( M ) = ( H  )  ( M  ) ( H  )
!                        RR
 
!              PARTITION (M)
!                                           ///       ///
!                                           /     *     /
!                                           / MBB * MBC /
!                                           /     *     /
!                                ( M )  =   / ********* /
!                                           /     *     /
!                                           / MCB * MCC /
!                                           /     *     /
!                                           ///       ///
!                                                       4 (3X3) MATRICES
!                        -1
!     ( M  ) = ( M  ) ( H  )
!        AI       AR
 
!              PARTITION (M  )              ///                 ///
!                          AI               /          *          /
!                               ( M  )  =   / M-BAR-AB * M-BAR-AC /
!                                  AI       /          *          /
!                                           ///                 ///
!                                                       2 (3X3) MATRICES
!                               T            T
!     ( MAB )  = (M-BAR-AB) - (S ) (MBB) - (S ) (MCB)
!                               B            C
 
!                               T            T
!     ( MAC )  = (M-BAR-AC) - (S ) (MBC) - (S ) (MCC)
!                               B            C
 
!                               T     T      T      T
!     ( MAA )  = (M-BAR-AA) - (S ) (M  ) - (S ) (MAC )
!                               B    AB      C
 
!                           - (M-BAR-AB) (S ) - (M-BAR-AC) (S )
!                                          B                 C
 
!                    T
!     ( MBA )  = (MAB )
 
!                    T
!     ( MCA )  = (MAC )
 
!HOOSE APPROPRIATE BLOCK OF A-ARRAY FOR STORAGE
 
!     (3X3)    STORED IN      (3X3)     STORED IN     (3X3)    STORED IN
!     (MAA)   A( 1... 9)      (MAB)   A(10)...8)      (MAC)   A(19...27)
!     (MBA)   A(28...36)      (MBB)   A(37)...45)     (MBC)   A(46...54)
!     (MCA)   A(55...63)      (MCB)   A(64...72)      (MCC)   A(73...81)
 
!       -1
!     (H  ) IS STORED AT A(100...135)
!     (S)   EQUIVALENCED A( 81... 99)
!     WORKING STORAGE IS A(181...216)
!     (M-BAR-AB) STORED UNTIL NO LONGER NEEDED IN A(163...171)
!     (M-BAR-AC) STORED UNTIL NO LONGER NEEDED IN A(172...180)
 
!               -1 T          -1
!OMPUTE (M) = (H  )  ((M  ) (H  ))
!                       RR
 
 CALL gmmatd(mrr(1),6,6,0,a(100),6,6,0,a(37))
 CALL gmmatd(a(100),6,6,1,a( 37),6,6,0,a( 1))
 
!REATE PARTITION OF 4 (3X3)
 DO  i=1,3
   a(i+36) = a(i   )
   a(i+39) = a(i+ 6)
   a(i+42) = a(i+12)
   
   a(i+45) = a(i+ 3)
   a(i+48) = a(i+ 9)
   a(i+51) = a(i+15)
   
   a(i+63) = a(i+18)
   a(i+66) = a(i+24)
   a(i+69) = a(i+30)
   
   a(i+72) = a(i+21)
   a(i+75) = a(i+27)
   a(i+78) = a(i+33)
 END DO
 
!OMPUTE                 -1
!       (M  ) = (M  ) (H  )    AND  PARTITION INTO 2 (3X3)  (M-BAR-AB)
!         AI      AR                                    AND (M-BAR-AC)
 
 CALL gmmatd(mar(1),3,6,0,a(100),6,6,0,a(181))
 DO  i=1,3
   a(i+162) = a(i+180)
   a(i+165) = a(i+186)
   a(i+168) = a(i+192)
   
   a(i+171) = a(i+183)
   a(i+174) = a(i+189)
   a(i+177) = a(i+195)
 END DO
!OMPUTE (MAB)
 CALL gmmatd(s( 1),3,3,1,a(37),3,3,0,a(181))
 CALL gmmatd(s(10),3,3,1,a(64),3,3,0,a(190))
 DO  i=1,9
   a(i+9)=a(i+162) - a(i+180) - a(i+189)
 END DO
!OMPUTE (MAC)
 CALL gmmatd(s( 1),3,3,1,a(46),3,3,0,a(181))
 CALL gmmatd(s(10),3,3,1,a(73),3,3,0,a(190))
 DO  i=1,9
   a(i+18) = a(i+171) - a(i+180) - a(i+189)
 END DO
!OMPUTE (MAA)
 CALL gmmatd(s(  1),3,3,1,a(10),3,3,1,a(181))
 CALL gmmatd(s( 10),3,3,1,a(19),3,3,1,a(190))
 CALL gmmatd(a(163),3,3,0,s( 1),3,3,0,a(199))
 CALL gmmatd(a(172),3,3,0,s(10),3,3,0,a(208))
 DO  i=1,9
   a(i) = mbaraa(i) - a(i+180) - a(i+189) - a(i+198) - a(i+207)
 END DO
!OMPUTE (MBA) AND (MCA)
 DO  i=1,3
   npt = 3 * i + 7
   a(i+27) = a(npt)
   a(i+30) = a(npt + 1)
   a(i+33) = a(npt + 2)
   
   a(i+54) = a(npt +  9)
   a(i+57) = a(npt + 10)
   a(i+60) = a(npt + 11)
 END DO
 
 RETURN
 
END SUBROUTINE mtrbsc
