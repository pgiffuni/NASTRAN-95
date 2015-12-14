SUBROUTINE pstrb1 (iopt)
     
!     THIS ROUTINE DOES SUB-CALCULATIONS FOR PLATE ELEMENTS IN PLA3
 
!     THIS ROUTINE IS SIMILAR TO STRBS1, BUT SINCE THE BASIC BENDING
!     TRIANGLE (IOPT = 0) IS NOT USED IN PLA, THE CORRESPONDING
!     EXECUTIABLE CODE FOR THAT CASE IS NOT USED.
 
!     PHASE ONE FOR STRESS RECOVERY
 
!              IOPT   = 0  (BASIC BENDING TRIANGLE)
!              IOPT   = 1  (SUB-CALCULATIONS FOR SQDPL1)
!              IOPT   = 2  (SUB-CALCULATIONS FOR STRPL1)
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
 
!     PLAMAT - ROTATES AND RETURNS GP MATRIX
!              MAT    - MATERIAL DATA ROUTINE
!              TRANSS - SINGLE PRECISION TRANSFORMATION SUPPLIER
!              INVERS - SINGLE PRECISION INVERSE ROUTINE
!              GMMATS - SINGLE PRECISION MATRIX MULTIPLY AND TRANSPOSE
!              MESAGE - ERROR MESSAGE WRITER
 
 
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER :: subsca,subscb
 REAL :: ks,j2x2
 DIMENSION       d(9),g2x2(4),j2x2(4),s(18),ecpt(1),g(9),hic(18),  &
     hib(18),tite(18),t(9),ks(30),hinv(36)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     t_sub_0, g sub e, sigten, sigcom, sigshe, g2x211, g2x212, g2x222
 COMMON /pla32s/ a(225),xsubb,xsubc,ysubc,e(18),temp,xbar,area,  &
     xcsq,ybar2,ycsq,ybar,xbsq,px2,xcyc,py2,pxy2,xbar3,  &
     ybar3,determ,prod9(9),temp9(9),nsized,dumdum(4),  &
     npivot,theta ,nsubc,ising,subsca,subscb,nerror,  &
     nbegin,ntyped,xc,yc,yc2,yc3,isub,xc3,dum55(26)
 COMMON /pla3es/ necpt(1),ngrid(3),angle,matid1,eye,matid2,t2,fmu,  &
     z11,z22,dummy1,x1,y1,z1,dummy2,x2,y2,z2,dummy3,x3, y3,z3,dumb(76),ph1out(200)
 EQUIVALENCE     (d(1),g(1),a(79)),(ecpt(1),necpt(1)),  &
     (ks(1),ph1out(1)),(g2x2(1),a(88)),(s(1),a(55)),  &
     (tite(1),a(127)),(j2x2(1),a(92)),(t(1),a(118)),  &
     (hib(1),a(109)),(hic(1),a(127)),(hinv(1),a(73))
 
!     ECPT LIST FOR BASIC BENDING TRIANGLE           NAME IN
!                                                    THIS
!     ECPT                                           ROUTINE   TYPE
!     ==========================================     ========  =======
!     ECPT( 1) = ELEMENT ID                          NECPT(1)  INTEGER
!     ECPT( 2) = GRID POINT A                        NGRID(1)  INTEGER
!     ECPT( 3) = GRID POINT B                        NGRID(2)  INTEGER
!     ECPT( 4) = GRID POINT C                        NGRID(3)  INTEGER
!     ECPT( 5) = THETA = ANGLE OF MATERIAL           ANGLE     REAL
!     ECPT( 6) = MATERIAL ID 1                       MATID1    INTEGER
!     ECPT( 7) = I = MOMENT OF INERTIA               EYE       REAL
!     ECPT( 8) = MATERIAL ID 2                       MATID2    INTEGER
!     ECPT( 9) = T2                                  T2        REAL
!     ECPT(10) = NON-STRUCTURAL-MASS                 FMU       REAL
!     ECPT(11) = Z1                                  Z11       REAL
!     ECPT(12) = Z2                                  Z22       REAL
!     ECPT(13) = COORD. SYSTEM ID 1                  NECPT(13) INTEGER
!     ECPT(14) = X1                                  X1        REAL
!     ECPT(15) = Y1                                  Y1        REAL
!     ECPT(16) = Z1                                  Z1        REAL
!     ECPT(17) = COORD. SYSTEM ID 2                  NECPT(17) INTEGER
!     ECPT(18) = X2                                  X2        REAL
!     ECPT(19) = Y2                                  Y2        REAL
!     ECPT(20) = Z2                                  Z2        REAL
!     ECPT(21) = COORD. SYSTEM ID 3                  NECPT(21) INTEGER
!     ECPT(22) = X3                                  X3        REAL
!     ECPT(23) = Y3                                  Y3        REAL
!     ECPT(24) = Z3                                  Z3        REAL
!     ECPT(25) = ELEMENT TEMPERATURE                 ELTEMP    REAL
 
 matid  = matid1
 inflag = -1
 
 CALL plamat
 
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
   d(i) = g(i)*eye
 END DO
 
 xbar = (xsubb+xsubc)/3.0
 ybar = ysubc/3.0
 
!     FORMING K  5X6 AND STORING TEMPORARILY IN PH1OUT OUTPUT SPACE.
!              S                             (EQUIVALENCED)
 
 xc3    = 3.0*xc
 yc3    = 3.0*yc
 yc2    = 2.0*yc
 ks( 1) = d(1)
 ks( 2) = d(3)
 ks( 3) = d(2)
 ks( 4) = d(1)*xc3
 ks( 5) = d(2)*xc + d(3)*yc2
 ks( 6) = d(2)*yc3
 ks( 7) = d(2)
 ks( 8) = d(6)
 ks( 9) = d(5)
 ks(10) = d(2)*xc3
 ks(11) = d(5)*xc + d(6)*yc2
 ks(12) = d(5)*yc3
 ks(13) = d(3)
 ks(14) = d(9)
 ks(15) = d(6)
 ks(16) = d(3)*xc3
 ks(17) = d(6)*xc + d(9)*yc2
 ks(18) = d(6)*yc3
 
!     ROWS 4 AND 5
 
 ks(19) = 0.0
 ks(20) = 0.0
 ks(21) = 0.0
 ks(22) =-d(1)*6.0
 ks(23) =-d(2)*2.0 - d(9)*4.0
 ks(24) =-d(6)*6.0
 ks(25) = 0.0
 ks(26) = 0.0
 ks(27) = 0.0
 ks(28) =-d(3)*6.0
 ks(29) =-d(6)*6.0
 ks(30) =-d(5)*6.0
 
!     MULTIPLY FIRST 3 ROWS BY 2.0
 
 DO  i = 1,18
   ks(i) = ks(i)*2.0
 END DO
 
 xcsq = xsubc**2
 ycsq = ysubc**2
 xbsq = xsubb**2
 xcyc = xsubc*ysubc
 
!     F1LL (HBAR) MATRIX STORING AT A(37) THRU A(72)
 
 DO  i = 37,72
   a(i)  = 0.0
 END DO
 a(37) = xbsq
 a(40) = xbsq*xsubb
 a(44) = xsubb
 a(49) =-2.0*xsubb
 a(52) =-3.0*xbsq
 a(55) = xcsq
 a(56) = xcyc
 a(57) = ycsq
 a(58) = xcsq*xsubc
 a(59) = ycsq*xsubc
 a(60) = ycsq*ysubc
 a(62) = xsubc
 a(63) = ysubc*2.0
 a(65) = xcyc *2.0
 a(66) = ycsq *3.0
 a(67) =-2.0*xsubc
 a(68) =-ysubc
 a(70) =-3.0*xcsq
 a(71) =-ycsq
 
 IF (t2 == 0.0) GO TO 110
 
!     ALL OF THE FOLLOWING OPERATIONS THROUGH STATEMENT LABEL 100
!     ARE NECESSARY IF T2 IS NON-ZERO.
 
!     GET THE G2X2 MATRIX
 
 matid  = matid2
 inflag = 3
 CALL mat (ecpt(1))
 IF (g2x211 == 0.0 .AND. g2x212 == 0.0 .AND. g2x222 == 0.0) GO TO 110
 g2x2(1) = g2x211*t2
 g2x2(2) = g2x212*t2
 g2x2(3) = g2x212*t2
 g2x2(4) = g2x222*t2
 
 determ  = g2x2(1)*g2x2(4) - g2x2(3)*g2x2(2)
 j2x2(1) = g2x2(4)/determ
 j2x2(2) =-g2x2(2)/determ
 j2x2(3) =-g2x2(3)/determ
 j2x2(4) = g2x2(1)/determ
 
!     (H  ) IS PARTITIONED INTO A LEFT AND RIGHT PORTION AND ONLY THE
!       YQ  RIGHT PORTION IS COMPUTED AND USED AS A  (2X3). THE LEFT
!           2X3 PORTION IS NULL.  THE RIGHT PORTION WILL BE STORED AT
!           A(73) THRU A(78) UNTIL NOT NEEDED ANY FURTHER.
 
 temp  =  2.0*d(2) + 4.0*d(9)
 a(73) = -6.0*(j2x2(1)*d(1) + j2x2(2)*d(3))
 a(74) = -j2x2(1)*temp + 6.0*j2x2(2)*d(6)
 a(75) = -6.0*(j2x2(1)*d(6) + j2x2(2)*d(5))
 a(76) = -6.0*(j2x2(2)*d(1) + j2x2(4)*d(3))
 a(77) = -j2x2(2)*temp + 6.0*j2x2(4)*d(6)
 a(78) = -6.0*(j2x2(2)*d(6) + j2x2(4)*d(5))
 
!     THE ABOVE 6 ELEMENTS NOW REPRESENT THE (H  ) MATRIX (2X3)
!                                              YQ
 
!     ADD TO 6 OF THE (HBAR) ELEMENTS THE RESULT OF(H  )(H  )
!                                                    UY   YQ
!     THE PRODUCT IS FORMED DIRECTLY IN THE ADDITION PROCESS BELOW.
!     NO (H  ) MATRIX IS ACTUALLY COMPUTED DIRECTLY.
!          UY
 
!     THE FOLLOWING IS THEN PER STEPS 6 AND 7 PAGE -16- MS-17.
 
 DO  i = 1,3
   a(i+39) = a(i+39) + xsubb*a(i+72)
   a(i+57) = a(i+57) + xsubc*a(i+72) + ysubc*a(i+75)
 END DO
 
!     THIS ENDS ADDED COMPUTATION FOR CASE OF T2 NOT ZERO
 
 110 CONTINUE
 
!     AT THIS POINT INVERT  (H) WHICH IS STORED AT A(37) THRU A(72)
!     STORE INVERSE BACK IN A(37) THRU A(72)
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUNTLY.
 
 ising = -1
 CALL invers (6,a(37),6,a(73),0,determ,ising,a(79))
 
!     CHECK TO SEE IF H WAS SINGULAR
 
 IF (ising /= 2) GO TO 120
 
!     ISING = 2 IMPLIES SINGULAR MATRIX THUS ERROR CONDITION.
 
 CALL mesage (-30,38,ecpt(1))
 
!     SAVE H-INVERSE IF TRI-PLATE IS CALLING
 
 120 DO  i = 1,36
   hinv(i) = a(i+36)
 END DO
 
!     FILL  S-MATRIX, EQUIVALENCED TO A(55).  (6X3)
 
 s( 1) = 1.0
 s( 2) = 0.0
 s( 3) =-xsubb
 s( 4) = 0.0
 s( 5) = 1.0
 s( 6) = 0.0
 s( 7) = 0.0
 s( 8) = 0.0
 s( 9) = 1.0
 s(10) = 1.0
 s(11) = ysubc
 s(12) =-xsubc
 s(13) = 0.0
 s(14) = 1.0
 s(15) = 0.0
 s(16) = 0.0
 s(17) = 0.0
 s(18) = 1.0
 
!     COMPUTE  S , S ,  AND S    NO TRANSFORMATIONS
!               A   B        C
 
!                -1
!     S  = - K  H  S ,   S  = K  H   ,   S  = K  H
!      A      S           B    S  IB      C    S  IC
 
!     S   COMPUTATION.
!      A
 
 CALL gmmats (hinv(1),6,6,0, s(1),6,3,0, a(16))
 
!     DIVIDE  H-INVERSE INTO A LEFT 6X3 AND RIGHT 6X3 PARTITION.
 
 i = 0
 j =-6
 150 j = j + 6
 k = 0
 160 k = k + 1
 i = i + 1
 isub = j + k
 hib(i) = hinv(isub    )
 hic(i) = hinv(isub + 3)
 IF (k <  3) GO TO 160
 IF (j < 30) GO TO 150
 
 CALL gmmats (ks(1),5,6,0, a(16),6,3,0, a(1))
 
!     MULTIPLY S SUB A BY -1
 
 DO  i = 1,15
   a(i) = -a(i)
 END DO
 
!     S  COMPUTATION
!      B
 
 CALL gmmats (ks,5,6,0, hib,6,3,0, a(16))
 
!     S  COMPUTATION
!      C
 
 CALL gmmats (ks,5,6,0, hic,6,3,0, a(31))
 
 RETURN
END SUBROUTINE pstrb1