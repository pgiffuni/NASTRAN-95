SUBROUTINE strsl1
     
!     OUTPUTS FROM THIS PHASE FOR USE IN PHASE II ARE THE FOLLOWING
 
!     1)  ELEMENT ID              WORDS    1    STORAGE IN PH1OUT  1
!     2)  SIX SILS                WORDS    6                     2-7
!     3)  3 MEMBRANE  THICKNESES  WORDS    3                     8-10
!     4)  3 BENDING  THICKNESES   WORDS    3                    11-13
!     5)  8  STRESS DATA POINTS   WORDS    8                    14-21
!     6)  4 NOS. STRESS MATRICES (6-5X6 EACH) WORDS 720         22-741
!     7)  S SUB T MATRIX          WORDS  4X3                   742-753
!     8)  ELEMENT ID              WORD     1                   754
!     9)  SIX SILS                WORDS    6                   755-760
!     10) ELEMENT TEMPERATURE     WORD     1                   761
!     11) 4 NOS. MEMBRANE STRESS MATRICES 4(6-3X3)             762-1193
 
!     ECPT LISTS
 
!     ECPT ( 1) = ELEMENT ID                                    INTEGER
!     ECPT ( 2) = SCALAR INDEX NUMBER FOR GRID POINT 1          INTEGER
!     ECPT ( 3) = SCALAR INDEX NUMBER FOR GRID POINT 2          INTEGER
!     ECPT ( 4) = SCALAR INDEX NUMBER FOR GRID POINT 3          INTEGER
!     ECPT ( 5) = SCALAR INDEX NUMBER FOR GRID POINT 4          INTEGER
!     ECPT ( 6) = SCALAR INDEX NUMBER FOR GRID POINT 5          INTEGER
!     ECPT ( 7) = SCALAR INDEX NUMBER FOR GRID POINT 6          INTEGER
!     ECPT ( 8) = THETA                                         REAL
!     ECPT ( 9) = MATERIAL ID 1                                 INTEGER
!     ECPT (10) = THICKNESS T1 AT GRID POINT G1
!     ECPT (11) = THICKNESS T3 AT GRID POINT G3
!     ECPT (12) = THICKNESS T5 AT GRID POINT G5
!     ECPT (13) = MATERIAL ID 2                                 INTEGER
!     ECPT (14) = THICKNESS TBEND1 FOR BENDING AT GRID POINT G1
!     ECPT (15) = THICKNESS TBEND3 FOR BENDING AT GRID POINT G3
!     ECPT (16) = THICKNESS TBEND5 FOR BENDING AT GRID POINT G5
!     ECPT (17) = MATERIAL ID 3                                 INTEGER
!     ECPT (18) = THICKNESS TSHR1 FOR TRANSVERSE SHEAR AT GRID POINT G1
!     ECPT (19) = THICKNESS TSHR3 FOR TRANSVERSE SHEAR AT GRID POINT G3
!     ECPT (20) = THICKNESS TSHR5 FOR TRANSVERSE SHEAR AT GRID POINT G5
!     ECPT (21) = NON-STRUCTURAL MASS                           REAL
!     ECPT (22) = DISTANCE Z11 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT (23) = DISTANCE Z21 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT (24) = DISTANCE Z13 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT (25) = DISTANCE Z23 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT (26) = DISTANCE 015 FOR STRESS CALCULATION  AT GRID POINT G5
!     ECPT (27) = DISTANCE Z25 FOR STRESS CALCULATION  AT GRID POINT G5
 
!     X1,Y1,Z1 FOR ALL SIX POINTS ARE IN NASTRAN BASIC SYSTEM
 
!     ECPT (28) = CO-ORDINATE SYSTEM ID FOR GRID A              INTEGER
!     ECPT (29) = CO-ORDINATE X1                                REAL
!     ECPT (30) = CO-ORDINATE Y1                                REAL
!     ECPT (31) = CO-ORDINATE Z1                                REAL
!     ECPT (32) = CO-ORDINATE SYSTEM ID FOR GRID B              INTEGER
!     ECPT (33) = CO-ORDINATE X1                                REAL
!     ECPT (34) = CO-ORDINATE Y1                                REAL
!     ECPT (35) = CO-ORDINATE Z1                                REAL
!     ECPT (36) = CO-ORDINATE SYSTEM ID FOR GRID C              INTEGE9
!     ECPT (37) = CO-ORDINATE X1                                REAL
!     ECPT (38) = CO-ORDINATE Y1                                REAL
!     ECPT (39) = CO-ORDINATE Z1                                REAL
!     ECPT (40) = CO-ORDINATE SYSTEM ID FOR GRID D              INTEGER
!     ECPT (41) = CO-ORDINATE X1                                REAL
!     ECPT (42) = CO-ORDINATE Y1                                REAL
!     ECPT (43) = CO-ORDINATE Z1                                REAL
!     ECPT (44) = CO-ORDINATE SYSTEM ID FOR GRID E              INTEGER
!     ECPT (45) = CO-ORDINATE X1                                REAL
!     ECPT (46) = CO-ORDINATE Y1                                REAL
!     ECPT (47) = CO-ORDINATE Z1                                REAL
!     ECPT (48) = CO-ORDINATE SYSTEM ID FOR GRID F              INTEGER
!     ECPT (49) = CO-ORDINATE X1                                REAL
!     ECPT (50) = CO-ORDINATE Y1                                REAL
!     ECPT (51) = CO-ORDINATE Z1                                REAL
!     ECPT (52) = ELEMENT TEMPERATURE  AT GRID POINTS G1        REAL
!     ECPT (53) = ELEMENT TEMPERATURE  AT GRID POINTS G2        REAL
!     ECPT (54) = ELEMENT TEMPERATURE  AT GRID POINTS G3        REAL
!     ECPT (55) = ELEMENT TEMPERATURE  AT GRID POINTS G4        REAL
!     ECPT (56) = ELEMENT TEMPERATURE  AT GRID POINTS G5        REAL
!     ECPT (57) = ELEMENT TEMPERATURE  AT GRID POINTS G6        REAL
 
 LOGICAL :: nots
 REAL :: j11,j12,j22,nsm,ivect(3),jvect(3),kvect(3)
 DIMENSION       NAME(2),INDEX(20,3),ics(6),nl(6),xc(6),yc(6),  &
     zc(6),qqq(20,20),qqqinv(360),ts6(40),ts7(60),  &
     e(18),v1(3),v2(3),v3(3),iest(42),e1(18),ph1ben(9),  &
     ph1shr(6),ph2(18),ph3(12),ph4(90),tmmm(36),  &
     q(6,6),ind(6,3),cab(3),ee( 30),ph1mem(6),eph1(15),  &
     si(9),emod(9),d(9),dph1(9),g(4),gph1(6),  &
     nph1ou(990),tm(96),tmqq(90),ee1(5,6),tmm(3,12),  &
     tmb(60),tmbq(54),trans(9),balotr(36)
 COMMON /matin / matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/ em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy,  &
     rj11,rj12,rj22
 COMMON /sdr2x5/ est(100),ph1out(1200),forvec(24),x,y,z,dista,  &
     distb,distc,a1,a2,a3,aa1,aa2,aa3,qqqinv,qq,tm,  &
     tmqq,ts6,ts7,q,ee,ee1,ph2,ph3,ph4,e,e1,xc,yc,zc,  &
     ph1mem,ph1ben,ph1shr,dph1,eph1,gph1,g,d,ics,nl, cab,trans,balotr,emod,si
 EQUIVALENCE     (a,dista),(b,distb),(c,distc),(v1(1),est(29)),  &
     (v2(1),est(37)),(v3(1),est(45)),(iest(1),est(1)),  &
     (tmmm(1),tmm(1,1)),(ph1out(1),qqq(1,1)),  &
     (ph1out(401),INDEX(1,1),ind(1,1)), (nph1ou(1),ph1out(1))
 DATA   degra  / 0.0174532925 /
 DATA   BLANK  , NAME / 4H    , 4HTRSH, 4HL     /
 
 nots   = .true.
 idele  = iest(1)
 DO  i = 1,6
   nl(i)  = iest(i+1)
 END DO
 thetam = est(8)
 matid1 = iest(9)
 tmem1  = est(10)
 tmem3  = est(11)
 tmem5  = est(12)
 matid2 = iest(13)
 tbend1 = (est(14)*12.0)**0.3333333333
 tbend3 = (est(15)*12.0)**0.3333333333
 tbend5 = (est(16)*12.0)**0.3333333333
 matid3 = iest(17)
 tshr1  = est(18)
 tshr3  = est(19)
 tshr5  = est(20)
 nsm    = est(21)
 j      = 0
 IF (tmem3 == 0.0 .OR. tmem3 == BLANK) tmem3 = tmem1
 IF (tmem5 == 0.0 .OR. tmem5 == BLANK) tmem5 = tmem1
 IF (tshr3 == 0.0 .OR. tshr3 == BLANK) tshr3 = tshr1
 IF (tshr5 == 0.0 .OR. tshr5 == BLANK) tshr5 = tshr1
 IF (tshr1 == 0.0) nots =.true.
 IF (tbend3 == 0.0 .OR. tbend3 == BLANK) tbend3 = tbend1
 IF (tbend5 == 0.0 .OR. tbend5 == BLANK) tbend5 = tbend1
 DO  i = 28,48,4
   j = j + 1
   ics(j) = iest(i )
   xc(j)  = est(i+1)
   yc(j)  = est(i+2)
   zc(j)  = est(i+3)
 END DO
 eltemp = est(52)
 theta1 = thetam*degra
 sinth  = SIN(theta1)
 costh  = COS(theta1)
 IF (ABS(sinth) <= 1.0E-06) sinth = 0.0
 
!     EVALUATE MATERIAL PROPERTIES
 
 matflg = 2
 matid  = matid1
 IF (matid1 == 0) GO TO 30
 CALL mat (idele)
 g11 = em(1)
 g12 = em(2)
 g13 = em(3)
 g22 = em(4)
 g23 = em(5)
 g33 = em(6)
 
 30 matid = matid2
 IF (matid2 == 0) GO TO 40
 matflg = 2
 CALL mat (idele)
 d11  = em(1)
 d12  = em(2)
 d13  = em(3)
 d21  = d12
 d22  = em(4)
 d23  = em(5)
 d31  = d13
 d32  = d23
 d33  = em(6)
 d(1) = d11
 d(2) = d12
 d(3) = d13
 d(4) = d21
 d(5) = d22
 d(6) = d23
 d(7) = d13
 d(8) = d23
 d(9) = d33
 d334 = d33*4.0
 d132 = d13*2.0
 d232 = d23*2.0
 j11  = 0.0
 j12  = 0.0
 j22  = 0.0
 IF (nots) GO TO 40
 CALL mat (idele)
 
!     CALCULATIONS FOR THE TRIANGLE
 
 40 CALL trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,iest(1),NAME)
 
!     EVALUATE THE CONSTANTS C1,C2,AND C3 IN THE LINEAR EQUATION FOR
!     THICKNESS VARIATION - MEMBRANE
 
 CALL af (f,1,a,b,c,c1,c2,c3,tmem1,tmem2,tmem3,1)
 cab(1) = c1
 cab(2) = c2
 cab(3) = c3
 
!     A1,A2,A3 ARE THE COEFFICIENTS OF LINEAR EQUATION FOR VARIATION
!     OF BENDING THICKNESSES
 
 CALL af (f,1,a,b,c,a1,a2,a3,tbend1,tbend3,tbend5,1)
 a1sq = a1*a1
 a2sq = a2*a2
 a3sq = a3*a3
 c1   = a1sq*a1
 c2   = 3.0*a1sq*a2
 c3   = 3.0*a1sq*a3
 c4   = 3.0*a1*a2sq
 c5   = 6.0*a1*a2*a3
 c6   = 3.0*a3sq*a1
 c7   = a2sq*a2
 c8   = 3.0*a2sq*a3
 c9   = 3.0*a2*a3sq
 c10  = a3*a3sq
 
!     AA1, AA2, AA3  ARE COEFFICIENTS IN THICKNESS VARIATION FOR
!     TRANSVERSE SHEAR
 
 CALL af (f,1,a,b,c,aa1,aa2,aa3,tshr1,tshr3,tshr5,1)
 
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e(i)  = 0.0
 END DO
 e( 1) = kvect(1)
 e( 4) = kvect(2)
 e( 7) = kvect(3)
 e(11) = ivect(1)
 e(14) = ivect(2)
 e(17) = ivect(3)
 e(12) = jvect(1)
 e(15) = jvect(2)
 e(18) = jvect(3)
 
!     CALCULATIONS FOR QMATRIX (QQQ) AND ITS INVERSE
 
 DO  i = 1,400
   qqq(i,1) = 0.0
 END DO
 DO  i = 1,6
   i1 = (i-1)*3 + 1
   i2 = (i-1)*3 + 2
   i3 = (i-1)*3 + 3
   qqq(i1, 1) = 1.0
   qqq(i1, 2) = xc(i)
   qqq(i1, 3) = yc(i)
   qqq(i1, 4) = xc(i)*xc(i)
   qqq(i1, 5) = xc(i)*yc(i)
   qqq(i1, 6) = yc(i)*yc(i)
   qqq(i1, 7) = qqq(i1, 4)*xc(i)
   qqq(i1, 8) = qqq(i1, 4)*yc(i)
   qqq(i1, 9) = qqq(i1, 5)*yc(i)
   qqq(i1,10) = qqq(i1, 6)*yc(i)
   qqq(i1,11) = qqq(i1, 7)*xc(i)
   qqq(i1,12) = qqq(i1, 7)*yc(i)
   qqq(i1,13) = qqq(i1, 8)*yc(i)
   qqq(i1,14) = qqq(i1, 9)*yc(i)
   qqq(i1,15) = qqq(i1,10)*yc(i)
   qqq(i1,16) = qqq(i1,11)*xc(i)
   qqq(i1,17) = qqq(i1,12)*yc(i)
   qqq(i1,18) = qqq(i1,13)*yc(i)
   qqq(i1,19) = qqq(i1,14)*yc(i)
   qqq(i1,20) = qqq(i1,15)*yc(i)
   qqq(i2, 3) = 1.0
   qqq(i2, 5) = xc(i)
   qqq(i2, 6) = yc(i)*2.0
   qqq(i2, 8) = qqq(i1, 4)
   qqq(i2, 9) = qqq(i1, 5)*2.0
   qqq(i2,10) = qqq(i1, 6)*3.0
   qqq(i2,12) = qqq(i1, 7)
   qqq(i2,13) = qqq(i1, 8)*2.0
   qqq(i2,14) = qqq(i1, 9)*3.0
   qqq(i2,15) = qqq(i1,10)*4.0
   qqq(i2,17) = qqq(i1,12)*2.0
   qqq(i2,18) = qqq(i1,13)*3.0
   qqq(i2,19) = qqq(i1,14)*4.0
   qqq(i2,20) = qqq(i1,15)*5.0
   qqq(i3, 2) =-1.0
   qqq(i3, 4) =-2.0*xc(i)
   qqq(i3, 5) =-yc(i)
   qqq(i3, 7) =-qqq(i1, 4)*3.0
   qqq(i3, 8) =-qqq(i1, 5)*2.0
   qqq(i3, 9) =-qqq(i1, 6)
   qqq(i3,11) =-qqq(i1, 7)*4.0
   qqq(i3,12) =-qqq(i1, 8)*3.0
   qqq(i3,13) =-qqq(i1, 9)*2.0
   qqq(i3,14) =-qqq(i1,10)
   qqq(i3,16) =-qqq(i1,11)*5.0
   qqq(i3,17) =-qqq(i1,13)*3.0
   qqq(i3,18) =-qqq(i1,14)*2.0
   qqq(i3,19) =-qqq(i1,15)
 END DO
 
 qqq(19,16) = 5.0*a**4*c
 qqq(19,17) = 3.0*a**2*c**3 - 2.0*a**4*c
 qqq(19,18) =-2.0*a*c**4 + 3.0*a**3*c**2
 qqq(19,19) = c**5 - 4.0*a**2*c**3
 qqq(19,20) = 5.0*a*c**4
 qqq(20,16) = 5.0*b**4*c
 qqq(20,17) = 3.0*b**2*c**3 - 2.0*b**4*c
 qqq(20,18) = 2.0*b*c**4 - 3.0*b**3*c**2
 qqq(20,19) = c**5 - 4.0*b**2*c**3
 qqq(20,20) =-5.0*b*c**4
 DO  i = 1,6
   DO  j = 1,6
     i1 = 3*(i-1) + 1
     q(i,j) = qqq(i1,j)
   END DO
 END DO
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (6,q,6,ts6(1),0,det,ising,ind)
 
!     FOURTH ARGUMENT IS A DUMMY LOCATION FOR INVERSE AND HENCE TS1(1)
!     IS U
 
!     SET ISING = -1 AGAIN.
 ising = -1
 CALL invers (20,qqq,20,ts6(1),0,determ,ising,INDEX)
 
!     ISING EQUAL TO 2 IMPLIES THAT QQQ IS SINGULAR
 
 
!     FIRST 18 COLUMNS OF QQQ INVERSE IS THE QQQINV FOR USE IN STIFFNESS
!     MATRIX CALCULATIONS
 
 h4 = q(4,1)*zc(1) + q(4,2)*zc(2) + q(4,3)*zc(3) + q(4,4)*zc(4) +  &
     q(4,5)*zc(5) + q(4,6)*zc(6)
 h5 = q(5,1)*zc(1) + q(5,2)*zc(2) + q(5,3)*zc(3) + q(5,4)*zc(4) +  &
     q(5,5)*zc(5) + q(5,6)*zc(6)
 h6 = q(6,1)*zc(1) + q(6,2)*zc(2) + q(6,3)*zc(3) + q(6,4)*zc(4) +  &
     q(6,5)*zc(5) + q(6,6)*zc(6)
 h4 = h4*2.0
 h6 = h6*2.0
 
!     H5 IS MULTIPLIED BY 2.0, SO THAT EXY=DU/DY + DV/DX - ZXY * W
 
 h5 = h5*2.0
 DO  i = 1,20
   DO  j = 1,18
     ij = (i-1)*18 + j
     qqqinv(ij) = qqq(i,j)
   END DO
 END DO
 DO  i = 1,36
   balotr(i) = 0.0
 END DO
 
 DO  i = 1,7
   nph1ou (i) = iest(i)
 END DO
 IF (matid2 == 0) GO TO 120
 GO TO 140
 120 DO  i = 2,7
   nph1ou (i) = 0
 END DO
 140 ph1out( 8) = tmem1
 ph1out( 9) = tmem3
 ph1out(10) = tmem5
 ph1out(11) = tbend1
 ph1out(12) = tbend3
 ph1out(13) = tbend5
 ph1out(14) = est(22)
 ph1out(15) = est(23)
 ph1out(16) = est(24)
 ph1out(17) = est(25)
 ph1out(18) = est(26)
 ph1out(19) = est(27)
 emod(1) = g11
 emod(2) = g12
 emod(3) = g13
 emod(4) = g12
 emod(5) = g22
 emod(6) = g23
 emod(7) = g13
 emod(8) = g23
 emod(9) = g33
 DO  i = 1,30
   ee(i)   = 0.0
 END DO
 ee( 1)  = ivect(1)
 ee( 2)  = jvect(1)
 ee( 3)  = kvect(1)
 ee( 6)  = ivect(2)
 ee( 7)  = jvect(2)
 ee( 8)  = kvect(2)
 ee(11)  = ivect(3)
 ee(12)  = jvect(3)
 ee(13)  = kvect(3)
 ee(19)  = ivect(1)
 ee(20)  = jvect(1)
 ee(24)  = ivect(2)
 ee(25)  = jvect(2)
 ee(29)  = ivect(3)
 ee(30)  = jvect(3)
 DO  jj = 1,4
   j = 2*jj - 1
   IF (jj == 4) GO TO 160
   x = xc(j)
   y = yc(j)
   GO TO 170
   160 x = (xc(1) + xc(3) + xc(5))/3.0
   y = (yc(1) + yc(3) + yc(5))/3.0
   ph1out(20) = (a1 + a2*x + a3*y)/2.0
   ph1out(21) =-ph1out(20)
   170 IF (matid2 == 0) GO TO 190
   DO  i = 1,60
     ts7(i) = 0.0
   END DO
   thk  = a1 + a2*x + a3*y
   thk1 = thk**3/12.0
   d(1) = d11*thk1
   d(2) = d12*thk1
   d(3) = d13*thk1
   d(4) = d(2)
   d(5) = d22*thk1
   d(6) = d23*thk1
   d(7) = d(3)
   d(8) = d(6)
   d(9) = d33*thk1
   x2   = x*x
   xy   = x*y
   y2   = y*y
   x3   = x2*x
   x2y  = x2*y
   xy2  = x*y2
   y3   = y2*y
   ts7( 4) = 2.0
   ts7( 7) = 6.0*x
   ts7( 8) = 2.0*y
   ts7(11) = 12.0*x2
   ts7(12) = 6.0*xy
   ts7(13) = 2.0*y2
   ts7(16) = 20.0*x3
   ts7(17) = 6.0*xy2
   ts7(18) = 2.0*y3
   ts7(26) = 2.0
   ts7(29) = 2.0*x
   ts7(30) = 6.0*y
   ts7(33) = 2.0*x2
   ts7(34) = ts7(12)
   ts7(35) = 12.0*y2
   ts7(37) = 2.0*x3
   ts7(38) = 6.0*x2y
   ts7(39) = 12.0*xy2
   ts7(40) = 20.0*y3
   ts7(45) = 2.0
   ts7(48) = 4.0*x
   ts7(49) = 4.0*y
   ts7(52) = 6.0*x2
   ts7(53) = 8.0*xy
   ts7(54) = 6.0*y2
   ts7(57) = 12.0*x2y
   ts7(58) = ts7(39)
   ts7(59) = 8.0*y3
   CALL gmmats (ts7,3,20,0,qqqinv,20,18,0,ph4(1))
   CALL strslv (ts6,nots)
   CALL gmmats (ts6,2,20,0,qqqinv,20,18,0,ph4(55))
   
   190 IF (matid1 == 0) GO TO 220
   DO  i = 1,36
     tmmm(i) = 0.0
   END DO
   DO  j = 1,6
     j1 = (j-1)*2 + 1
     j2 = j1 + 1
     tmm(1,j1) = q(2,j) + 2.0*x*q(4,j) + y*q(5,j)
     tmm(2,j2) = q(3,j) + x*q(5,j) + 2.0*y*q(6,j)
     tmm(3,j1) = tmm(2,j2)
     tmm(3,j2) = tmm(1,j1)
   END DO
   x4   = x3*x
   x3y  = x3*y
   x2y2 = x2*y2
   xy3  = x*y3
   y4   = y*y3
   x5   = x4*x
   x3y2 = x3*y2
   x2y3 = x2*y3
   xy4  = x*y4
   y5   = y*y4
   tmb( 1) = -h4
   tmb( 2) = -h4*x
   tmb( 3) = -h4*y
   tmb( 4) = -h4*x2
   tmb( 5) = -h4*xy
   tmb( 6) = -h4*y2
   tmb( 7) = -h4*x3
   tmb( 8) = -h4*x2y
   tmb( 9) = -h4*xy2
   tmb(10) = -h4*y3
   tmb(11) = -h4*x4
   tmb(12) = -h4*x3y
   tmb(13) = -h4*x2y2
   tmb(14) = -h4*xy3
   tmb(15) = -h4*y4
   tmb(16) = -h4*x5
   tmb(17) = -h4*x3y2
   tmb(18) = -h4*x2y3
   tmb(19) = -h4*xy4
   tmb(20) = -h4*y5
   tmb(21) = -h6
   tmb(22) = -h6*x
   tmb(23) = -h6*y
   tmb(24) = -h6*x2
   tmb(25) = -h6*xy
   tmb(26) = -h6*y2
   tmb(27) = -h6*x3
   tmb(28) = -h6*x2y
   tmb(29) = -h6*xy2
   tmb(30) = -h6*y3
   tmb(31) = -h6*x4
   tmb(32) = -h6*x3y
   tmb(33) = -h6*x2y2
   tmb(34) = -h6*xy3
   tmb(35) = -h6*y4
   tmb(36) = -h6*x5
   tmb(37) = -h6*x3y2
   tmb(38) = -h6*x2y3
   tmb(39) = -h6*xy4
   tmb(40) = -h6*y5
   tmb(41) = -h5
   tmb(42) = -h5*x
   tmb(43) = -h5*y
   tmb(44) = -h5*x2
   tmb(45) = -h5*xy
   tmb(46) = -h5*y2
   tmb(47) = -h5*x3
   tmb(48) = -h5*x2y
   tmb(49) = -h5*xy2
   tmb(50) = -h5*y3
   tmb(51) = -h5*x4
   tmb(52) = -h5*x3y
   tmb(53) = -h5*x2y2
   tmb(54) = -h5*xy3
   tmb(55) = -h5*y4
   tmb(56) = -h5*x5
   tmb(57) = -h5*x3y2
   tmb(58) = -h5*x2y3
   tmb(59) = -h5*xy4
   tmb(60) = -h5*y5
   CALL gmmats (tmb,3,20,0, qqqinv,20,18,0, tmbq)
   
   220 DO  ii = 1,6
     IF (ics(ii) == 0) GO TO 240
     CALL transs (iest(4*ii+24),trans)
     DO  j = 1,3
       l = 6*(j-1) + 1
       m = 3*(j-1) + 1
       balotr(l   ) = trans(m  )
       balotr(l+1 ) = trans(m+1)
       balotr(l+2 ) = trans(m+2)
       balotr(l+21) = trans(m  )
       balotr(l+22) = trans(m+1)
       balotr(l+23) = trans(m+2)
     END DO
     CALL gmmats (e,6,3,+1, balotr,6,6,0, e1)
     GO TO 260
     240 DO  i = 1,3
       DO  j = 1,6
         i1 = (i-1)*6 + j
         j1 = (j-1)*3 + i
         e1(i1) = e(j1)
       END DO
     END DO
     260 IF (matid2 == 0) GO TO 300
     kz = (ii-1)*3 + 1
     ph1ben(1) = ph4(kz   )
     ph1ben(2) = ph4(kz+ 1)
     ph1ben(3) = ph4(kz+ 2)
     ph1ben(4) = ph4(kz+18)
     ph1ben(5) = ph4(kz+19)
     ph1ben(6) = ph4(kz+20)
     ph1ben(7) = ph4(kz+36)
     ph1ben(8) = ph4(kz+37)
     ph1ben(9) = ph4(kz+38)
     CALL gmmats (d,3,3,0, ph1ben,3,3,0, dph1)
     CALL gmmats (dph1,3,3,0, e1,3,6,0, ph2)
     mz = (ii-1)*3 + 55
     ph1shr(1) = ph4(mz   )
     ph1shr(2) = ph4(mz+ 1)
     ph1shr(3) = ph4(mz+ 2)
     ph1shr(4) = ph4(mz+18)
     ph1shr(5) = ph4(mz+19)
     ph1shr(6) = ph4(mz+20)
     IF (nots) GO TO 270
     thk  = aa1 + aa2*x + aa3*y
     g(1) = em(6)*thk
     g(2) = 0.0
     g(3) = 0.0
     g(4) = g(1)
     CALL gmmats (g,2,2,0, ph1shr,2,3,0, gph1)
     GO TO 280
     270 gph1(1) = ph1shr(1)
     gph1(2) = ph1shr(2)
     gph1(3) = ph1shr(3)
     gph1(4) = ph1shr(4)
     gph1(5) = ph1shr(5)
     gph1(6) = ph1shr(6)
     280 CALL gmmats (gph1,2,3,0, e1,3,6,0, ph3)
     DO  i = 1,3
       DO  j = 1,6
         i1 = (i-1)*6 + j
         i2 = i1 + 18
         j1 = (ii-1)*30 + (jj-1)*180 + i1 + 21
         j2 = j1 + 18
         ph1out(j1) = ph2(i1)
         IF (i /= 3) ph1out(j2) = ph3(i1)
       END DO
     END DO
     
     300 IF (matid1 == 0) CYCLE
     DO  i = 1,3
       DO  j = 1,2
         ji = (i-1)*5 + j
         ij = (j-1)*3 + i + (ii-1)*6
         tm(ji) = tmmm(ij)
       END DO
     END DO
     DO  i = 1,3
       DO  j = 1,3
         ji = (i-1)*5  + j + 2
         ij = (i-1)*18 + j + (ii-1)*3
         tm(ji) = tmbq(ij)
       END DO
     END DO
     IF (ics(ii) /= 0) CALL gmmats (ee,6,5,+1, balotr,6,6,0, ee1)
     ij1 = (jj-1)*108 + (ii-1)*18 + 762
     CALL  gmmats (emod,3,3,0, tm(1),3,5,0, eph1)
     IF (ics(ii) == 0) CALL gmmats (eph1,3,5,0,ee,6,5,+1,ph1out(ij1))
     IF (ics(ii) /= 0) CALL gmmats (eph1,3,5,0,ee1,5,6,0,ph1out(ij1))
   END DO
 END DO
 
 jst = 742 + (jj-1)*3
 IF (matid2 /= 0) CALL gmmats (d,3,3,0,alf(1),3,1,0,ph1out(jst))
 IF (matid1 /= 0) CALL gmmats (emod,3,3,0,alf(1),3,1,0, ph1out(1194))
 IF (matid1 == 0) GO TO 360
 DO  i = 1,7
   nph1ou(753+i) = iest(i)
 END DO
 GO TO 380
 360 DO  i = 1,7
   nph1ou(753+i) = 0
 END DO
 380 ph1out(761) = tref
 RETURN
END SUBROUTINE strsl1
