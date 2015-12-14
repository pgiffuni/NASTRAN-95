SUBROUTINE dtshls
     
!     ECPT ENTRIES
 
!     ECPT( 1) = ELEMENT ID                                    INTEGER
!     ECPT( 2) = SCALAR INDEX NUMBER FOR GRID POINT 1          INTEGER
!     ECPT( 3) = SCALAR INDEX NUMBER FOR GRID POINT 2          INTEGER
!     ECPT( 4) = SCALAR INDEX NUMBER FOR GRID POINT 3          INTEGER
!     ECPT( 5) = SCALAR INDEX NUMBER FOR GRID POINT 4          INTEGER
!     ECPT( 6) = SCALAR INDEX NUMBER FOR GRID POINT 5          INTEGER
!     ECPT( 7) = SCALAR INDEX NUMBER FOR GRID POINT 6          INTEGER
!     ECPT( 8) = THETA                                         REAL
!     ECPT( 9) = MATERIAL  ID 1                                INTEGER
!     ECPT(10) = THICKNESS T1 AT GRID POINT G1
!     ECPT(11) = THICKNESS T3 AT GRID POINT G3
!     ECPT(12) = THICKNESS T5 AT GRID POINT G5
!     ECPT(13) = MATERIAL  ID 2                                INTEGER
!     ECPT(14) = THICKNESS TBEND1 FOR BENDING AT GRID POINT G1
!     ECPT(15) = THICKNESS TBEND3 FOR BENDING AT GRID POINT G3
!     ECPT(16) = THICKNESS TBEND5 FOR BENDING AT GRID POINT G5
!     ECPT(17) = MATERIAL  ID 3                                INTEGER
!     ECPT(18) = THICKNESS TSHR1 FOR TRANSVERSE SHEAR AT GRID POINT G1
!     ECPT(19) = THICKNESS TSHR3 FOR TRANSVERSE SHEAR AT GRID POINT G3
!     ECPT(20) = THICKNESS TSHR5 FOR TRANSVERSE SHEAR AT GRID POINT G5
!     ECPT(21) = NON-STRUCTURAL MASS                           REAL
!     ECPT(22) = DISTANCE Z11 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT(23) = DISTANCE Z21 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT(24) = DISTANCE Z13 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT(25) = DISTANCE Z23 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT(26) = DISTANCE Z15 FOR STRESS CALCULATION  AT GRID POINT G5
!     ECPT(27) = DISTANCE Z25 FOR STRESS CALCULATION  AT GRID POINT G5
 
!     X1,Y1,Z1 FOR ALL SIX POINTS ARE  IN NASTRAN BASIC SYSTEM
 
!     ECPT(28) = COORDINATE SYSTEM ID FOR GRID A               INTEGER
!     ECPT(29) = COORDINATE X1                                 REAL
!     ECPT(30) = COORDINATE Y1                                 REAL
!     ECPT(31) = COORDINATE Z1                                 REAL
!     ECPT(32) = COORDINATE SYSTEM ID FOR GRID B               INTEGER
!     ECPT(33) = COORDINATE X1                                 REAL
!     ECPT(34) = COORDINATE Y1                                 REAL
!     ECPT(35) = COORDINATE Z1                                 REAL
!     ECPT(36) = COORDINATE SYSTEM ID FOR GRID C               INTEGER
!     ECPT(37) = COORDINATE X1                                 REAL
!     ECPT(38) = COORDINATE Y1                                 REAL
!     ECPT(39) = COORDINATE Z1                                 REAL
!     ECPT(40) = COORDINATE SYSTEM ID FOR GRID D               INTEGER
!     ECPT(41) = COORDINATE X1                                 REAL
!     ECPT(42) = COORDINATE Y1                                 REAL
!     ECPT(43) = COORDINATE Z1                                 REAL
!     ECPT(44) = COORDINATE SYSTEM ID FOR GRID E               INTEGER
!     ECPT(45) = COORDINATE X1                                 REAL
!     ECPT(46) = COORDINATE Y1                                 REAL
!     ECPT(47) = COORDINATE Z1                                 REAL
!     ECPT(48) = COORDINATE SYSTEM ID FOR GRID F               INTEGER
!     ECPT(49) = COORDINATE X1                                 REAL
!     ECPT(50) = COORDINATE Y1                                 REAL
!     ECPT(51) = COORDINATE Z1                                 REAL
!     EST (52) = ELEMENT  TEMPERATURE
!     EST (53) = ENFORCED ELEMENT DEFORMATION (NOT USED)
!     EST (54) = LOADING  TEMPERATURE
!     EST (55) TO EST (90) = GLOBAL DISPLACEMENT VECTOR
!                REPLACES ECPT(65) TO ECPT(100) DESCRIBED BELOW
!     ECPT(65) = U1-DISP FOR X1
!     ECPT(66) = V1-DISP FOR Y1
!     ECPT(67) = W1-DISP FOR Z1
!     ECPT(68) = ALFA1-ROTATION FOR X1
!     ECPT(69) = BETA1-ROTATION FOR Y1
!     ECPT(70) = GAMA1-ROTATION FOR Z1
!     ECPT(71) = U2-DISP FOR X2
!     ECPT(72) = V2-DISP FOR Y2
!     ECPT(73) = W2-DISP FOR Z2
!     ECPT(74) = ALFA2-ROTATION FOR X2
!     ECPT(75) = BETA2-ROTATION FOR Y2
!     ECPT(76) = GAMA2-ROTATION FOR Z2
!     ECPT(77) = U3-DISP FOR X3
!     ECPT(78) = V3-DISP FOR Y3
!     ECPT(79) = W3-DISP FOR Z3
!     ECPT(80) = ALFA3-ROTATION FOR X3
!     ECPT(81) = BETA3-ROTATION FOR Y3
!     ECPT(82) = GAMA3-ROTATION FOR Z3
!     ECPT(83) = U4-DISP FOR X4
!     ECPT(84) = V4-DISP FOR Y4
!     ECPT(85) = W4-DISP FOR Z4
!     ECPT(86) = ALFA4-ROTATION FOR X4
!     ECPT(87) = BETA4-ROTATION FOR Y4
!     ECPT(88) = GAMA4-ROTATION FOR Z4
!     ECPT(89) = U5-DISP FOR X5
!     ECPT(90) = V5-DISP FOR Y5
!     ECPT(91) = W5-DISP FOR Z5
!     ECPT(92) = ALFA5-ROTATION FOR X5
!     ECPT(93) = BETA5-ROTATION FOR Y5
!     ECPT(94) = GAMA5-ROTATION FOR Z5
!     ECPT(95) = U6-DISP FOR X6
!     ECPT(96) = V6-DISP FOR Y6
!     ECPT(97) = W6-DISP FOR Z6
!     ECPT(98) = ALFA6-ROTATION FOR X6
!     ECPT(99) = BETA6-ROTATION FOR Y6
!     ECPT(100)= GAMA6-ROTATION FOR Z6
 
!     RK AND SK ARE EXPONENTS IN THICKNESS VARIATION
 
 LOGICAL :: nots,unimem,uniben,nogo
 INTEGER :: rk(3),sk(3),rl(3),sl(3),xu(32),yu(32),xv(32),  &
     yv(32),xw(32),yw(32),sil(6),sil1,sil2, rr,rr0,rr1,ss,ss0,ss1
 REAL :: j11,j12,j22,nsm,ivect(3),jvect(3),kvect(3),xc(6),  &
     yc(6),zc(6),f(18,18)
 REAL :: trand(9),balotr(36),ksub(36),ksubt(36)
 REAL :: d334,d132,d232,rmx,rnx,rmnx,rmx1,rnx1,rmy,rny,  &
     rmny,rmy1,rny1,x,y,qqq(20,20),cmt(1296),  &
     ctm(36,36),cms(900),cm1(30,30),cab(3),csub(5,5),  &
     csubt(6,5),c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,  &
     h4,h5,h6,sb1,sb2,sb3,sb4,sb5,sb6,sb7,sb8,sb9,  &
     rix,riy,rjx,rjy,rkx,rky,rlx,rly,ee(30),q(6,6),  &
     qqqinv(360),qks(960),kshl(1024),mshl(1024)
 REAL :: sb10,sb11,sb12,sb13,sb14,sb15,sb16,sb17,sb18,sb19  &
     ,                sb20,sb21,sb22,sb23,sb24,sb25,sb26,sb27,sb28,sb29  &
     ,                sb30,sb31,sb32,sb33,sb34,sb35,sb36,sb37,sb38,sb39  &
     ,                sb40,cc(10),st
 DOUBLE PRECISION :: ksubd(36)
 DIMENSION        ind(6,3),el(3),fl(3),gl(3),NAME(2),INDEX(20,3),  &
     ics(6),iest(100),nl(6),sigx(3),sigy(3),sigxy(3),  &
     es(6),stress(3),str(3),vec(3),ph1out(250),  &
     tm(3,12),emod(9),tmmm(36),trans(9),eph1(6), ee1(6),nsil(6),ti(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /system/  ibuf,ioutpt
 COMMON /ds1aet/  est(100)
 COMMON /ds1aaa/  npvt,icstm,ncstm
 COMMON /ds1adp/  f
 COMMON /matin /  matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/  em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy,  &
     rj11,rj12,rj22
 EQUIVALENCE      (c1,cc(1)),(c2,cc(2)),(c3,cc(3)),(c4,cc(4)),  &
     (c5,cc(5)),(c6,cc(6)),(c7,cc(7)),(c8,cc(8)), (c9,cc(9)),(c10,cc(10)),  &
     (nsil(1),ph1out(2)),(tm(1,1),tmmm(1)),(a,dista),  &
     (b,distb),(c,distc),(iest(1),est(1)),  &
     (cm1(1,1),cms(1)),(thk1,tbend1),(thk2,tbend3),  &
     (thk3,tbend5),(cmt(1025),qqqinv(1)),  &
     (ctm(1,1),cmt(1),kshl(1),mshl(1),qqq(1,1)),  &
     (cmt(437),ph1out(1)),(cmt(687),INDEX(1,1)),  &
     (cmt(747),ind(1,1)),(ti(1),est(65))
 DATA    rk    /  0,1,0 /, rl / 0,1,0 /, sk / 0,0,1 /, sl / 0,0,1/,  &
     xu    /  0,1,0,2,1,0,26*0    /, yu / 0,0,1,0,1,2,26*0   /,  &
     xv    /  6*0,0,1,0,2,1,0,20*0/, yv /6*0,0,0,1,0,1,2,20*0/,  &
     xw    /  12*0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0   /,  &
     yw    /  12*0,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5   /,  &
     BLANK ,  NAME / 4H    ,4HDTSH,4HLD    /, degra /  0.0174532925 /
 
 nots   =.false.
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
 tbend1 = (est(14)*12.0)**0.333333333333
 tbend3 = (est(15)*12.0)**0.333333333333
 tbend5 = (est(16)*12.0)**0.333333333333
 matid3 = iest(17)
 tshr1  = est(18)
 tshr3  = est(19)
 tshr5  = est(20)
 nsm    = est(21)
 j      = 0
 DO  i = 28,48,4
   j      = j + 1
   ics(j) = iest(i)
   xc(j)  = est(i+1)
   yc(j)  = est(i+2)
   zc(j)  = est(i+3)
 END DO
 
!     IF TMEM3 OR TMEM5 EQUAL TO ZERO OR BLANK, THEY WILL BE
!     SET EQUAL TO TMEM1 SO ALSO FOR TSHR3,TSHR5,TBEND3 AND TBEND5
 
 IF (tmem3 == 0.0 .OR. tmem3 == BLANK) tmem3 = tmem1
 IF (tmem5 == 0.0 .OR. tmem5 == BLANK) tmem5 = tmem1
 IF (tshr3 == 0.0 .OR. tshr3 == BLANK) tshr3 = tshr1
 IF (tshr5 == 0.0 .OR. tshr5 == BLANK) tshr5 = tshr1
 tshr   = (tshr1+tshr3+tshr5)/3.0
 IF (tshr1 == 0.0) nots =.true.
 IF (tbend3 == 0.0 .OR. tbend3 == BLANK) tbend3 = tbend1
 IF (tbend5 == 0.0 .OR. tbend5 == BLANK) tbend5 = tbend1
 eltemp = est(52)
 theta1 = thetam*degra
 sinth  = SIN(theta1)
 costh  = COS(theta1)
 IF (ABS(sinth) <= 1.0E-06) sinth = 0.0
 
!     EVALUTE MATERIAL PROPERTIES
 
 matflg = 2
 matid  = matid1
 IF (matid1 <= 0) GO TO 670
 CALL mat (idele)
 
 matflg = 2
 matid  = matid2
 CALL mat (idele)
 d13    = em(3)
 d23    = em(5)
 d33    = em(6)
 j11    = 0.0
 j12    = 0.0
 j22    = 0.0
 IF (nots) GO TO 30
 matflg = 3
 matid  = matid3
 CALL mat (idele)
 j11    = 1.0/(rj11*tshr)
 j12    = 0.0
 j22    = 1.0/(rj22*tshr)
 30 CONTINUE
 
!     CALCULATIONS FOR THE TRIANGLE
 
 CALL trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,iest(1),NAME)
 
!     CALCULATIONS FOR QMATRIX (QQQ) AND ITS INVERSE
 
 DO  i = 1,20
   DO  j = 1,20
     qqq(i,j) = 0.0
   END DO
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
   i1 = (i-1)*3 + 1
   DO  j = 1,6
     q(i,j) = qqq(i1,j)
   END DO
 END DO
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (6,q,6,qqqinv(1),0,determ,ising,ind)
 IF (ising == 2) GO TO 660
 
!     FOURTH ARGUMENT IS A DUMMY LOCATION FOR INVERSE AND HENCE TS1(1)
!     IS U
 
!     AGAIN RESET ISING TO -1
 
 ising = -1
 CALL invers (20,qqq,20,qqqinv(1),0,determ,ising,INDEX)
 
!     ISING EQUAL TO 2 IMPLIES THAT QQQ IS SINGULAR
 
 IF (ising == 2) GO TO 660
 
!     FIRST 18 COLUMNS OF QQQ INVERSE IS THE QQQINV FOR USE IN STIFFNESS
!     MA CALCULATIONS
 
 DO  i = 1,20
   DO  j = 1,18
     ij = (i-1)*18 + j
     qqqinv(ij) = qqq(i,j)
   END DO
 END DO
 
!     START EXECUTION FOR STIFFNESS MATRIX CALCULATION
 
!     CM IS STIFFNESS MATRIX IN ELEMENT COORDINATES
 
!     OBTAIN MEMBRANE STRESSES
 
!     RELEVANT PORTION OF STRESS ROUTINE OF TRIM6 IS CODED HERE
 
!     TRANSFORMATION MATRIX BETWEEN ELEMENT AND BASIC COORDINATES
 
 es(1) = ivect(1)
 es(2) = jvect(1)
 es(3) = ivect(2)
 es(4) = jvect(2)
 es(5) = ivect(3)
 es(6) = jvect(3)
 DO  i = 1,9
   balotr(i) = 0.0
 END DO
 
 DO  i = 1,7
   ph1out(i) = est(i)
 END DO
 ph1out( 8) = est(10)
 ph1out( 9) = est(11)
 ph1out(10) = est(12)
 ph1out(11) = tref
 emod(1) = em(1)
 emod(2) = em(2)
 emod(3) = em(3)
 emod(4) = em(2)
 emod(5) = em(4)
 emod(6) = em(5)
 emod(7) = em(3)
 emod(8) = em(5)
 emod(9) = em(6)
 
 CALL gmmats (emod,3,3,0,alf(1),3,1,0,ph1out(228))
 DO  jj = 1,3
   j = 2*jj - 1
   x = xc(j)
   y = yc(j)
   DO  i = 1,36
     tmmm(i) = 0.0
   END DO
   
!     TM MATRIX IS THE PRODUCT OF B AND Q INVERSE MATRICES
   
   DO  j = 1,6
     j1 = (j-1)*2 + 1
     j2 = j1 + 1
     tm(1,j1) = q(2,j) + 2.0*x*q(4,j) + y*q(5,j)
     tm(2,j2) = q(3,j) + x*q(5,j) + 2.0*y*q(6,j)
     tm(3,j1) = tm(2,j2)
     tm(3,j2) = tm(1,j1)
   END DO
   
!     ZERO STRESS VECTOR STORAGE
   
   DO  i = 1,3
     stress(i) = 0.0
   END DO
   
   DO  ii = 1,6
     ij1 = (jj-1)*54 + (ii-1)*9 + 12
     IF (ics(ii) == 0) GO TO 140
     CALL transs (iest(4*ii+24),trans)
     CALL gmmats (es,3,2,+1,trans,3,3,0,ee1)
     GO TO 160
     140 CONTINUE
     DO  i = 1,3
       DO  j = 1,2
         i1 = (i-1)*2 + j
         j1 = (j-1)*3 + i
         ee1(j1) = es(i1)
       END DO
     END DO
     160 CONTINUE
     mz = (ii-1)*6 + 1
     CALL gmmats (emod,3,3,0,tmmm(mz),2,3,+1,eph1)
     CALL gmmats (eph1,3,2,0,ee1,2,3,0,ph1out(ij1))
     
!     POINTER TO I-TH SIL IN PH1OUT
     
     npoint = 55 + (ii-1)*6
     
!     POINTER TO  3X3 S SUB I MATRIX
     
     npt1 = 12 + (ii-1)*9 + (jj-1)*54
     
     CALL gmmats (ph1out(npt1),3,3,0,est(npoint),3,1,0,vec(1))
     DO  j = 1,3
       stress(j) = stress(j) + vec(j)
       str(j) = stress(j)
     END DO
   END DO
   IF (iest(54) == -1) GO TO 200
   tem = est(54) - ph1out(11)
   DO  i = 1,3
     stress(i) = stress(i) - ph1out(227+i)*tem
     str(i) = stress(i)
   END DO
   200 CONTINUE
   sigx(jj)  = stress(1)
   sigy(jj)  = stress(2)
   sigxy(jj) = stress(3)
 END DO
 
!     EL, FL, GL ARE COEFFICIENTS IN LINEAR VARIATION OF SIGX, SIGY,
!     SIGXY RESPECTIVELY OVER THE ELEMENT
 
 
 el(1) = (sigx(1)*a + sigx(2)*b)/(a+b)
 el(2) = (sigx(2) - sigx(1))/(a+b)
 el(3) = (sigx(3) - el(1))/c
 fl(1) = (sigy(1)*a + sigy(2)*b)/(a+b)
 fl(2) = (sigy(2) - sigy(1))/(a+b)
 fl(3) = (sigy(3) - fl(1))/c
 gl(1) = (sigxy(1)*a + sigxy(2)*b)/(a+b)
 gl(2) = (sigxy(2) - sigxy(1))/(a+b)
 gl(3) = (sigxy(3) - gl(1))/c
 
!     EVALUATE THE CONSTANTS C1,C2,AND C3 IN THE LINEAR EQUATION FOR
!     THICKNESS VARIATION
 
 CALL af (f,18,a,b,c,cab1,cab2,cab3,tmem1,tmem3,tmem5,0)
 cab(1) = cab1
 cab(2) = cab2
 cab(3) = cab3
 unimem =.false.
 uniben =.false.
 
 d334 = d33*4.0
 d132 = d13*2.0
 d232 = d23*2.0
 
!     A1,A2,A3 ARE THE COEFFICIENTS OF LINEAR EQUATION FOR VARIATION
!     OF BENDING THICKNESSES
 
 CALL af (f,18,a,b,c,a1,a2,a3,thk1,thk2,thk3,0)
 IF (ABS(cab2) <= 1.e-6 .AND. ABS(cab3) <= 1.e-6) unimem =.true.
 IF (ABS(a2) <= 1.0E-06 .AND. ABS(a3) <= 1.0E-06) uniben =.true.
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
 h4   = q(4,1)*zc(1) + q(4,2)*zc(2) + q(4,3)*zc(3) + q(4,4)*zc(4) +  &
     q(4,5)*zc(5) + q(4,6)*zc(6)
 h5   = q(5,1)*zc(1) + q(5,2)*zc(2) + q(5,3)*zc(3) + q(5,4)*zc(4) +  &
     q(5,5)*zc(5) + q(5,6)*zc(6)
 h6   = q(6,1)*zc(1) + q(6,2)*zc(2) + q(6,3)*zc(3) + q(6,4)*zc(4) +  &
     q(6,5)*zc(5) + q(6,6)*zc(6)
 h4   = h4*2.0
 h6   = h6*2.0
 
!     H5 IS MULTIPLIED BY 2.0, SO THAT EXY=DU/DY + DV/DX - ZXY*W
 
 h5   = h5*2.0
 
 DO  i = 1,32
   ix   = xu(i)
   rix  = ix
   jx   = yu(i)
   rjx  = jx
   kx   = xv(i)
   rkx  = kx
   lx   = yv(i)
   rlx  = lx
   mx   = xw(i)
   rmx  = mx
   nx   = yw(i)
   rnx  = nx
   rmnx = rmx*rnx
   rmx1 = rmx*(rmx-1.0)
   rnx1 = rnx*(rnx-1.0)
   
   DO  j = i,32
     ij   = (i-1)*32 + j
     ji   = (j-1)*32 + i
     iy   = xu(j)
     riy  = iy
     jy   = yu(j)
     rjy  = jy
     ky   = xv(j)
     rky  = ky
     ly   = yv(j)
     rly  = ly
     my   = xw(j)
     rmy  = my
     ny   = yw(j)
     rny  = ny
     rmny = rmy*rny
     rmy1 = rmy*(rmy-1.0)
     rny1 = rny*(rny-1.0)
     st   = 0.0
     DO  k = 1,3
       DO  l = 1,3
         rr    = rk(k) + rl(l)
         rr0   = rk(k) + rl(l) - 1
         rr1   = rk(k) + rl(l) + 1
         ss    = sk(k) + sl(l)
         ss0   = sk(k) + sl(l) - 1
         ss1   = sk(k) + sl(l) + 1
         mm    = mx + my
         mmrr0 = mm + rr0
         mmrr1 = mm + rr1
         nn    = nx + ny
         nnss1 = nn + ss1
         nnss0 = nn + ss0
         mmrr  = mm + rr
         nnss  = nn + ss
         kk    = kx + ky
         kkrr0 = kk + rr0
         ll    = lx + ly
         llss1 = ll + ss1
         ii    = ix + iy
         jj    = jx + jy
         iirr1 = ii + rr1
         jjss0 = jj + ss0
         ki    = kx + iy
         kirr  = ki + rr
         lj    = lx + jy
         ljss  = lj + ss
         ik    = ix + ky
         ikrr  = ik + rr
         jl    = jx + ly
         jlss  = jl + ss
         km    = kx + my
         kmrr  = km + rr
         ln    = lx + ny
         lnss1 = ln + ss1
         im    = ix + my
         imrr1 = im + rr1
         jn    = jx + ny
         jnss  = jn + ss
         kkrr  = kk + rr
         llss  = ll + ss
         kirr1 = ki + rr1
         ljss0 = lj + ss0
         mk    = mx + ky
         mkrr  = mk + rr
         nlss1 = nx + ly + ss1
         mi    = mx + iy
         mirr1 = mi + rr1
         nj    = nx + jy
         njss  = nj + ss
         ikrr0 = ik + rr0
         jlss1 = jl + ss1
         iirr  = ii + rr
         jjss  = jj + ss
         ikrr1 = ik + rr1
         jlss0 = jl + ss0
         lnss1 = ln + ss1
         kirr0 = ki + rr0
         ljss1 = lj + ss1
         sb1   = 0.0
         sb2   = 0.0
         sb3   = 0.0
         sb4   = 0.0
         sb5   = 0.0
         sb6   = 0.0
         sb7   = 0.0
         sb8   = 0.0
         sb9   = 0.0
         sb10  = 0.0
         sb11  = 0.0
         sb12  = 0.0
         sb13  = 0.0
         sb14  = 0.0
         sb15  = 0.0
         sb16  = 0.0
         sb17  = 0.0
         sb18  = 0.0
         sb19  = 0.0
         sb20  = 0.0
         sb21  = 0.0
         sb22  = 0.0
         sb23  = 0.0
         sb24  = 0.0
         sb25  = 0.0
         sb26  = 0.0
         sb27  = 0.0
         sb28  = 0.0
         sb29  = 0.0
         sb30  = 0.0
         sb31  = 0.0
         sb32  = 0.0
         sb33  = 0.0
         sb34  = 0.0
         sb35  = 0.0
         sb36  = 0.0
         sb37  = 0.0
         sb38  = 0.0
         sb39  = 0.0
         sb40  = 0.0
         IF (mmrr0 > 0) sb1 = cab(k)*el(l)*rmx*rmy*f(mmrr0,nnss1)
         IF (nnss0 > 0) sb2 = cab(k)*fl(l)*rnx*rny*f(mmrr1,nnss0)
         IF (mmrr > 0 .AND. nnss > 0) sb3 = cab(k)*gl(l)*rnx*rmy*  &
             f(mmrr,nnss)
         IF (mmrr > 0 .AND. nnss > 0) sb4 = cab(k)*gl(l)*rmx*rny*  &
             f(mmrr,nnss)
         IF (kkrr0 > 0) sb5 = cab(k)*el(l)*rkx*rky*f(kkrr0,llss1)
         IF (jjss0 > 0) sb6 = cab(k)*el(l)*rjx*rjy*f(iirr1,jjss0)
         IF (kirr > 0 .AND. ljss > 0) sb7 = cab(k)*el(l)*rkx*rjy*  &
             f(kirr,ljss)
         IF (ikrr > 0 .AND. jlss > 0) sb8 = cab(k)*el(l)*rjx*rky*  &
             f(ikrr,jlss)
         IF (kirr > 0 .AND. ljss > 0) sb9 = cab(k)*el(l)*rkx*rjy*  &
             f(kirr,ljss)
         IF (kkrr0 > 0) sb10 = cab(k)*el(l)*rkx*rky*f(kkrr0,llss1)
         IF (kmrr  > 0) sb11 = cab(k)*el(l)*rkx*h5*f(kmrr,lnss1)
         IF (jjss0 > 0) sb12 = cab(k)*el(l)*rjx*rjy*f(iirr1,jjss0)
         IF (ikrr > 0 .AND. jlss > 0) sb13 = cab(k)*el(l)*rjx*rky*  &
             f(ikrr,jlss)
         IF (jnss  > 0) sb14 = cab(k)*el(l)*rjx*h5*f(imrr1,jnss)
         IF (kkrr0 > 0) sb15 = cab(k)*fl(l)*rkx*rky*f(kkrr0,llss1)
         IF (kirr > 0 .AND. ljss > 0) sb16 = cab(k)*fl(l)*rkx*rjy*  &
             f(kirr,ljss)
         IF (jjss0 > 0) sb17 = cab(k)*fl(l)*rjx*rjy*f(iirr1,jjss0)
         IF (ikrr > 0 .AND. jlss > 0) sb18 = cab(k)*fl(l)*rjx*rky*  &
             f(ikrr,jlss)
         IF (kirr > 0 .AND. ljss > 0) sb19 = cab(k)*fl(l)*rkx*rjy*  &
             f(kirr,ljss)
         IF (kkrr0 > 0) sb20 = cab(k)*fl(l)*rkx*rky*f(kkrr0,llss1)
         IF (kmrr  > 0) sb21 = cab(k)*fl(l)*rkx*h5*f(kmrr,lnss1)
         IF (jjss0 > 0) sb22 = cab(k)*fl(l)*rjx*rjy*f(iirr1,jjss0)
         IF (ikrr > 0 .AND. jlss > 0) sb23 = cab(k)*fl(l)*rjx*rky*  &
             f(ikrr,jlss)
         IF (jnss  > 0) sb24 = cab(k)*fl(l)*rjx*h5*f(imrr1,jnss)
         IF (kkrr > 0 .AND. llss > 0) sb25 = cab(k)*gl(l)*rlx*rky*  &
             f(kkrr,llss)
         IF (kkrr > 0 .AND. llss > 0) sb26 = cab(k)*gl(l)*rkx*rly*  &
             f(kkrr,llss)
         IF (ljss0 > 0) sb27 = cab(k)*gl(l)*rlx*rjy*f(kirr1,ljss0)
         IF (jlss0 > 0) sb28 = cab(k)*gl(l)*rjx*rly*f(ikrr1,jlss0)
         IF (mkrr  > 0) sb29 = cab(k)*gl(l)*rky*h6*f(mkrr,nlss1)
         IF (kmrr  > 0) sb30 = cab(k)*gl(l)*rkx*h6*f(kmrr,lnss1)
         IF (njss  > 0) sb31 = cab(k)*gl(l)*rjy*h6*f(mirr1,njss)
         IF (jnss  > 0) sb32 = cab(k)*gl(l)*rjx*h6*f(imrr1,jnss)
         IF (ikrr0 > 0) sb33 = cab(k)*gl(l)*rix*rky*f(ikrr0,jlss1)
         IF (kirr0 > 0) sb34 = cab(k)*gl(l)*rkx*riy*f(kirr0,ljss1)
         IF (iirr > 0 .AND. jjss > 0) sb35 = cab(k)*gl(l)*rix*rjy*  &
             f(iirr,jjss)
         IF (iirr > 0 .AND. jjss > 0) sb36 = cab(k)*gl(l)*rjx*riy*  &
             f(iirr,jjss)
         IF (mkrr  > 0) sb37 = cab(k)*gl(l)*rky*h4*f(mkrr,nlss1)
         IF (kmrr  > 0) sb38 = cab(k)*gl(l)*rkx*h4*f(kmrr,lnss1)
         IF (njss  > 0) sb39 = cab(k)*gl(l)*rjy*h4*f(mirr1,njss)
         IF (jnss  > 0) sb40 = cab(k)*gl(l)*rjx*h4*f(imrr1,jnss)
         st = st + sb1 + sb2 + sb3 + sb4 +  &
             0.25*(sb5+sb6-sb7-sb8) + (sb9+sb10-sb11-sb12-sb13+sb14) +  &
             0.25*(sb15-sb16+sb17-sb18) + (-sb19-sb20+sb21+sb22+sb23-sb24)  &
             + 0.5*(sb25+sb26-sb27-sb28-sb29-sb30+sb31+sb32) +  &
             0.5*(-sb33-sb34+sb35+sb36+sb37+sb38-sb39-sb40)
       END DO
       IF (unimem) EXIT
     END DO
     240 CONTINUE
     kshl(ij) = st
     kshl(ji) = kshl(ij)
   END DO
 END DO
 
!     IF NO TRANSVERSE SHEAR GO TO 230
 
!     IF TSHR EQUAL TO ZERO OR MATID3 EQUAL TO ZERO , SKIP THESE
!     CALCULATION
 
 IF (nots) GO TO 270
 
!     CURRENTLY, TRANSVERSE SHEAR CALCULATIONS ARE NOT CODED FOR SHELL
!     ELEMENT WHEN IT IS CODED, CALL THE ROUTINE HERE
 
 270 CONTINUE
 CALL gmmats (q,6,6,0,kshl(1),6,32,0,qks(1))
 CALL gmmats (q,6,6,0,kshl(193),6,32,0,qks(193))
 CALL gmmats (qqqinv,20,18,+1,kshl(385),20,32,0,qks(385))
 DO  i = 1,30
   DO  j = 1,6
     ij = (i-1)*32 + j
     ji = (i-1)*6  + j
     kshl(    ji) = qks(  ij)
     kshl(180+ji) = qks(6+ij)
   END DO
 END DO
 DO  i = 1,30
   DO  j = 1,20
     ij = (i-1)*32 + j + 12
     ji = (i-1)*20 + j + 360
     kshl(ji) = qks(ij)
   END DO
 END DO
 CALL gmmats (kshl(  1),30,6 ,0,q,6,6,1 ,qks(  1))
 CALL gmmats (kshl(181),30,6 ,0,q,6,6,1 ,qks(181))
 CALL gmmats (kshl(361),30,20,0,qqqinv,20,18,0,qks(361))
 DO  i = 1,30
   DO  j = 1,6
     ij = (i-1)*30 + j
     ji = (i-1)*6  + j
     cms(ij  ) = qks(ji    )
     cms(ij+6) = qks(ji+180)
   END DO
 END DO
 DO  i = 1,30
   DO  j = 1,18
     ij = (i-1)*30 + j + 12
     ji = (i-1)*18 + j + 360
     cms(ij) = qks(ji)
   END DO
 END DO
 DO  i = 1,30
   ee(i) = 0.0D0
 END DO
 ee( 1) = ivect(1)
 ee( 2) = jvect(1)
 ee( 3) = kvect(1)
 ee( 6) = ivect(2)
 ee( 7) = jvect(2)
 ee( 8) = kvect(2)
 ee(11) = ivect(3)
 ee(12) = jvect(3)
 ee(13) = kvect(3)
 ee(19) = ivect(1)
 ee(20) = jvect(1)
 ee(24) = ivect(2)
 ee(25) = jvect(2)
 ee(29) = ivect(3)
 ee(30) = jvect(3)
 DO  k = 1,6
   DO  i = 1,2
     k1 = 6*(i-1) + k
     i1 = 5*(k-1) + i
     DO  j = 1,30
       ctm (i1,j) = cm1(k1,j)
     END DO
   END DO
 END DO
 DO  k = 1,6
   DO  i = 1,3
     i2 = 5*(k-1)  + i + 2
     k2 = 12 + (k-1)*3 + i
     DO  j = 1,30
       ctm (i2,j) = cm1(k2,j)
     END DO
   END DO
 END DO
 DO  k = 1,6
   DO  i = 1,2
     k1 = 6*(i-1) + k
     i1 = 5*(k-1) + i
     DO  j = 1,30
       cm1(j,i1) = ctm (j,k1)
     END DO
   END DO
 END DO
 DO  k = 1,6
   DO  i = 1,3
     i2 = 5*(k-1)  + i + 2
     k2 = 12 + (k-1)*3 + i
     DO  j = 1,30
       cm1(j,i2) = ctm(j,k2)
     END DO
   END DO
 END DO
 
!     LOCATE THE TRANSFORMATION MATRICES FROM BASIC TO LOCAL (THAT IS
!     COORDINATE AT ANY GRID POINT IN WHICH DISPLACEMENT AND STRESSES
!     ARE R
!     - NOT NEEDED IF FIELD 7 IN GRID CARD IS ZERO)
 
!     TRANSFORM STIFFNESS MATRIX FROM ELEMENT COORDINATES TO BASIC
!     COORDINATE
 
!     TRANSFORM STIFFNESS MATRIX FROM BASIC COORDINATES TO GLOBAL (DISP)
!     COORDINATES
 
!     INSERT THE 6X6 SUBMATRIX  INTO KGG MATRIX
 
 DO  i = 1,1296
   cmt(i) = 0.0
 END DO
 DO  i = 1,6
   sil(i) = i
 END DO
 DO  i = 1,6
   IF (npvt /= iest(i+1)) CYCLE
   npivot = i
   GO TO  520
 END DO
 nogo = .true.
 WRITE (ioutpt,720) sfm,iest(1)
 RETURN
 
 520 CONTINUE
 i    = npivot
 sil1 = sil(npivot)
 DO  j = 1,6
   sil2 = sil(j)
   DO  ii = 1,36
     balotr(ii) = 0.0
     ksub(ii)   = 0.0
   END DO
   DO  k = 1,5
     k1 = (sil1-1)*5 + k
     DO  l = 1,5
       l1 = (sil2-1)*5 + l
       csub(k,l) = cm1(k1,l1)
     END DO
   END DO
   CALL gmmats (ee,6,5,0,csub,5,5,0,csubt)
   CALL gmmats (csubt,6,5,0,ee,6,5,+1,ksubt)
   DO  k = 1,6
     DO  l = 1,6
       k1 = (k-1)*6 + l
       l1 = (l-1)*6 + k
       ksub(l1) = ksubt(k1)
     END DO
   END DO
   
!     TRANSFORM THE KSUB(36) FROM BASIC TO DISPLACEMENT COORDINATES
   
   IF (nl(sil1) == 0 .OR. ics(sil1) == 0) GO TO 590
   CALL transs (iest(4*sil1+24),trand)
   DO  jj = 1,3
     l = 6*(jj-1) + 1
     m = 3*(jj-1) + 1
     balotr(l   ) = trand(m  )
     balotr(l+1 ) = trand(m+1)
     balotr(l+2 ) = trand(m+2)
     balotr(l+21) = trand(m  )
     balotr(l+22) = trand(m+1)
     balotr(l+23) = trand(m+2)
   END DO
   CALL gmmats (balotr(1),6,6,1,ksub(1),6,6,0,ksubt)
   DO  k = 1,36
     ksub(k) = ksubt(k)
   END DO
   590 CONTINUE
   IF (nl(sil2) == 0 .OR. ics(sil2) == 0) GO TO 630
   IF (j == i) GO TO 610
   CALL transs (iest(4*sil2+24),trand)
   DO  jj = 1,3
     l = 6*(jj-1) + 1
     m = 3*(jj-1) + 1
     balotr(l   ) = trand(m  )
     balotr(l+1 ) = trand(m+1)
     balotr(l+2 ) = trand(m+2)
     balotr(l+21) = trand(m  )
     balotr(l+22) = trand(m+1)
     balotr(l+23) = trand(m+2)
   END DO
   610 CONTINUE
   CALL gmmats (ksub(1),6,6,0,balotr(1),6,6,0,ksubt)
   DO  k = 1,36
     ksub(k) = ksubt(k)
   END DO
   630 CONTINUE
   DO  ijk = 1,36
     ksubd(ijk) = DBLE(ksub(ijk))
   END DO
   CALL ds1b (ksubd(1),iest(j+1))
 END DO
 GO TO 730
 660 CONTINUE
 nogo =.true.
 WRITE (ioutpt,700) ufm,iest(1)
 RETURN
 
 670 CONTINUE
 nogo =.true.
 WRITE (ioutpt,710) ufm,iest(1)
 RETURN
 
 700 FORMAT (a23,' 2416, MATRIX RELATING GENERALIZED PARAMETERS AND ',  &
     'GRID POINT DISPLACEMENTS IS SINGULAR.', /26X,  &
     'CHECK COORDINATES OF ELEMENT  TRSHL WITH ID =',i9,1H.)
 710 FORMAT (a23,' 2418, MATERIAL ID FOR MEMBRANE EFFECTS IS LESS ',  &
     'THAN OR EQUAL TO ZERO FOR TRSHL ELEMENT WITH ID =',i9,1H.)
 720 FORMAT (a25,' 2419, PIVOT POINT IS NOT EQUAL TO TRSHL ELEMENT ',  &
     'GRID POINTS FOR ELEMENT ID =',i9,1H.)
 730 CONTINUE
 RETURN
END SUBROUTINE dtshls
