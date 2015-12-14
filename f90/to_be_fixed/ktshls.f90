SUBROUTINE ktshls
     
!     ECPT ENTRIES
 
!     ECPT( 1) = ELEMENT ID                                     INTEGER
!     ECPT( 2) = SCALAR INDEX NUMBER FOR GRID POINT 1           INTEGER
!     ECPT( 3) = SCALAR INDEX NUMBER FOR GRID POINT 2           INTEGER
!     ECPT( 4) = SCALAR INDEX NUMBER FOR GRID POINT 3           INTEGER
!     ECPT( 5) = SCALAR INDEX NUMBER FOR GRID POINT 4           INTEGER
!     ECPT( 6) = SCALAR INDEX NUMBER FOR GRID POINT 5           INTEGER
!     ECPT( 7) = SCALAR INDEX NUMBER FOR GRID POINT 6           INTEGER
!     ECPT( 8) = THETA                                          REAL
!     ECPT( 9) = MATERIAL ID 1                                  INTEGER
!     ECPT(10) = THICKNESS T1 AT GRID POINT G1
!     ECPT(11) = THICKNESS T3 AT GRID POINT G3
!     ECPT(12) = THICKNESS T5 AT GRID POINT G5
!     ECPT(13) = MATERIAL ID 2                                  INTEGER
!     ECPT(14) = THICKNESS TBEND1 FOR BENDING AT GRID POINT G1
!     ECPT(15) = THICKNESS TBEND3 FOR BENDING AT GRID POINT G3
!     ECPT(16) = THICKNESS TBEND5 FOR BENDING AT GRID POINT G5
!     ECPT(17) = MATERIAL ID 3                                  INTEGER
!     ECPT(18) = THICKNESS TSHR1 FOR TRANSVERSE SHEAR AT GRID POINT G1
!     ECPT(19) = THICKNESS TSHR3 FOR TRANSVERSE SHEAR AT GRID POINT G3
!     ECPT(20) = THICKNESS TSHR5 FOR TRANSVERSE SHEAR AT GRID POINT G5
!     ECPT(21) = NON-STRUCTURAL MASS                            REAL
!     ECPT(22) = DISTANCE Z11 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT(23) = DISTANCE Z21 FOR STRESS CALCULATION  AT GRID POINT G1
!     ECPT(24) = DISTANCE Z13 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT(25) = DISTANCE Z23 FOR STRESS CALCULATION  AT GRID POINT G3
!     ECPT(26) = DISTANCE Z15 FOR STRESS CALCULATION  AT GRID POINT G5
!     ECPT(27) = DISTANCE Z25 FOR STRESS CALCULATION  AT GRID POINT G5
 
!     X1,Y1,Z1 FOR ALL SIX POINTS ARE  IN NASTRAN BASIC SYSTEM
 
!     ECPT(28) = COORDINATE SYSTEM ID FOR GRID A                INTEGER
!     ECPT(29) = COORDINATE X1                                  REAL
!     ECPT(30) = COORDINATE Y1                                  REAL
!     ECPT(31) = COORDINATE Z1                                  REAL
!     ECPT(32) = COORDINATE SYSTEM ID FOR GRID B                INTEGER
!     ECPT(33) = COORDINATE X1                                  REAL
!     ECPT(34) = COORDINATE Y1                                  REAL
!     ECPT(35) = COORDINATE Z1                                  REAL
!     ECPT(36) = COORDINATE SYSTEM ID FOR GRID C                INTEGER
!     ECPT(37) = COORDINATE X1                                  REAL
!     ECPT(38) = COORDINATE Y1                                  REAL
!     ECPT(39) = COORDINATE Z1                                  REAL
!     ECPT(40) = COORDINATE SYSTEM ID FOR GRID D                INTEGER
!     ECPT(41) = COORDINATE X1                                  REAL
!     ECPT(42) = COORDINATE Y1                                  REAL
!     ECPT(43) = COORDINATE Z1                                  REAL
!     ECPT(44) = COORDINATE SYSTEM ID FOR GRID E                INTEGER
!     ECPT(45) = COORDINATE X1                                  REAL
!     ECPT(46) = COORDINATE Y1                                  REAL
!     ECPT(47) = COORDINATE Z1                                  REAL
!     ECPT(48) = COORDINATE SYSTEM ID FOR GRID F                INTEGER
!     ECPT(49) = COORDINATE X1                                  REAL
!     ECPT(50) = COORDINATE Y1                                  REAL
!     ECPT(51) = COORDINATE Z1                                  REAL
!     EST (52) = ELEMENT TEMPERATURE
 
 LOGICAL :: imass,nots,nogo,unimem,uniben
 INTEGER :: xu(32),yu(32),xv(32),yv(32),xw(32),yw(32),  &
     sil(6),sil1,sil2,SAVE(6),xthk(10),ythk(10),  &
     rk(3),sk(3),eltype,elid,estid,dict(15),small(6)
!                     RK AND SK ARE EXPONENTS IN THICKNESS VARIATION
 REAL :: j11,j12,j22,nsm,mshl(1024),kshl(1024),ksub(6,6),  &
     ksubt(6,6),ivect(3),jvect(3),kvect(3),ksup,ksupt
 DIMENSION       INDEX(20,3),ics(6),iest(42),nl(6),ind(6,3),  &
     qqqinv(360),xc(6),yc(6),zc(6),qqq(20,20),  &
     cmt(1296),cms(900),cm1(30,30),trand(9),balotr(36),  &
     cc(10),cab(3),qks(960),ksup(36),ksupt(36), ctm(36,36),NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nok,nom,nob
 COMMON /emgest/ est(100)
 COMMON /emgdic/ eltype,ldict,nlocs,elid,estid
 COMMON /sma1dp/ f(14,14),q(6,6),ee(30),csubt(6,5),csub(5,5)
 COMMON /sma2dp/ trand,balotr,ksub,ksubt,fac,xc,yc,zc,ivect,jvect,  &
     kvect,cc,cab,dict,sil,SAVE,small,INDEX,ics,nl
 COMMON /sma1cl/ kdummy(22),knogo
 COMMON /emgprm/ ixtra,izr,nzr,dumy(12),kmbgg(3),iprec,nogo
 COMMON /system/ ksystm(65)
 COMMON /matin / matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/ em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy,  &
     rj11,rj12,rj22
 
!     SMA1 WORKING STORAGE
 
!     EQUIVALENCE IECPT WITH ECPT IN COMMON BLOCK /SMA1ET/ SINCE ECPT IS
!     A MIXED INTEGER AND REAL ARRAY
 
 EQUIVALENCE    (c1,cc(1)), (c2,cc(2)), (c3,cc(3)), (c4,cc(4)),  &
     (c5,cc(5)), (c6,cc(6)), (c7,cc(7)), (c8,cc(8)), (c9,cc(9)), (c10,cc(10)),  &
     (ksub(1,1),ksup(1)), (ksubt(1,1),ksupt(1)),  &
     (cmt(1),ctm(1,1)),   (qks(1),cmt(1025))
 EQUIVALENCE    (a,dista), (b,distb), (c,distc), (iest(1),est(1))
 EQUIVALENCE    (cmt(1),kshl(1),mshl(1),qqq(1,1))
 EQUIVALENCE    (ksystm(2),ioutpt)
 EQUIVALENCE    (thk1,tbend1), (thk2,tbend3), (thk3,tbend5)
 EQUIVALENCE    (cm1(1,1),cms(1)), (ind(1,1),INDEX(1,1))
 DATA   xu    / 0,1,0,2,1,0,26*0     /, yu    / 0,0,1,0,1,2,26*0     /,  &
     xv    / 6*0,0,1,0,2,1,0,20*0 /, yv    / 6*0,0,0,1,0,1,2,20*0 /,  &
     xw    / 12*0,0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/  &
     yw    / 12*0,0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/
 DATA   BLANK , NAME  / 4H    , 4HTRSH, 4HL     /
 DATA   rk    / 0,1,0 /
 DATA   sk    / 0,0,1 /
 DATA   degra / 0.0174532925        /
 DATA   xthk  / 0,1,0,2,1,0,3,2,1,0 /
 DATA   ythk  / 0,0,1,0,1,2,0,1,2,3 /
 
 dict(1) = estid
 
!     COMPONENT CODE,ICODE,IS  111111  AND HAS A VALUE OF 63
 
 icode  = 63
 ndof   = 36
 nsq    = ndof**2
 dict(2)= 1
 dict(3)= ndof
 dict(4)= icode
 dict(5)= gsube
 nots   =.false.
 imass  =.false.
 IF (nom > 0) imass =.true.
 ipass  = 1
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
 DO  i = 28,48,4
   j      = j + 1
   ics(j) = iest(i )
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
 IF (tshr1 == 0.0) nots=.true.
 tshr = (tshr1 + tshr3 + tshr5)/3.0
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
 IF  (matid1 == 0) GO TO 30
 CALL mat (idele)
 
 g11 = em(1)
 g12 = em(2)
 g13 = em(3)
 g22 = em(4)
 g23 = em(5)
 g33 = em(6)
 30 CONTINUE
 matflg = 2
 matid  = matid2
 IF (matid2 == 0) GO TO 40
 CALL mat (idele)
 d11 = em(1)
 d12 = em(2)
 d13 = em(3)
 d22 = em(4)
 d23 = em(5)
 d33 = em(6)
 j11 = 0.0
 j12 = 0.0
 j22 = 0.0
 IF (nots) GO TO 40
 matflg = 3
 matid  = matid3
 CALL mat (idele)
 j11 = 1.0/(rj11*tshr)
 j12 = 0.0
 j22 = 1.0/(rj22*tshr)
 40 CONTINUE
 
!     CALCULATIONS FOR THE TRIANGLE
 
 CALL trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,iest(1),NAME)
 
!     COMPUTE THE AREA INTEGRATION FUNCTION F
 
 CALL af (f,14,a,b,c,0,0,0,0,0,0,-1)
 
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
   DO  j = 1,6
     i1 = (i-1)*3 + 1
     q(i,j) = qqq(i1,j)
   END DO
 END DO
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (6,q,6,balotr(1),0,determ,ising,ind)
 IF (ising == 2) GO TO 700
 
!     FOURTH ARGUMENT IS A DUMMY LOCATION FOR INVERSE AND HENCE TS1(1)
!     IS U
 
 ising = -1
 CALL invers (20,qqq,20,balotr(1),0,determ,ising,INDEX)
 
!     ISING EQUAL TO 2 IMPLIES THAT QQQ IS SINGULAR
 
 IF (ising == 2) GO TO 700
 
!     FIRST 18 COLUMNS OF QQQ INVERSE IS THE QQQINV FOR USE IN STIFFNESS
!     CALCULATIONS
 
 
 DO  i = 1,20
   DO  j = 1,18
     ijk = (i-1)*18 + j
     qqqinv(ijk) = qqq(i,j)
   END DO
 END DO
 
!     START EXECUTION FOR STIFFNESS MATRIX CALCULATION
 
!     CM IS STIFFNESS MATRIX IN ELEMENT COORDINATES
 
 90 CONTINUE
 
!     EVALUATE THE CONSTANTS C1,C2,AND C3 IN THE LINEAR EQUATION FOR
!     THICKNESS VARIATION - MEMBRANE
 
 CALL af (f,14,a,b,c,c1,c2,c3,tmem1,tmem3,tmem5,1)
 cab(1) = c1
 cab(2) = c2
 cab(3) = c3
 area = f(1,1)
 vol  = c1*f(1,1) + c2*f(2,1) + c3*f(1,2)
 
 
 d334 = d33*4.0
 d132 = d13*2.0
 d232 = d23*2.0
 
!     A1,A2,A3 ARE THE COEFFICIENTS OF LINEAR EQUATION FOR VARIATION
!     OF BENDING THICKNESSES
 
 CALL af (f,14,a,b,c,a1,a2,a3,thk1,thk2,thk3,1)
 unimem =.false.
 uniben =.false.
 IF (ABS(c2) <= 1.0E-06 .AND. ABS(c3) <= 1.0E-06) unimem =.true.
 IF (ABS(a2) <= 1.0E-06 .AND. ABS(a3) <= 1.0E-06) uniben =.true.
 a1sq= a1*a1
 a2sq= a2*a2
 a3sq= a3*a3
 c1  = a1sq*a1
 c2  = 3.0*a1sq*a2
 c3  = 3.0*a1sq*a3
 c4  = 3.0*a1*a2sq
 c5  = 6.0*a1*a2*a3
 c6  = 3.0*a3sq*a1
 c7  = a2sq*a2
 c8  = 3.0*a2sq*a3
 c9  = 3.0*a2*a3sq
 c10 = a3*a3sq
 
!     AA1, AA2, AA3  ARE COEFFICIENTS IN THICKNESS VARIATION FOR
!     TRANSVERSE SHEAR
 
 
!    (POSSIBLY AN ERROR HERE - AA1,AA2, AND AA3 ARE NOT USED IN PROGRAM)
!     CALL AF (F,14,A,B,C,AA1,AA2,AA3,TSHR1,TSHR3,TSHR5,1)
 
 h4 = q(4,1)*zc(1) + q(4,2)*zc(2) + q(4,3)*zc(3) + q(4,4)*zc(4) +  &
     q(4,5)*zc(5) + q(4,6)*zc(6)
 h5 = q(5,1)*zc(1) + q(5,2)*zc(2) + q(5,3)*zc(3) + q(5,4)*zc(4) +  &
     q(5,5)*zc(5) + q(5,6)*zc(6)
 h6 = q(6,1)*zc(1) + q(6,2)*zc(2) + q(6,3)*zc(3) + q(6,4)*zc(4) +  &
     q(6,5)*zc(5) + q(6,6)*zc(6)
 h4 = h4*2.0
 h6 = h6*2.0
 
!     H5 IS MULTIPLIED BY 2.0, SO THAT EXY=DU/DY + DV/DX - ZXY*W
 
 h5 = h5*2.0
 
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
   ixp1 = ix + 1
   jxp1 = jx + 1
   kxp1 = kx + 1
   lxp1 = lx + 1
   mxp1 = mx + 1
   nxp1 = nx + 1
   DO  j = i,32
     ij  = (i-1)*32 + j
     ji  = (j-1)*32 + i
     iy  = xu(j)
     riy = iy
     jy  = yu(j)
     rjy = jy
     ky  = xv(j)
     rky = ky
     ly  = yv(j)
     rly = ly
     my  = xw(j)
     rmy = my
     ny  = yw(j)
     rny = ny
     rmny= rmy*rny
     rmy1= rmy*(rmy-1.0)
     rny1= rny*(rny-1.0)
     mx0 = mx + my
     mx1 = mx + my - 1
     mx2 = mx + my - 2
     mx3 = mx + my - 3
     nx0 = nx + ny
     nx1 = nx + ny - 1
     nx2 = nx + ny - 2
     nx3 = nx + ny - 3
     my1 = mx + my + 1
     ny1 = nx + ny + 1
     ix0 = ix + iy
     ix1 = ix0- 1
     ix01= ix0+ 1
     jx0 = jx + jy
     jx1 = jx0- 1
     jx01= jx0+ 1
     kx0 = kx + ky
     kx1 = kx0- 1
     kx01= kx0+ 1
     lx0 = lx + ly
     lx1 = lx0- 1
     lx01= lx0+ 1
     IF (ipass == 1) GO TO 110
     ix011 = ix01 + 1
     jx011 = jx01 + 1
     rho   = rhoy*1.0
     IF (j > 12) GO TO 100
     mshl(ij) = rho*(cab(1)*f(ix01,jx01) + cab(2)*f(ix011,jx01) +  &
         cab(3)*f(ix01,jx011)) + nsm*f(ix01,jx01)
     mshl(ji) = mshl (ij)
     100 CONTINUE
     mx01 = mx0  + 1
     nx01 = nx0  + 1
     mx011= mx01 + 1
     nx011= nx01 + 1
     mshl(ij) = rho*(a1*f(mx01,nx01) + a2*f(mx011,nx01) +  &
         a3*f(mx01,nx011)) + nsm*f(mx01,nx01)
     mshl(ji) = mshl(ij)
     GO TO 210
     110 CONTINUE
     st = 0.0
     IF (i <= 12 .AND. j > 12) GO TO 160
     IF (i > 12) GO TO 140
     DO  k = 1,3
       ixr1  = ix1  + rk(k)
       jxs01 = jx01 + sk(k)
       lxs1  = lx1  + sk(k)
       kxr01 = kx01 + rk(k)
       ixr01 = ix01 + rk(k)
       jxs1  = jx1  + sk(k)
       kxr1  = kx1  + rk(k)
       lxs01 = lx01 + sk(k)
       iykx1 = iy + kx + rk(k)
       jylx1 = jy + lx + sk(k)
       ixky1 = ix + ky + rk(k)
       jxly1 = jx + ly + sk(k)
       ixiy0 = ix + iy + rk(k)
       jxjy0 = jx + jy + sk(k)
       iykx2 = iykx1 - 1
       jylx0 = jylx1 + 1
       ixky2 = ixky1 - 1
       jxly0 = jxly1 + 1
       kxky0 = kx + ky + rk(k)
       lxly0 = lx + ly + sk(k)
       ixky0 = ix + ky + rk(k) + 1
       jxly2 = jxly1 - 1
       iykx0 = iy + kx + rk(k) + 1
       jylx2 = jylx1 - 1
       st11  = 0.0
       st22  = 0.0
       st331 = 0.0
       st332 = 0.0
       st121 = 0.0
       st122 = 0.0
       st131 = 0.0
       st132 = 0.0
       st133 = 0.0
       st231 = 0.0
       st232 = 0.0
       st233 = 0.0
       IF (ixr1 > 0) st11  = g11*rix*riy*f(ixr1,jxs01)
       IF (lxs1 > 0) st22  = g22*rlx*rly*f(kxr01,lxs1)
       IF (jxs1 > 0) st331 = g33*rjx*rjy*f(ixr01,jxs1)
       IF (kxr1 > 0) st332 = g33*rkx*rky*f(kxr1,lxs01)
       IF (ixky1 > 0 .AND. jxly1 > 0) st121 = (g33*rjx*rky +  &
           g12*rix*rly)*f(ixky1,jxly1)
       IF (iykx1 > 0 .AND. jylx1 > 0) st122 = (g33*rjy*rkx +  &
           g12*riy*rlx)*f(iykx1,jylx1)
       IF (ixiy0 > 0 .AND. jxjy0 > 0) st131 = g13*(riy*rjx +  &
           rix*rjy)*f(ixiy0,jxjy0)
       IF (iykx2 > 0) st132 = g13*riy*rkx*f(iykx2,jylx0)
       IF (ixky2 > 0) st133 = g13*rix*rky*f(ixky2,jxly0)
       IF (kxky0 > 0 .AND. lxly0 > 0) st231 = g23*(rkx*rly +  &
           rky*rlx)*f(kxky0,lxly0)
       IF (jxly2 > 0) st232 = g23*rjx*rly*f(ixky0,jxly2)
       IF (jylx2 > 0) st233 = g23*rjy*rlx*f(iykx0,jylx2)
       
       st1 = (st11  + st22  + st331 + st332 + st121 + st122 + st131 +  &
           st132 + st133 + st231 + st232 + st233)* cab(k)
       st  = st + st1
       IF (unimem) EXIT
     END DO
     130 CONTINUE
     GO TO 200
     140 CONTINUE
     st = 0.0
     DO  k = 1,10
       mx3x = mx3 + xthk(k)
       ny1y = ny1 + ythk(k)
       my1x = my1 + xthk(k)
       nx3y = nx3 + ythk(k)
       mx1x = mx1 + xthk(k)
       nx1y = nx1 + ythk(k)
       mx2x = mx2 + xthk(k)
       nx0y = nx0 + ythk(k)
       mx0x = mx0 + xthk(k)
       nx2y = nx2 + ythk(k)
       s11  = 0.0
       s22  = 0.0
       s33  = 0.0
       s13  = 0.0
       s23  = 0.0
       IF (mx3x > 0) s11 = d11*rmx1*rmy1*cc(k)*f(mx3x,ny1y)
       IF (nx3y > 0) s22 = d22*rnx1*rny1*cc(k)*f(my1x,nx3y)
       IF (mx1x > 0 .AND. nx1y > 0) s33 = (d334*rmnx*rmny +  &
           d12*(rmx1*rny1 + rmy1*rnx1))*cc(k)*f(mx1x,nx1y)
       IF (mx2x > 0 .AND. nx0y > 0) s13 = d132*(rmx1*rmny +  &
           rmnx*rmy1)*cc(k)*f(mx2x,nx0y)
       IF (mx0x > 0 .AND. nx2y > 0) s23 = d232*(rmnx*rny1 +  &
           rnx1*rmny)*cc(k)*f(mx0x,nx2y)
       st = st + (s11 + s22 + s33 + s13 + s23)/12.0
       IF (uniben) EXIT
     END DO
     160 CONTINUE
     sb 7 = 0.0
     sb 9 = 0.0
     sb10 = 0.0
     sb18 = 0.0
     sb21 = 0.0
     sb26 = 0.0
     sb28 = 0.0
     sb31 = 0.0
     sb36 = 0.0
     sb38 = 0.0
     DO  k = 1,3
       ixmyr = ix + my + rk(k)
       jxnys1= jx + ny + sk(k) + 1
       sb1  = 0.0
       sb2  = 0.0
       sb3  = 0.0
       sb4  = 0.0
       sb5  = 0.0
       sb6  = 0.0
       sb8  = 0.0
       sb11 = 0.0
       sb12 = 0.0
       sb13 = 0.0
       sb14 = 0.0
       sb15 = 0.0
       sb16 = 0.0
       sb17 = 0.0
       sb19 = 0.0
       sb20 = 0.0
       sb22 = 0.0
       sb23 = 0.0
       sb24 = 0.0
       sb25 = 0.0
       sb27 = 0.0
       sb29 = 0.0
       sb30 = 0.0
       sb32 = 0.0
       sb33 = 0.0
       sb34 = 0.0
       sb35 = 0.0
       sb37 = 0.0
       sb39 = 0.0
       sb40 = 0.0
       IF (ixmyr > 0)  sb 1 =-g11*rix*h4*cab(k)*f(ixmyr,jxnys1)
       iymxr  = iy + mx + rk(k)
       jynxs1 = jy + nx + sk(k) + 1
       IF (iymxr > 0)  sb 2 =-g11*riy*h4*cab(k)*f(iymxr,jynxs1)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb 3 = g11*h4**2*cab(k)*f(mxmyr1,nxnys1)
       kxmyr1 = kx + my + rk(k) + 1
       lxnys  = lx + ny + sk(k)
       IF (lxnys > 0)  sb 4 =-g22*rlx*h6*cab(k)*f(kxmyr1,lxnys)
       mxkyr1 = mx + ky + rk(k) + 1
       nxlys  = nx + ly + sk(k)
       IF (nxlys > 0)  sb 5 =-g22*rly*h6*cab(k)*f(mxkyr1,nxlys)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb 6 = g22*h6**2*cab(k)*f(mxmyr1,nxnys1)
       ixmyr1 = ix + my + rk(k) + 1
       jxnys  = jx + ny + sk(k)
       IF(jxnys > 0)  sb 8 =-g33*rjx*h5*cab(k)*f(ixmyr1,jxnys)
       kxmyr  = kx + my + rk(k)
       lxnys1 = lx + ny + sk(k) + 1
       IF (kxmyr > 0)  sb11 =-g33*rkx*h5*cab(k)*f(kxmyr,lxnys1)
       mxiyr1 = mx + iy + rk(k) + 1
       nxjys  = nx + jy + sk(k)
       IF (nxjys > 0)  sb12 =-g33*rjy*h5*cab(k)*f(mxiyr1,nxjys)
       mxkyr  = mx + ky + rk(k)
       nxlys1 = nx + ly + sk(k) + 1
       IF (mxkyr > 0)  sb13 =-g33*rky*h5*cab(k)*f(mxkyr,nxlys1)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb14 = g33*h5**2*cab(k)*f(mxmyr1,nxnys1)
       ixmyr  = ix + my + rk(k)
       jxnys1 = jx + ny + sk(k) + 1
       IF (ixmyr > 0)  sb15 =-g12*rix*h6*cab(k)*f(ixmyr,jxnys1)
       mxkyr1 = mx + ky + rk(k) + 1
       nxlys  = nx + ly + sk(k)
       IF (nxlys > 0)  sb16 =-g12*rly*h4*cab(k)*f(mxkyr1,nxlys)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb17 = 2*g12*h4*h6*cab(k)*f(mxmyr1,nxnys1)
       kxmyr1 = kx + my + rk(k) + 1
       lxnys  = lx + ny + sk(k)
       IF (lxnys > 0)  sb19 =-g12*rlx*h4*cab(k)*f(kxmyr1,lxnys)
       mxiyr  = mx + iy + rk(k)
       nxjys1 = nx + jy + sk(k) + 1
       IF (mxiyr > 0)  sb20 =-g12*riy*h6*cab(k)*f(mxiyr,nxjys1)
       ixmyr  = ix + my + rk(k)
       jxnys1 = jx + ny + sk(k) + 1
       IF (ixmyr > 0)  sb22 =-g13*rix*h5*cab(k)*f(ixmyr,jxnys1)
       mxiyr1 = mx + iy + rk(k) + 1
       nxjys  = nx + jy + sk(k)
       IF (nxjys > 0)  sb23 =-g13*rjy*h4*cab(k)*f(mxiyr1,nxjys)
       mxkyr  = mx + ky + rk(k)
       nxlys1 = nx + ly + sk(k) + 1
       IF (mxkyr > 0)  sb24 =-g13*rky*h4*cab(k)*f(mxkyr,nxlys1)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb25 = 2*g13*h4*h5*cab(k)*f(mxmyr1,nxnys1)
       ixmyr1 = ix + my + rk(k) + 1
       jxnys  = jx + ny + sk(k)
       IF (jxnys > 0)  sb27 =-g13*rjx*h4*cab(k)*f(ixmyr1,jxnys)
       kxmyr  = kx + my + rk(k)
       lxnys1 = lx + ny + sk(k) + 1
       IF (kxmyr > 0)  sb29 =-g13*rkx*h4*cab(k)*f(kxmyr,lxnys1)
       mxiyr  = mx + iy + rk(k)
       nxjys1 = nx + jy + sk(k) + 1
       IF (mxiyr > 0)  sb30 =-g13*riy*h5*cab(k)*f(mxiyr,nxjys1)
       kxmyr1 = kx + my + rk(k) + 1
       lxnys  = lx + ny + sk(k)
       IF (lxnys > 0)  sb32 =-g23*rlx*h5*cab(k)*f(kxmyr1,lxnys)
       mxiyr1 = mx + iy + rk(k) + 1
       nxjys  = nx + jy + sk(k)
       IF (nxjys > 0)  sb33 =-g23*rjy*h6*cab(k)*f(mxiyr1,nxjys)
       mxkyr  = mx + ky + rk(k)
       nxlys1 = nx + ly + sk(k) + 1
       IF (mxkyr > 0)  sb34 =-g23*rky*h6*cab(k)*f(mxkyr,nxlys1)
       mxmyr1 = mx + my + rk(k) + 1
       nxnys1 = nx + ny + sk(k) + 1
       sb35 = 2*g23*h5*h6*cab(k)*f(mxmyr1,nxnys1)
       ixmyr1 = ix + my + rk(k) + 1
       jxnys  = jx + ny + sk(k)
       IF (jxnys > 0)  sb37 =-g23*rjx*h6*cab(k)*f(ixmyr1,jxnys)
       kxmyr  = kx + my + rk(k)
       lxnys1 = lx + ny + sk(k) + 1
       IF (kxmyr > 0)  sb39 =-g23*rkx*h6*cab(k)*f(kxmyr,lxnys1)
       mxkyr1 = mx + ky + rk(k) + 1
       nxlys  = nx + ly + sk(k)
       IF (nxlys. GT. 0)  sb40 =-g23*rly*h5*cab(k)*f(mxkyr1,nxlys)
       sb41 = sb3 + sb6 + sb14 + sb17 + sb25 + sb35
       IF (i <= 12) sb41 = 0.0
       st =  st   + sb1  + sb2  + sb4  + sb5  + sb7  + sb8  + sb9  +  &
           sb10 + sb11 + sb12 + sb13 + sb15 + sb16 + sb18 + sb19 +  &
           sb20 + sb21 + sb22 + sb23 + sb24 + sb26 + sb27 + sb28 +  &
           sb29 + sb30 + sb31 + sb32 + sb33 + sb34 + sb36 + sb37 +  &
           sb38 + sb39 + sb40 + sb41
       IF (unimem) EXIT
     END DO
     190 CONTINUE
     200 CONTINUE
     kshl(ij) = st
     kshl(ji) = kshl(ij)
     210 CONTINUE
   END DO
 END DO
 IF (ipass == 2) GO TO 240
 
!    CURRENTLY,TRANSVERSE SHEAR CALCULATIONS ARE NOT CODED FOR SHELL
!    ELEMENT WHEN IT IS CODED, CALL THE ROUTINE HERE
 
 240 CONTINUE
 
!     (QQQINV) TRANSPOSE (KTR3)  (QQQINV)
 
 CALL gmmats (q,6,6,0, kshl(1),6,32,0, qks(1))
 CALL gmmats (q,6,6,0, kshl(193),6,32,0, qks(193))
 CALL gmmats (qqqinv,20,18,+1, kshl(385),20,32,0, qks(385))
 DO  i = 1,30
   DO  j = 1,6
     ij =(i-1)*32 + j
     ji =(i-1)*6  + j
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
 CALL gmmats (kshl(1  ),30,6 ,0,q,6,6,1 ,qks(1  ))
 CALL gmmats (kshl(181),30,6 ,0,q,6,6,1 ,qks(181))
 CALL gmmats (kshl(361),30,20,0, qqqinv,20,18,0, qks(361))
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
   ee(i) = 0.0
 END DO
 ee( 1) =  ivect(1)
 ee( 2) =  jvect(1)
 ee( 3) =  kvect(1)
 ee( 6) =  ivect(2)
 ee( 7) =  jvect(2)
 ee( 8) =  kvect(2)
 ee(11) =  ivect(3)
 ee(12) =  jvect(3)
 ee(13) =  kvect(3)
 ee(19) =  ivect(1)
 ee(20) =  jvect(1)
 ee(24) =  ivect(2)
 ee(25) =  jvect(2)
 ee(29) =  ivect(3)
 ee(30) =  jvect(3)
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
     i2 =  5*(k-1) + i + 2
     k2 = 12+(k-1)*3 + i
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
       cm1(j,i1) = ctm(j,k1)
     END DO
   END DO
 END DO
 DO  k = 1,6
   DO  i = 1,3
     i2 =  5*(k-1) + i + 2
     k2 = 12+(k-1)*3 + i
     DO  j = 1,30
       cm1(j,i2) = ctm(j,k2)
     END DO
   END DO
 END DO
 DO  i = 1,1296
   cmt(i) = 0.0
 END DO
 
!     LUMPED MASS COMPUTATION
 
 IF (ipass /= 2) GO TO 490
 470 amass = (rhoy*vol + nsm*area)/6.
 DO  i = 1,1296,37
   cmt(i) = amass
 END DO
 ipass = 2
 GO TO 690
 
!     LOCATE THE TRANSFORMATION MATRICES FROM BASIC TO LOCAL (THAT IS
!     COORDINATE AT ANY GRID POINT IN WHICH DISPLACEMENT AND STRESSES
!     ARE R
!     - NOT NEEDED IF FIELD 7 IN GRID CARD IS ZERO)
 
!     TRANSFORM STIFFNESS MATRIX FROM ELEMENT COORDINATES TO BASIC
!     COORDINATES
 
!     TRANSFORM STIFFNESS MATRIX FROM BASIC COORDINAYES TO GLOBAL (DISP)
!     COORDINATES
 
!     INSERT THE 6X6 SUBMATRIX  INTO KGG MATRIX
 
 490 DO  i = 1,6
   SAVE(i) = nl(i)
 END DO
 DO  i = 1,6
   small(i) = i
   ismall = nl(i)
   DO  j = 1,6
     IF (ismall <= nl(j)) GO TO 510
     small(i) = j
     ismall = nl(j)
     510 CONTINUE
   END DO
   ism = small(i)
   nl(ism) = 1000000
 END DO
 DO  i = 1,6
   nl(i) = SAVE(i)
 END DO
 DO  i = 1,6
   sil1 = small(i)
   DO  j = i,6
     sil2 = small(j)
     DO  ii = 1,36
       balotr(ii) = 0.0
       ksup(ii)   = 0.0
     END DO
     DO  k = 1,5
       k1 = (sil1-1)*5 + k
       DO  l = 1,5
         l1 = (sil2-1)*5 + l
         csub(k,l) = cm1(k1,l1)
       END DO
     END DO
     CALL gmmats (ee,6,5,0, csub,5,5,0, csubt)
     CALL gmmats (csubt,6,5,0, ee,6,5,+1, ksupt)
     DO  k = 1,6
       DO  l = 1,6
         k1 =(k-1)*6 + l
         l1 =(l-1)*6 + k
         ksup(l1) = ksupt(k1)
       END DO
     END DO
     
!     TRANSFORM THE KSUP(36) FROM BASIC TO DISPLACEMENT COORDINATES
     
     IF (nl(sil1) == 0 .OR. ics(sil1) == 0) GO TO 610
     jj = 4*sil1 + 24
     CALL transs (iest(jj),trand)
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
     CALL gmmats (balotr(1),6,6,1, ksup(1),6,6,0, ksupt)
     DO  k = 1,36
       ksup(k) = ksupt(k)
     END DO
     610 CONTINUE
     IF (nl(sil2) == 0 .OR. ics(sil2) == 0) GO TO 650
     IF (j == i) GO TO 630
     jj = 4*sil2 + 24
     CALL transs (iest(jj),trand)
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
     630 CONTINUE
     CALL gmmats (ksup(1),6,6,0, balotr(1),6,6,0, ksupt)
     DO  k = 1,36
       ksup(k) = ksupt(k)
     END DO
     650 CONTINUE
     DO  ii = 1,6
       DO  jj = 1,6
         i1 = (i-1)*6 + ii
         j1 = (j-1)*6 + jj
         ctm(i1,j1) = ksub(jj,ii)
         ctm(j1,i1) = ksub(jj,ii)
       END DO
     END DO
   END DO
 END DO
 690 CALL emgout (cmt(1),cmt(1),1296,1,dict,ipass,iprec)
 IF (.NOT.imass .OR. ipass >= 2) RETURN
 
!     TO TO 295 TO COMPUTE LUMPED MASS MATRIX
!     GO TO 211 TO COMPUTE CONSIST. MASS MATRIX (THIS PATH DOES NOT
!     WROK)
 
 ipass = 3
 SELECT CASE ( ipass )
   CASE (    1)
     GO TO 720
   CASE (    2)
     GO TO 90
   CASE (    3)
     GO TO 470
 END SELECT
 
!     ERROR
 
 700 CONTINUE
 nogo  = .true.
 knogo = 1
 WRITE  (ioutpt,710) ufm,iest(1)
 710 FORMAT (a23,' 2416, MATRIX RELATING GENERALIZED PARAMETERS AND ',  &
     'GRID POINT DISPLACEMENTS IS SINGULAR.', //26X,  &
     'CHECK COORDINATES OF ELEMENT  TRSHL WITH ID',i9,1H.)
 720 CONTINUE
 RETURN
END SUBROUTINE ktshls
