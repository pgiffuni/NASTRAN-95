SUBROUTINE strp11
     
!     PHASE 1 STRESS DATA RECOVERY FOR CTRPLT1 - HIGHER ORDER PLATE
!     ELEMENT
 
!     OUTPUTS FROM THIS PHASE FOR USE IN PHASE II ARE THE FOLLOWING
 
!     1) ELEMENT ID              WORDS    1     STORAGE IN PH1OUT  1
!     2) SIX SILS                WORDS    6                      2-7
!     3) BENDING THICKNESSES     WORDS    3                      8-10
!     4) STRESS POINTS           WORDS    8                     11-18
!     5) 4 NOS. 6 5X6 S MATRICES WORDS    720                   19-738
!     6) 3X1 S SUB T MATRIX      WORDS    3                    739-741
 
!     ECPT ENTRIES
!     AS IN STIFFNESS ROUTINE KTRPL1
 
 LOGICAL :: nots
 REAL :: j11,j12,j22,nsm,ivect,jvect,kvect
 DOUBLE PRECISION :: determ
 DIMENSION        NAME(2),INDEX(20,3),ics(6),nl(6),q(6,6),ind(6,3),  &
     emod(9),xc(6),yc(6),zc(6),qqq(20,20),qqqinv(360),  &
     ts6(40),ts7(59),iest(42),ivect(3),jvect(3),  &
     kvect(3),e(18),v1(3),v2(3),v3(3),e1(18),  &
     ph1ben(9),ph1shr(6),ph2(18),ph3(12),ph4(90),  &
     trans(9),balotr(36),d(9),dph1(9),g(4),gph1(6), nph1ou(990)
 COMMON /sdr2x5/  est(100),ph1out(990),forvec(24),  &
     x,y,z,dista,distb,distc,a1,a2,a3,b1,b2,b3,  &
     qqqinv,ts6,ts7,ph2,ph3,ph4,q,e,e1,trans,balotr
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin /  matid,matflg,eltemp,pla34,sinth,costh
 COMMON /matout/  em(6),rhoy,alf(3),tref,gsube,sigty,sigcy,sigsy,  &
     rj11,rj12,rj22
 
 
!     EQUIVALENCE IECPT WITH ECPT IN COMMON BLOCK /SMA1ET/ SINCE ECPT IS
!     A MIXED INTEGER AND REAL ARRAY
 
 EQUIVALENCE     (a,dista), (b,distb), (c,distc),  &
     (v1(1),est(19)),(v2(1),est(23)),(v3(1),est(27)), (iest(1),est(1)),  &
     (d11,em(1)),(d12,em(2)), (d13,em(3)), (d22,em(4)),(d23,em(5)), (d33,em(6))
 EQUIVALENCE     (nph1ou(1),ph1out(1))
 EQUIVALENCE     (ph1out(401),INDEX(1,1),ind(1,1))
 EQUIVALENCE     (ph1out(1),qqq(1,1))
 DATA  degra  /  0.0174532925            /
 DATA  BLANK  ,  NAME / 4H    , 4HCTRP, 4HLT1   /
 
 nots  =.false.
 idele = iest(1)
 DO  i = 1,6
   nl(i) = iest(i+1)
 END DO
 thetam = est(8)
 matid1 = iest(9)
 tmem1  = (est(10)*12.0)**0.333333333333
 tmem3  = (est(11)*12.0)**0.333333333333
 tmem5  = (est(12)*12.0)**0.333333333333
 matid2 = iest(13)
 tshr1  = est(14)
 tshr3  = est(15)
 tshr5  = est(16)
 nsm    = est(17)
 j      = 0
 DO  i = 24,44,4
   j      = j + 1
   ics(j) = iest(i)
   xc(j)  = est(i+1)
   yc(j)  = est(i+2)
   zc(j)  = est(i+3)
 END DO
 
!     IF TMEM3 OR TMEM5 IS ZERO OR BLANK, THEY WILL BE SET EQUAL TO
!     TMEM1
!     SO ALSO FOR TEMP3 OR TEMP5
 
 IF (tmem3 == 0.0 .OR. tmem3 == BLANK) tmem3 = tmem1
 IF (tmem5 == 0.0 .OR. tmem5 == BLANK) tmem5 = tmem1
 IF (tshr3 == 0.0 .OR. tshr3 == BLANK) tshr3 = tshr1
 IF (tshr5 == 0.0 .OR. tshr5 == BLANK) tshr5 = tshr1
 IF (tshr1 == 0.0) nots = .true.
 eltemp = est(48)
 theta1 = thetam*degra
 sinth  = SIN(theta1)
 costh  = COS(theta1)
 IF (ABS(sinth) <= 1.0E-06) sinth = 0.0
 
!     EVALUATE MATERIAL PROPERTIES
 
 matflg = 2
 matid  = matid1
 CALL mat (idele)
 
 emod(1) = d11
 emod(2) = d12
 emod(3) = d13
 emod(4) = d12
 emod(5) = d22
 emod(6) = d23
 emod(7) = d13
 emod(8) = d23
 emod(9) = d33
 matid   = matid2
 matflg  = 3
 j11     = 0.0
 j12     = 0.0
 j22     = 0.0
 IF (nots) GO TO 146
 CALL mat (idele)
 146 CONTINUE
 
!     CALCULATIONS FOR THE TRIANGLE
 
 CALL trif (xc,yc,zc,ivect,jvect,kvect,a,b,c,iest(1),NAME)
 CALL af (f,1,a,b,c,a1,a2,a3,tmem1,tmem3,tmem5,1)
 CALL af (f,1,a,b,c,b1,b2,b3,tshr1,tshr3,tshr5,1)
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e( i) = 0.0
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
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 
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
   
!     IF NO TRANSVERSE SHEAR GO TO 113
   
   IF (nots) GO TO 1137
   x = xc(i)
   y = yc(i)
   CALL strpts (ts6,nots)
   DO  jj = 1,20
     qqq(i2,jj) = qqq(i2,jj) - ts6(20+jj)
     qqq(i3,jj) = qqq(i3,jj) + ts6(   jj)
   END DO
   1137 CONTINUE
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
 
!     FOURTH ARGUMENT IS A DUMMY LOCATION FOR INVERSE AND HENCE TS1(1)
!     IS U
 
 CALL invers (20,qqq,20,ts6(1),0,determ,ising,INDEX)
 
!     ISING EQUAL TO 2 IMPLIES THAT QQQ IS SINGULAR
 
!     FIRST 18 COLUMNS OF QQQ INVERSE IS THE QQQINV FOR USE IN STIFFNESS
!     MATRIX CALCULATIONS
 
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
   ph1out(i) = est(i)
 END DO
 ph1out( 8) = tmem1
 ph1out( 9) = tmem3
 ph1out(10) = tmem5
 ph1out(11) = est(18)
 ph1out(12) = est(19)
 ph1out(13) = est(20)
 ph1out(14) = est(21)
 ph1out(15) = est(22)
 ph1out(16) = est(23)
 DO  jj = 1,4
   jj1 = jj*2 - 1
   IF (jj /= 4) x = xc(jj1)
   IF (jj /= 4) y = yc(jj1)
   IF (jj == 4) x = (xc(1)+xc(3)+xc(5))/3.0
   IF (jj == 4) y = (yc(1)+yc(3)+yc(5))/3.0
   IF (jj == 4) ph1out(17) = (a1+a2*x+a3*y)/2.0
   IF( jj == 4) ph1out(18) = -ph1out(17)
   DO  i = 1,60
     ts7(i) = 0.0
   END DO
   ai = ph1out(7+jj)**3/12.0
   IF (jj == 4) ai = ph1out(17)**3/1.5
   DO  i = 1,9
     d(i) = emod(i)*ai
   END DO
   x2  = x*x
   xy  = x*y
   y2  = y*y
   x3  = x2*x
   x2y = x2*y
   xy2 = x*y2
   y3  = y2*y
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
   CALL gmmats (ts7,3,20,0, qqqinv,20,18,0, ph4(1))
   CALL strpts (ts6,nots)
   CALL gmmats (ts6,2,20,0, qqqinv,20,18,0, ph4(55))
   DO  ii = 1,6
     IF (ics(ii) == 0) GO TO 130
     j = 4*ii + 20
     CALL transs (iest(j),trans)
     DO  j = 1,3
       l = 6*(j-1) + 1
       m = 3*(j-1) + 1
       balotr(l   ) = trans(m  )
       balotr(l+ 1) = trans(m+1)
       balotr(l+ 2) = trans(m+2)
       balotr(l+21) = trans(m  )
       balotr(l+22) = trans(m+1)
       balotr(l+23) = trans(m+2)
     END DO
     CALL gmmats (e,6,3,+1, balotr,6,6,0, e1)
     GO TO 133
     130 CONTINUE
     DO  i = 1,3
       DO  j = 1,6
         i1 = (i-1)*6 + j
         j1 = (j-1)*3 + i
         e1(i1) = e(j1)
       END DO
     END DO
     133 CONTINUE
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
     IF (nots) GO TO 166
     thk  = b1 + b2*x + b3*y
     g(1) = em(6)*thk
     g(2) = 0.0
     g(3) = 0.0
     g(4) = g(1)
     CALL gmmats (g,2,2,0, ph1shr,2,3,0, gph1)
     GO TO 168
     166 CONTINUE
     gph1(1) = ph1shr(1)
     gph1(2) = ph1shr(2)
     gph1(3) = ph1shr(3)
     gph1(4) = ph1shr(4)
     gph1(5) = ph1shr(5)
     gph1(6) = ph1shr(6)
     168 CONTINUE
     CALL gmmats (gph1,2,3,0, e1,3,6,0, ph3)
     DO  i = 1,3
       DO  j = 1,6
         i1 = (i-1)*6 + j
         i2 = i1 + 18
         j1 = (ii-1)*30 + (jj-1)*180 + i1 + 18
         j2 = j1 + 18
         ph1out(j1) = ph2(i1)
         IF (i /= 3) ph1out(j2) = ph3(i1)
       END DO
     END DO
   END DO
   jj1 = (jj-1)*3 + 1
   CALL gmmats (d,3,3,0, alf,3,1,0, ph1out(738+jj1))
 END DO
 RETURN
END SUBROUTINE strp11
