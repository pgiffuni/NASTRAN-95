SUBROUTINE dtrbsc (iopt,npivot)
     
!     IOPT = 1  IMPLIES THAT A CLOUGH TRIANGLE IS CALLING
!     IOPT = 2  IMPLIES THAT A QUADRILATERAL IS CALLING
 
!     ECPT LISTS OF NECESSARY VARIABLES
 
!     POSITION     TRIA1      QUAD1
!     ========     =====      =====
!     ECPT(51)      EID        EID
!     ECPT(52)      SIL1       SIL1
!     ECPT(53)      SIL2       SIL2
!     ECPT(54)      SIL3       SIL3
!     ECPT(55)      THETA      SIL4
!     ECPT(56)      MATID1     THETA
!     ECPT(57)      T1         MATID1
!     ECPT(58)      MATID2     T1
!     ECPT(59)      EYE        MATID2
!     ECPT(60)      MATID3     EYE
!     ECPT(61)      T2         MATID3
!     ECPT(62)      NSMASS     T2
!       :
!     ECT.
 
 
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN)                      :: npivot
 DOUBLE PRECISION :: g,g2x2,ar,eye,xbar,ybar,xcsq,ycsq,xbsq,xcyc,px2,  &
     py2,pxy2,xbar3,ybar3,ybar2,t2,r,sp,t,u,r2,s2,di, dumdp,c,a,d,s,j2x2,determ
 DOUBLE PRECISION :: xb2,xc2,yc2,xbc,sx,sy,sxy,xsubb,xsubc,ysubc
 DIMENSION        d(9),di(5,5),j2x2(4),s(18),necpt(51),a(144)
 COMMON /matin /  matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/  g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alph12,  &
     tsubd,gsube,sigten,sigcon,sigshe,g2x211,g2x212, g2x222
 COMMON /ds1adp/  g(9),g2x2(4),ar,eye,xbar,ybar,xcsq,ycsq,xbsq,xcyc  &
     ,                px2,py2,pxy2,xbar3,ybar3,ybar2,t2,r,sp,t,u,r2,s2,  &
     dumdp(81),c(24,3),sx,sy,sxy,xsubb,xsubc,ysubc
 COMMON /ds1aet/  ecpt(100)
 EQUIVALENCE      (a(1),d(1),g(1)),(necpt(1),ecpt(1)),  &
     (j2x2(1),dumdp(1)),(di(1,1),g(1))
 
!//////
!     CALL BUG (4HTBIG,30,SX,12)
!//////
!     IF NO TRANSVERSE SHEAR FLEXIBILITY EXISTS THE H-INVERSE IS
!     CALCULATED DIRECTLY.  TEST AS FOLLOWS
 
 IF (ecpt(iopt+60) /= 0.0 .AND. necpt(iopt+59) /= 0) GO TO 30
 
!     THE H-INVERSE MATRIX IS GENERATED IN TWO PARTITIONS
!         HB IS IN POSITIONS C(7,2) TO C(24,2)
!         HC IS IN POSITIONS C(7,3) TO C(24,3)
 
 10 nohyq = 1
 r  = 1.0/xsubb
 sp = 1.0/ysubc
 t  = sp*xsubc
 u  = r*r*sp*t
 r2 = r*r
 s2 = sp**2
 
 DO  i = 1,72
   c(i,1) = 0.0D0
 END DO
 
 c(7 ,2) = 3.0D0*r2
 c(9 ,2) = r
 c(11,2) = r
 c(13,2) =-c(7,2)*t**2
 c(14,2) =-r*t
 c(15,2) = c(14,2)*t
 c(16,2) =-2.0D0*r2*r
 c(18,2) =-r2
 c(19,2) =-6.0D0*r*u*(xsubb-xsubc)
 c(20,2) =-r*sp
 c(21,2) = u*(3.0D0*xsubc -2.0D0*xsubb)
 c(22,2) = r*t*u*(6.0D0*xsubb - 4.0D0*xsubc)
 c(23,2) = r*sp*t
 c(24,2) = 2.0D0*t*u*(xsubb - xsubc)
 
 c(13,3) = 3.0D0*s2
 c(14,3) =-sp
 c(15,3) = sp*t
 c(21,3) =-s2
 c(22,3) =-2.0D0*s2*sp
 c(23,3) = s2
 GO TO 110
 
!     THE  MATERIAL COEFFICIENTS FOR TRANSVERSE SHEAR ARE CALCULATE HERE
!     AND THE H-INVERSE MATRIX IS GENERATED THE NORMAL WAY
 
!     GET THE G2X2 MATRIX
 
 30 matid  = necpt(iopt+59)
 inflag = 3
 CALL mat (ecpt(51))
 IF (g2x211 == 0. .AND. g2x212 == 0. .AND. g2x222 == 0.) GO TO 10
 t2      = ecpt(iopt+60)
 g2x2(1) = g2x211*t2
 g2x2(2) = g2x212*t2
 g2x2(3) = g2x212*t2
 g2x2(4) = g2x222*t2
 
 determ  = g2x2(1)*g2x2(4) - g2x2(3)*g2x2(2)
 j2x2(1) = g2x2(4)/determ
 j2x2(2) =-g2x2(2)/determ
 j2x2(3) =-g2x2(3)/determ
 j2x2(4) = g2x2(1)/determ
 
!     SETTING UP G MATRIX
 
 inflag = 2
 matid  = necpt(iopt+57)
 CALL mat (necpt(51))
 
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
 
 eye = ecpt(iopt+58)
 DO  i = 1,9
   d(i) = g(i)*eye
 END DO
 
!     (H  ) IS PARTITIONED INTO A LEFT AND RIGHT PORTION AND ONLY THE
!       YQ  RIGHT PORTION IS COMPUTED AND USED AS A (2X3). THE LEFT
!           2X3 PORTION IS NULL.  THE RIGHT PORTION WILL BE STORED AT
!           A(73) THRU A(78) UNTIL NOT NEEDED ANY FURTHER.
 
 
 
 temp  = 2.0D0*d(2) + 4.0D0*d(9)
 a(73) =-6.0D0*(j2x2(1)*d(1) + j2x2(2)*d(3))
 a(74) =-j2x2(1)*temp - 6.0D0*j2x2(2)*d(6)
 a(75) =-6.0D0*(j2x2(1)*d(6) + j2x2(2)*d(5))
 a(76) =-6.0D0*(j2x2(2)*d(1) + j2x2(4)*d(3))
 a(77) =-j2x2(2)*temp - 6.0D0*j2x2(4)*d(6)
 a(78) =-6.0D0*(j2x2(2)*d(6) + j2x2(4)*d(5))
 
!     THE ABOVE 6 ELEMENTS NOW REPRESENT THE (H  ) MATRIX (2X3)
!                                              YQ
 
 xbar = (xsubb + xsubc)/3.0D0
 ybar = ysubc/3.0D0
 
 xcsq  = xsubc**2
 ycsq  = ysubc**2
 xbsq  = xsubb**2
 xcyc  = xsubc*ysubc
 px2   = (xbsq + xsubb*xsubc + xcsq)/6.0D0
 py2   = ycsq/6.0D0
 pxy2  = ysubc*(xsubb + 2.0D0*xsubc)/12.0D0
 xbar3 = 3.0D0*xbar
 ybar3 = 3.0D0*ybar
 ybar2 = 2.0D0*ybar
 
!     F1LL (HBAR) MATRIX STORING AT A(37) THRU A(72)
 
 DO  i = 37,72
   a(i)  = 0.0D0
 END DO
 
 a(37) = xbsq
 a(40) = xbsq*xsubb
 a(44) = xsubb
 a(49) =-2.0D0*xsubb
 a(52) =-3.0D0*xbsq
 a(55) = xcsq
 a(56) = xcyc
 a(57) = ycsq
 a(58) = xcsq*xsubc
 a(59) = ycsq*xsubc
 a(60) = ycsq*ysubc
 a(62) = xsubc
 a(63) = ysubc*2.0D0
 a(65) = xcyc *2.0D0
 a(66) = ycsq *3.0D0
 a(67) =-2.0D0*xsubc
 a(68) =-ysubc
 a(70) =-3.0D0*xcsq
 a(71) =-ycsq
 
!     ADD TO 6 OF THE (HBAR) ELEMENTS THE RESULT OF (H  )(H  )
!                                                     UY   YQ
!     THE PRODUCT IS FORMED DIRECTLY IN THE ADDITION PROCESS BELOW.
!     NO (H  ) MATRIX IS ACTUALLY COMPUTED DIRECTLY.
!          UY
 
!     THE FOLLOWING IS THEN PER STEPS 6 AND 7 PAGE -16- MS-17.
 
 DO  i = 1,3
   a(i+39) = a(i+39) + xsubb*a(i+72)
   a(i+57) = a(i+57) + xsubc*a(i+72) + ysubc*a(i+75)
 END DO
 
!     AT THIS POINT INVERT  (H) WHICH IS STORED AT A(37) THRU A(72)
!     STORE INVERSE BACK IN A(37) THRU A(72)
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL inverd (6,a(37),6,a(73),0,determ,ising,a(79))
 
!     CHECK TO SEE IF H WAS SINGULAR
 
!     ISING = 2 IMPLIES SINGULAR MATRIX THUS ERROR CONDITION.
 IF(ising /= 2) GO TO 90
 CALL mesage (-30,33,ecpt(1))
 RETURN
 
!     PARTITION H-INVERSE AND STORE IN C2 AND C3 LOCATIONS 7 THRU 24
 
 90 DO  i = 1,6
   ih = 6*i -6
   ic = 3*i -3
   
   DO  j = 1,3
     jh= ih + j + 36
     jc= ic + j + 6
     c(jc,2) = a(jh)
     c(jc,3) = a(jh+3)
   END DO
 END DO
 nohyq = 0
 
!     THIS ENDS ADDED COMPUTATION FOR CASE OF T2 NOT ZERO
 
!     THE C1, C2, AND C3 MATRICES ARE GENERATED WITH THE FOLLOWING CODE
!     FIRST GENERATE THE S MATRICES IN POSITIONS 1 THRU 9 AND 10 THRU 18
 
 110 DO  i = 1,18
   s(i) = 0.0
 END DO
 DO  i = 1,9,4
   s(i  ) = 1.0
   s(i+9) = 1.0
 END DO
 s(  3) =-xsubb
 s( 11) = ysubc
 s( 12) =-xsubc
 
!     COMPUTE HA  AND STORE IN  CA, POSITIONS 7 THRU 24
 
!         HA =  -(HB TIMES SB + HC TIMES SC)
 
 CALL gmmatd (c(7,2),6,3,0, s(1),3,3,0, a(37))
 CALL gmmatd (c(7,3),6,3,0, s(10),3,3,0, a(55))
 
 DO  i = 1,18
   
   c(i+6,1) = -a(i+36) - a(i+54)
 END DO
 
!     COMPUTE  HYQ TIMES HX  AND STORE IN CX POSITIONS 1 THRU 6
!     (THE FIRST THREE COLUMNS OF HYQ ARE NULL)
 
 IF (nohyq == 1) GO TO 160
 
 DO  i = 1,3
   CALL gmmatd (a(73),2,3,0, c(16,i),3,3,0, c(1,i))
 END DO
 
 160 c(3,1) = c(3,1) - 1.0D0
 c(5,1) = c(5,1) + 1.0D0
 
!     THE INTEGRALS FOR THE  KDQQ MATRIX ARE GENERATED HERE
 
 yc2 = ysubc**2
 xb2 = xsubb**2
 xc2 = xsubc**2
 xbc = xsubb*xsubc
 
 di(1,1) = 1.0D0
 di(1,2) = ysubc/3.0D0
 di(1,3) = yc2/6.0D0
 di(1,4) = yc2*ysubc/10.0D0
 di(1,5) = yc2**2/15.0D0
 di(2,1) = (xsubb + xsubc)/3.0D0
 di(2,2) = ysubc*(xsubb + 2.0D0*xsubc)/12.0D0
 di(2,3) = di(1,3)*(xsubb + 3.0D0*xsubc)/5.0D0
 di(2,4) = di(1,4)*(xsubb + 4.0D0*xsubc)/6.0D0
 di(3,1) = (xb2 +xbc + xc2)/6.0D0
 di(3,2) = di(1,2)*(xb2 + 2.0D0*xbc + 3.0D0*xc2)/10.0D0
 di(3,3) = di(1,3)*(xb2 + 3.0D0*xbc + 6.0D0*xc2)/15.0D0
 di(4,1) = (xsubb + xsubc)*(xb2 + xc2)/10.0D0
 di(4,2) = di(1,2)*((xsubb + 2.0D0*xsubc)*xb2 +  &
     (3.0D0*xsubb + 4.0D0*xsubc)*xc2)/20.0D0
 di(5,1) = (xb2*xb2 + xb2*xbc + xbc*xbc + xbc*xc2 + xc2*xc2)/15.0
 
 ar = xsubb*ysubc*DBLE(ecpt(iopt+56))/2.0D0
 DO  i = 1,5
   ic = 6 - i
   DO  j = 1,ic
     di(i,j) = di(i,j)*ar
   END DO
 END DO
 
!     THE ABOVE INTEGRALS  D(I,J) CORRESPOND TO THE DOCUMENTED
!     VALUES  I(I-1,J-1).  ZERO INDICES DONT ALWAYS COMPILE.
 
!     THE DIFFERENTIAL STIFFNESS MATRIX IN GENERALIZED COORDINATES IS
!     CREATED BELOW AT POSITIONS A(28) TO A(91)
 
 a(28) = sx*di(1,1)
 a(29) = sxy*di(1,1)
 a(30) = 2.0D0*sx*di(2,1)
 a(31) = sx*di(1,2) + sxy*di(2,1)
 a(32) = 2.0D0*sxy*di(1,2)
 a(33) = 3.0D0*sx *di(3,1)
 a(34) = sx*di(1,3) + 2.0*sxy*di(2,2)
 a(35) = 3.0D0*sxy*di(1,3)
 
 a(37) = sy*di(1,1)
 a(38) = 2.0D0*sxy*di(2,1)
 a(39) = sxy*di(1,2) + sy*di(2,1)
 a(40) = 2.0D0*sy*di(1,2)
 a(41) = 3.0D0*sxy*di(3,1)
 a(42) = sxy*di(1,3) + 2.0D0*sy*di(2,2)
 a(43) = 3.0D0*sy*di(1,3)
 
 a(46) = 4.0D0*sx*di(3,1)
 a(47) = 2.0D0*(sx*di(2,2) + sxy*di(3,1))
 a(48) = 4.0D0*sxy*di(2,2)
 a(49) = 6.0D0*sx*di(4,1)
 a(50) = 2.0D0*(sx*di(2,3) + 2.0D0*sxy*di(3,2))
 a(51) = 6.0D0*sxy*di(2,3)
 
 a(55) = sx*di(1,3) + 2.0D0*sxy*di(2,2)+sy*di(3,1)
 a(56) = 2.0D0*(sxy*di(1,3) + sy*di(2,2))
 a(57) = 3.0D0*(sx* di(3,2) + sxy*di(4,1))
 a(58) = sx*di(1,4) + 3.0D0*sxy*di(2,3) + 2.0D0*sy*di(3,2)
 a(59) = 3.0D0*(sxy*di(1,4) + sy*di(2,3))
 
 a(64) = 4.0D0*sy*di(1,3)
 a(65) = 6.0D0*sxy*di(3,2)
 a(66) = 2.0D0*(sxy*di(1,4) + 2.0D0*sy*di(2,3))
 a(67) = 6.0D0*sy*di(1,4)
 
 a(73) = 9.0D0*sx*di(5,1)
 a(74) = 3.0D0*(sx*di(3,3) + 2.0D0*sxy*di(4,2))
 a(75) = 9.0D0*sxy*di(3,3)
 
 a(82) = sx*di(1,5) + 4.0D0*sxy*di(2,4) + 4.0D0*sy*di(3,3)
 a(83) = 3.0D0*sxy*di(1,5) + 6.0D0*sy*di(2,4)
 
 a(91) = 9.0D0*sy*di(1,5)
 
!     FILL IN SYMMETRIC TERMS
 
 DO  i = 2,8
   ih = i - 1
   DO  j = 1,ih
     ic =  8*(i-1) + j
     jc =  8*(j-1) + i
     a(ic+27) = a(jc+27)
   END DO
 END DO
 
!     AT THIS STAGE THE 3X3 MATRIX PARTITIONS MAY BE GENERATED
!     THE ACTUAL MATRICES DEPEND ON IOPT
 
 ic = npivot
 IF (ic == 0) GO TO 200
 CALL gmmatd (c(1,ic),8,3,1, a(28),8,8,0, a(92))
 DO  i = 1,3
   ih=  9*(i-1) + 1
   CALL gmmatd (a(92),3,8,0, c(1,i),8,3,0, a(ih))
 END DO
!//////
!     CALL BUG (4HTBKD,300,A,54)
!//////
 
!     AT THIS STAGE THE QUADRILATERAL CALCULATIONS ARE COMPLETE
 
 200 IF (iopt == 2) RETURN
 
!     THE TRIANGLE SUBROUTINE  MUST RETURN THE FOLLOWING DATA
!         KAC,KBC,KCC  IN POSITIONS  A(28) THRU A(54) -I=NPIVOT
!             S        IN POSITIONS  A(55) THRU A(72)
!           H-INVERSE  IN POSITIONS  A(73) THRU A(108)
 
 CALL gmmatd (a(28),8,8,0, c(1,3),8,3,0, a(92))
 DO  i = 1,3
   ih = 28 + 9*(i-1)
   CALL gmmatd (c(1,i),8,3,1, a(92),8,3,0, a(ih))
 END DO
 
!     RECALCULATE THE S MATRIX (IT WAS DESTROYED) -
!     PLACE IN A(55 THRU 72)
 
 DO  i = 1,18
   a(i+54) = 0.0
 END DO
 DO  i = 1,9,4
   a(i+54) = 1.0
   a(i+63) = 1.0
 END DO
 a(57) =-xsubb
 a(65) = ysubc
 a(66) =-xsubc
 
!     EXTRACT THE H-INVERSE MATRIX FROM THE C MATRICES
!     STORE AT POSITIONS A(73) THRU A(108)
 
 DO  i = 1,6
   ih = 6*i - 6
   ic = 3*i - 3
   
   DO  j = 1,3
     jh = ih + j + 72
     jc = ic + j + 6
     a(jh  ) = c(jc,2)
     a(jh+3) = c(jc,3)
   END DO
 END DO
 RETURN
END SUBROUTINE dtrbsc
