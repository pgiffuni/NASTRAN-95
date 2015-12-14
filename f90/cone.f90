SUBROUTINE cone (ti,z)
     
!     THIS ROUTINE COMPUTES THE THERMAL LOADS ON A AXISYMMETRIC CONE
 
 
 REAL, INTENT(IN)                         :: ti(2)
 REAL, INTENT(OUT)                        :: z(1)
 REAL :: i00      ,i10
 REAL :: i01      ,i11      ,i21      ,i31      ,i41
 REAL :: i02      ,i12      ,i22
 REAL :: i03      ,i13      ,i23      ,i33
 REAL :: pa(8)    ,xi(6)    ,ecpt(35)
 REAL :: n2d33    ,nsp      ,ncp      ,nspopi
 REAL :: eht(96)  ,huq(100) ,hyq(10)
 COMMON /condas/   pi       ,twopi    ,radeg    ,degra    , s4pisq
 COMMON /trimex/   mecpt(35)
 COMMON /matin /   matid    ,inflag   ,temp     ,stress   ,sinth  , costh
 COMMON /matout/   g11,g12  ,g13,g22  ,g23,g33  ,rho,alph1,alph2  ,  &
     alph3    ,tsub0    ,gsube    ,sigten   ,sigcom ,  &
     sigshe   ,g2x211   ,g2x212   ,g2x222
 EQUIVALENCE      (ecpt( 1) ,mecpt(1)),(ecpt( 9), ts    )
 EQUIVALENCE      (ecpt(28) ,ra      ),(ecpt(32), rb    )
 EQUIVALENCE      (ecpt(29) ,za      ),(ecpt(33), zb    )
 EQUIVALENCE      (ecpt( 6) ,matid2  ),(ecpt( 8), matid3)
 EQUIVALENCE      (gshear   ,g12     )
 DATA    one   /  1.0   /
 
!     DEFINITION OF VARIABLES
 
!        ECPT  ENTRIES FOR CONE
 
!     ECPT(1)  INTEGER   ELEMENT ID = 1000*ELID + HARMONIC
!     ECPT(2)  INTEGER   SIL  A
!     ECPT(3)  INTEGER   SIL  B
!     ECPT(4)  INTEGER   MAT ID  1
!     ECPT(5)  REAL      T    MEMBRANE THICKNESS
!     ECPT(6)  INTEGER   MAT ID  2
!     ECPT(7)  REAL      MOMENT OF INERTIA
!     ECPT(8)  INTEGER   MAT ID 3
!     ECPT(9)  REAL      SHEAR THICKNESS
!     ECPT(10) REAL      NON -STRUCTRAL  MASS
!     ECPT(11) REAL      Z1
!     ECPT(12) REAL      Z2
!     ECPT(13) REAL      PHI 1
!     ECPT(14) REAL          2
!     ECPT(15) REAL          3
!     ECPT(16) REAL          4
!     ECPT(17) REAL          5
!     ECPT(18) REAL          6
!     ECPT(19) REAL          7
!     ECPT(20) REAL          8
!     ECPT(21) REAL          9
!     ECPT(22) REAL         10
!     ECPT(23) REAL         11
!     ECPT(24) REAL         12
!     ECPT(25) REAL         13
!     ECPT(26) REAL         14
!     ECPT(27) INTEGER   COORDINANT SYSTEM FOR POINT  A
!     ECPT(28) REAL      R   (A)
!     ECPT(29) REAL      Z   (A)
!     ECPT(30) REAL      NULL
!     ECPT(31) INTEGER   COORDINANT SYSTEM FOR POINT B
!     ECPT(32) REAL      R    (B)
!     ECPT(33) REAL      Z    (B)
!     ECPT(34) REAL      NULL
!     ECPT(35) REAL      TEMPERATURE OF MATERIAL
 
!     XL       LENGTH  BETWEEN  POINTS
!     SP       SINE  OF  PHI
!     CP       COSINE OF PHI
!     I-S      INTEGRAL  FROM  PAGE 46 MS,28
!     MATID    MATERIAL ID  (MAT 1 CARD)
!     INFLAG   OPTION  2  OF MAT ROUTINE
!     TEMP     MATERIAL TEMPERATURE
!     SINTH    0.0  DUMMY
!     COSTH    1.0  DUMMY
!     XN       HARMONIC NUMBER
!     PA(8)    TOTAL LOAD VECTOR
!     XI(6)    CYLINDRICAL LOAD
 
 
!     IF MEMBRANE THICKNESS = 0, THEN LOAD IS ZERO
 
 IF (ecpt(5) == 0.0) GO TO 160
 
!     COMPUTE  L, SINPHI, COSPHI
 
 rbma = rb - ra
 zbma = zb - za
 xl2  = rbma**2 + zbma**2
 xl   = SQRT(xl2)
 IF (xl == 0.0) GO TO 160
 sp = rbma/xl
 cp = zbma/xl
 
!     COMPUTE  I-S
 
 xl4 = xl2*xl2
 rav = (ra + rb)*0.5
 i00 = xl *rav
 i10 = xl2*(ra + 2.0*rb)/6.0
 i01 = xl
 i11 = xl2/2.0
 i21 = xl2*xl/3.0
 i31 = xl4/4.0
 i41 = xl4*xl/5.0
 
!     SET UP FOR MAT ROUTINE
 
 matid = mecpt(4)
 inflag= 2
 temp  = ecpt(35)
 sinth = 0.0
 costh = 1.0
 CALL mat (mecpt(1))
 
!     COMPUTE COEFICCIENTS
 
 f  = (g12*alph2 + g22*alph1)*ecpt(5)*pi
 ff = (g11*alph2 + g12*alph1)*ecpt(5)*pi
 
!     COMPUTE  A
 
 a = (ti(1)-tsub0)*f
 
!     COMPUTE  B
 
 b = (ti(2)-ti(1))/xl*f
 
!     COMPUTE  C
 
 c = (ti(1)-tsub0)*ff
 
!     COMPUTE  D
 
 d = (ti(2)-ti(1))/xl*ff
 
!     DECODE  N
 
 ixn = mecpt(1)/1000
 xn  = mecpt(1) - ixn*1000 - 1
 
!     COMPUTE  PA
 
 f  = i01*a + i11*b
 ff = i11*a + i21*b
 pa(1) = xn*f
 pa(2) = xn*ff
 pa(3) = sp*f
 pa(4) = sp*ff + i00*c + i10*d
 pa(5) = cp*f
 pa(6) = cp*ff
 pa(7) = cp*(i21*a + i31*b)
 pa(8) = cp*(i31*a + i41*b)
 
!     CHECK HARMONIC NO.  IF(XN = 0.0) DOUBLE PA VECTOR
 
 IF (xn /= 0.0) GO TO 30
 DO  i = 1,8
   pa(i) = 2.0*pa(i)
 END DO
 
!     OMPUTE TRANSFORMATION MATRIX HUQ. SEE MS-28, PP. 15, 16, 24, 25
 
 30 DO  i = 1,100
   huq(i) = 0.0
 END DO
 huq(  1) = one
 huq( 13) = one
 huq( 25) = one
 huq( 36) = one
 huq( 41) = cp/ra
 huq( 45) = xn/ra
 huq( 49) = one
 huq( 51) = one
 huq( 52) = xl
 huq( 63) = one
 huq( 64) = xl
 huq( 75) = one
 huq( 76) = xl
 huq( 77) = xl2
 huq( 78) = huq(77)*xl
 huq( 86) = one
 huq( 87) = 2.0*xl
 huq( 88) = 3.0*huq(77)
 huq( 91) = cp/rb
 huq( 92) = huq(91)*xl
 huq( 95) = xn/rb
 huq( 96) = huq(95)*xl
 huq( 97) = huq(95)*xl2
 huq( 98) = huq(96)*xl2
 huq( 99) = one
 huq(100) = xl
 
!     CHCEK IF HYQ VECTOR NEEDED
 
 IF (matid2 == 0   .OR. matid3 == 0  ) GO TO 60
 IF (ecpt(7) == 0.0 .OR. ecpt(9) == 0.0) GO TO 60
 
!     FORM  (D) = I*(G)
 
 d11 = ecpt(7)*g11
 d12 = ecpt(7)*g12
 d22 = ecpt(7)*g22
 d33 = ecpt(7)*g33
 
!     PICK UP GSHEAR FROM MAT
 
 inflag = 1
 matid  = matid3
 temp   = ecpt(35)
 CALL mat (mecpt(1))
 IF (gshear == 0.0) GO TO 60
 
!     COMPUTE INTEGRALS
 
 b  = sp
 b2 = b*b
 b3 = b*b2
 b4 = b*b3
 rlog = ALOG(rb/ra)
 rasq = ra*ra
 rbma2  = rbma*rav
 orbora = one/rb - one/ra
 twora  = ra + ra
 
!     IF SP = 0 EVALUATE INTEGRALS DIFFERENTLY
 
 IF (sp /= 0.0) GO TO 45
 temp1= rav*rav
 temp3= xl2*xl
 i02  = xl/rav
 i12  = xl2/(2.0*rav)
 i22  = temp3/(3.0*rav)
 i03  = xl/temp1
 i13  = xl2/(2.0*temp1)
 i23  = temp3/(3.0*temp1)
 i33  = (xl2*xl2)/(4.0*temp1)
 GO TO 49
 45 CONTINUE
 i02 = rlog/b
 i12 = (rbma - ra*rlog)/b2
 i22 = (rbma2 - twora*rbma + rasq*rlog)/b3
 i03 =-orbora/b
 i13 = (rlog + ra*orbora)/b2
 i23 = (rbma - twora*rlog - rasq*orbora)/b3
 i33 = (rbma2 - 3.0*ra*rbma + 3.0*rasq*rlog + rasq*ra*orbora)/b4
 
!     COMPUTE HYQ
 
 49 CONTINUE
 cp2 = cp*cp
 sp2 = sp*sp
 xn2 = xn*xn
 opi = one/pi
 n2d33  = xn2*d33
 sp2d22 = sp2*d22
 oq  = xl*ts*gshear*rav + i02*(n2d33 + sp2d22)*opi
 oq  = one/oq
 nsp = xn*sp
 ncp = xn*cp
 nspopi  = nsp*opi
 twod33  = 2.0*d33
 temp1   = d12*orbora
 temp2   = nspopi*(d22 + d33)
 temp3   = xn*nspopi*(twod33 + d22)
 temp4   = oq*0.5*n2d33*cp*opi
 temp5   = opi*(xn2*twod33 + sp2d22)
 temp6   = d12*xn2*xl2/rb
 temp7   = nspopi*cp*0.5
 hyq( 1) = oq*(temp1*ncp - temp7*i03*(d33 + 2.0*d22))
 hyq( 2) = oq*(ncp*xl/rb*d12 - temp7*i13*(3.0*d33 + d22)  &
     + 1.5*ncp*opi*i02*d33)
 hyq( 3) = temp4*i03
 hyq( 4) = temp4*i13
 hyq( 5) = oq*(temp1*xn2  -  temp3*i03)
 hyq( 6) = oq*(d12*xn2*xl/rb - temp3*i13 + temp5*i02)
 hyq( 7) = oq*(2.0*d11*(ra-rb) + temp6 + 2.0*i12*temp5 - temp3*i23)
 hyq( 8) = oq*(-d11*6.*xl*rb + temp6*xl + 3.*i22*temp5 - temp3*i33)
 hyq( 9) =-oq*temp2*i02
 hyq(10) = oq*(xn*xl*(d12 + d33) - temp2*i12)
 DO  i = 1,10
   huq(i+30) = huq(i+30) - hyq(i)
   huq(i+80) = huq(i+80) - hyq(i)
 END DO
 
 itest = 1
 GO TO 61
 60 itest = 0
 huq(41) = 0.0
 huq(45) = 0.0
 huq(91) = 0.0
 huq(92) = 0.0
 huq(95) = 0.0
 huq(96) = 0.0
 huq(97) = 0.0
 huq(98) = 0.0
 huq(99) = 0.0
 61 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (10,huq(1),10,dum,0,determ,ising,eht(1))
 IF (ising == 2) CALL mesage (-30,40,mecpt(1))
 IF (itest /= 0) GO TO 62
 huq( 85) = 0.0
 huq(100) = 0.0
 62 CONTINUE
 
!     COMPLETE SOLUTION
 
!     FIRST OBTAIN PRODUCTS
!                       T
!        EHAT  =  (E)(H  )      AND STORE AT EHT(1) . . . EHT(48)
!                      A
 
!                       T
!        EHBT  =  (E)(H  )      AND STORE AT EHT(49). . . EHT(96)
!                      B
!                                /
!              WHERE  (HUQ) = (HA/HB)
!                                /
!              AND
!                             0    CP   SP   0    0
 
!                             1    0    0    0    0
 
!                             0    CP  -SP   0    0
!                  E MATRIX =
!                             0    0    0    0    SP
 
!                             0    0    0    1    0
 
!                             0    0    0    0    CP
 
 inc1 = 0
 inc2 = 0
 110 DO  i = 1,8
   krow = i + inc1
   ncol = (i-1)*10 + inc2
   eht(krow   ) = sp*huq(ncol+2) + cp*huq(ncol+3)
   eht(krow+ 8) =    huq(ncol+1)
   eht(krow+16) = cp*huq(ncol+2) - sp*huq(ncol+3)
   eht(krow+24) = sp*huq(ncol+5)
   eht(krow+32) =    huq(ncol+4)
   eht(krow+40) = cp*huq(ncol+5)
 END DO
 IF (inc1 > 0) GO TO 130
 inc1 = 48
 inc2 = 5
 GO TO 110
 
!     PERFORM TRANSFORMATION OF LOAD VECTOR
 
 130 DO  j = 1,2
   CALL gmmats (eht(48*j-47),6,8,0,pa(1),8,1,0,xi(1))
   k = mecpt(j+1) - 1
   DO  i = 1,6
     k = k + 1
     z(k) = z(k) + xi(i)
   END DO
 END DO
 
 160 RETURN
 
END SUBROUTINE cone
