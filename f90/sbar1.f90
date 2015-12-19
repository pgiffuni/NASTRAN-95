SUBROUTINE sbar1
     
!     THIS ROUTINE IS PHASE 1 OF STRESS DATA RECOVERY FOR THE BAR
!     ELEMENT. MUCH OF THE CODE WAS LIFTED FROM THE KBAR SUBROUTINE.
 
!     ECPT FOR THE BAR
 
!     ECPT( 1)  -  IELID       ELEMENT ID. NUMBER
!     ECPT( 2)  -  ISILNO(2)   * SCALAR INDEX NOS. OF THE GRID POINTS
!     ECPT( 3)  -    ...       *
!     ECPT( 4)  -  SMALLV(3)   $ REFERENCE VECTOR
!     ECPT( 5)  -    ...       $
!     ECPT( 6)  -    ...       $
!     ECPT( 7)  -  ICSSV       COOR. SYS. ID FOR SMALLV VECTOR
!     ECPT( 8)  -  IPINFL(2)   * PIN FLAGS
!     ECPT( 9)  -    ...       *
!     ECPT(10)  -  ZA(3)       $ OFFSET VECTOR FOR POINT A
!     ECPT(11)  -    ...       $
!     ECPT(12)  -    ...       $
!     ECPT(13)  -  ZB(3)       * OFFSET VECTOR FOR POINT B
!     ECPT(14)  -    ...       *
!     ECPT(15)  -    ...       *
!     ECPT(16)  -  IMATID      MATERIAL ID.
!     ECPT(17)  -  A           CROSS-SECTIONAL AREA
!     ECPT(18)  -  I1          $ AREA MOMENTS OF INERTIA
!     ECPT(19)  -  I2          $
!     ECPT(20)  -  FJ          POLAR MOMENT OF INERTIA
!     ECPT(21)  -  NSM         NON-STRUCTURAL MASS
!     ECPT(22)  -  FE          FORCE ELEMENT DESCRIPTIONS (FORCE METHOD)
!     ECPT(23)  -  C1          * STRESS RECOVERY COEFFICIENTS
!     ECPT(24)  -  C2          *
!     ECPT(25)  -  D1          *
!     ECPT(26)  -  D2          *
!     ECPT(27)  -  F1          *
!     ECPT(28)  -  F2          *
!     ECPT(29)  -  G1          *
!     ECPT(30)  -  G2          *
!     ECPT(31)  -  K1          $ AREA FACTORS FOR SHEAR
!     ECPT(32)  -  K2          $
!     ECPT(33)  -  I12         AREA MOMENT OF INERTIA
!     ECPT(34)  -  MCSIDA      COOR. SYS. ID. FOR GRID POINT A
!     ECPT(35)  -  GPA(3)      * BASIC COORDINATES FOR GRID POINT A
!     ECPT(36)  -    ...       *
!     ECPT(37)  -    ...       *
!     ECPT(38)  -  MCSIDB      COOR. SYS. ID. FOR GRID POINT B
!     ECPT(39)  -  GPB(3)      $ BASIC COORDINATES FOR GRID POINT B
!     ECPT(40)  -    ...       $
!     ECPT(41)  -    ...       $
!     ECPT(42)  -  ELTEMP      AVG. ELEMENT TEMPERATURE
 
 LOGICAL :: abasic,bbasic,basic,aofset,bofset,offset
 REAL :: l,lsq,lcube,i1,i2,k1,k2,ke,kep,i12,nsm,lr1,lr2,lb, l2b3,l2b6,hut(36)
 DIMENSION       veci(3),vecj(3),veck(3),ecpt(100),iecpt(100),  &
     ipin(10),ta(18),tb(9),smalv0(6)
 
!     SDR2 PHASE I INPUT AND OUTPUT COMMON BLOCK
 
 COMMON /sdr2x5/ ielid,isilno(2),smallv(3),icssv,ipinfl(2),za(3),  &
     zb(3),imatid,a,i1,i2,fj,nsm,fe,c1,c2,d1,d2,f1,f2,  &
     g1,g2,k1,k2,i12,mcsida,gpa(3),mcsidb,gpb(3), tempel,dum3(58)
 COMMON /sdr2x5/ jelid,jsilno(2),sa(36),sb(36),out(19),therm(30)
 
!     SDR2 SCRATCH BLOCK
 
 COMMON /sdr2x6/ ke(144),kep(144),dela(6),delb(6)
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e,g,nu,rho,alpha,t_sub_0,g_sub_e,sigt, sigc,sigs
 EQUIVALENCE     (ielid,ecpt(1),iecpt(1)),(ta(10),tb(1))
 
 
!     SET UP POINTERS TO COOR. SYS. IDS., OFFSET VECTORS, AND PIN FLAGS.
!     ICSIDA AND ICSIDB ARE COOR. SYS. IDS.
 
 jcsida = 34
 jcsidb = 38
 jofsta = 10
 jofstb = 13
 jpina  =  8
 jpinb  =  9
 icsida = iecpt(34)
 icsidb = iecpt(38)
 
!     NORMALIZE THE REFERENCE VECTOR WHICH LIES IN THE FIRST PRINCIPAL
!     AXIS PLANE  (FMMS - 36 P. 4)
 
 fl = 0.0
 DO  i = 1,3
   fl = fl + smallv(i)**2
 END DO
 fl = SQRT(fl)
 DO  i = 1,3
   smallv(i) = smallv(i)/fl
 END DO
 
!     DETERMINE IF POINT A AND B ARE IN BASIC COORDINATES OR NOT.
 
 abasic = .true.
 bbasic = .true.
 IF (icsida /= 0) abasic = .false.
 IF (icsidb /= 0) bbasic = .false.
 
!     COMPUTE THE TRANSFORMATION MATRICES TA AND TB IF NECESSARY
 
 IF (.NOT.abasic) CALL transs (ecpt(jcsida),ta)
 IF (.NOT.bbasic) CALL transs (ecpt(jcsidb),tb)
 
!     DETERMINE IF WE HAVE NON-ZERO OFFSET VECTORS.
 
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
 
!     FORM THE CENTER AXIS OF THE BEAM WITHOUT OFFSETS.
 
 100 veci(1) = ecpt(jcsida+1) - ecpt(jcsidb+1)
 veci(2) = ecpt(jcsida+2) - ecpt(jcsidb+2)
 veci(3) = ecpt(jcsida+3) - ecpt(jcsidb+3)
 
!     TRANSFORM THE OFFSET VECTORS IF NECESSARY
 
 IF (.NOT.aofset .AND. .NOT.bofset) GO TO 150
 
!     TRANSFORM THE OFFSET VECTOR FOR POINT A IF NECESSARY.
 
 idela = 1
 j = jofsta - 1
 DO  i = 1,3
   j = j + 1
   dela(i) = ecpt(j)
 END DO
 IF (abasic) GO TO 120
 idela = 4
 CALL gmmats (ta,3,3,0, dela(1),3,1,0, dela(4))
 
!     TRANSFORM THE OFFSET VECTOR FOR POINT B IF NECESSARY
 
 120 idelb = 1
 j = jofstb - 1
 DO  i = 1,3
   j = j + 1
   delb(i) = ecpt(j)
 END DO
 IF (bbasic) GO TO 140
 idelb = 4
 CALL gmmats (tb,3,3,0, delb(1),3,1,0, delb(4))
 
!     SINCE THERE WAS AT LEAST ONE NON-ZERO OFFSET VECTOR RECOMPUTE VECI
 
 140 veci(1) = veci(1) + dela(idela  ) - delb(idelb  )
 veci(2) = veci(2) + dela(idela+1) - delb(idelb+1)
 veci(3) = veci(3) + dela(idela+2) - delb(idelb+2)
 
!     COMPUTE THE LENGTH OF THE BIG V (VECI) VECTOR AND NORMALIZE
 
 150 veci(1) = -veci(1)
 veci(2) = -veci(2)
 veci(3) = -veci(3)
 fl = SQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 DO  i = 1,3
   veci(i) = veci(i)/fl
 END DO
 
!     COMPUTE THE SMALL V SUB 0 VECTOR, SMALV0.  ** CHECK THIS LOGIC **
 
 DO  i = 1,3
   smalv0(i) = smallv(i)
 END DO
 isv = 1
 IF (icssv == 0) GO TO 180
 isv = 4
 CALL gmmats (ta,3,3,0, smalv0(1),3,1,0, smalv0(4))
 
!     COMPUTE THE K VECTOR, VECK = VECI  X  SMALV0, AND NORMALIZE
 
 180 veck(1) = veci(2)*smalv0(isv+2) - veci(3)*smalv0(isv+1)
 veck(2) = veci(3)*smalv0(isv  ) - veci(1)*smalv0(isv+2)
 veck(3) = veci(1)*smalv0(isv+1) - veci(2)*smalv0(isv  )
 fll     = SQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 veck(1) = veck(1)/fll
 veck(2) = veck(2)/fll
 veck(3) = veck(3)/fll
 
!     COMPUTE THE J VECTOR, VECJ = VECK  X  VECI, AND NORMALIZE
 
 vecj(1) = veck(2)*veci(3) - veck(3)*veci(2)
 vecj(2) = veck(3)*veci(1) - veck(1)*veci(3)
 vecj(3) = veck(1)*veci(2) - veck(2)*veci(1)
 fll     = SQRT(vecj(1)**2 + vecj(2)**2 + vecj(3)**2)
 vecj(1) = vecj(1)/fll
 vecj(2) = vecj(2)/fll
 vecj(3) = vecj(3)/fll
 
!     CALL MAT TO GET MATERIAL PROPERTIES.
 
 matidc = imatid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 
!     SET UP INTERMEDIATE VARIABLES FOR ELEMENT STIFFNESS MATRIX
!     CALCULATION
 
 l    = fl
 lsq  = l**2
 lcube= lsq*l
 ei1  = e*i1
 ei2  = e*i2
 IF (k1 == 0.0 .OR. i12 /= 0.0) GO TO 210
 gak1 = g*a*k1
 r1   = (12.0*ei1*gak1)/(gak1*lcube + 12.0*l*ei1)
 GO TO 220
 210 r1   =  12.0*ei1/lcube
 220 IF (k2 == 0.0 .OR. i12 /= 0.0) GO TO 230
 gak2 = g*a*k2
 r2   = (12.0*ei2*gak2)/(gak2*lcube + 12.0*l*ei2)
 GO TO 240
 230 r2   = 12.0*ei2/lcube
 
!     COMPUTE THE -SMALL- K-S, SK1, SK2, SK3 AND SK4
 
 240 sk1 = 0.25*r1*lsq + ei1/l
 sk2 = 0.25*r2*lsq + ei2/l
 sk3 = 0.25*r1*lsq - ei1/l
 sk4 = 0.25*r2*lsq - ei2/l
 
!     COMPUTE THE TERMS THAT WILL BE NEEDED FOR THE 12 X 12 MATRIX KE
 
 ael = a*e /l
 lr1 = l*r1/2.0
 lr2 = l*r2/2.0
 gjl = g*fj/l
 
!     CONSTRUCT THE 12 X 12 MATRIX KE
 
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
 beta    = 12.0*e*i12/lcube
 lb      = l  *beta/2.0
 l2b3    = lsq*beta/3.0
 l2b6    = lsq*beta/6.0
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
 
!     DETERMINE IF THERE ARE NON-ZERO PIN FLAGS.
 
 255 ka = iecpt(jpina)
 kb = iecpt(jpinb)
 IF (ka == 0 .AND. kb == 0) GO TO 325
 
!     SET UP THE IPIN ARRAY
 
 DO  i = 1,5
   ipin(i  ) = MOD(ka,10)
   ipin(i+5) = MOD(kb,10) + 6
   IF (ipin(i+5) == 6) ipin(i+5) = 0
   ka = ka/10
   kb = kb/10
 END DO
 
!     ALTER KE MATRIX DUE TO PIN FLAGS.
 
 DO  i = 1,10
   IF (ipin(i) == 0) CYCLE
   ii = 13*ipin(i) - 12
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
     ji = 12*(j-1) + ipin(i)
     ij = 12*(ipin(i)-1) + j
     DO  ll = 1,12
       jll = 12*(j-1) + ll
       ill = 12*(ipin(i)-1) + ll
       kep(jll) = ke(jll) - (ke(ill)/ke(ii))*ke(ji)
     END DO
     kep(ij) = 0.0
     kep(ji) = 0.0
   END DO
   DO  k = 1,144
     ke(k) = kep(k)
   END DO
 END DO
 
!            E
!     STORE K   IN KEP(1),...,KEP(36) AND
!            AA
 
!            E
!     STORE K   IN KEP(37),...,KEP(72)
!            AB
 
 325 j   = 0
 DO  i = 1,72,12
   low = i
   lim = low + 5
   DO  k = low,lim
     j = j + 1
     kep(j   ) = ke(k  )
     kep(j+36) = ke(k+6)
   END DO
 END DO
 
!     COMPUTE THERMAL MATRIX
 
 DO  i = 1,30
   hut(i)  = 0.0
 END DO
 alphal  = alpha*l
 alpl6   = alphal*l/6.0
 alpl3   = alpl6*2.0
 alpl2   = alphal/2.0
 hut( 1) = alphal
 hut( 7) = alpl6
 hut( 8) = alpl3
 hut(14) = alpl6
 hut(15) = alpl3
 hut(24) = alpl2
 hut(25) = alpl2
 hut(27) =-alpl2
 hut(28) =-alpl2
 CALL gmmats (kep(1),6,6,0, hut,6,5,0, therm(1))
 
!                                                            T
!     STORE VECI, VECJ, VECK IN KE(1),...,KE(9) FORMING THE A  MATRIX.
 
 ke(1) = veci(1)
 ke(2) = veci(2)
 ke(3) = veci(3)
 ke(4) = vecj(1)
 ke(5) = vecj(2)
 ke(6) = vecj(3)
 ke(7) = veck(1)
 ke(8) = veck(2)
 ke(9) = veck(3)
 
!     SET POINTERS SO THAT WE WILL BE WORKING WITH POINT A.
 
 basic  = abasic
 jcsid  = jcsida
 offset = aofset
 jofset = jofsta
 iwbeg  = 0
 ikel   = 1
 iab    = 1
 INDEX  = isilno(1)
 
!     ZERO OUT THE ARRAY WHERE THE 3 X 3 MATRIX AND THE W  AND W  6 X 6
!     MATRICES WILL RESIDE.                              A      B
 
 DO  i = 28,108
   ke(i) = 0.0
 END DO
 
!     SET UP THE -G- MATRIX. IG POINTS TO THE BEGINNING OF THE G MATRIX.
!     G = AT X TI
 
 360 ig = 1
 IF (basic) GO TO 370
 CALL transs (ecpt(jcsid),ke(10))
 CALL gmmats (ke(1),3,3,0, ke(10),3,3,0, ke(19))
 ig = 19
 
!     IF THERE IS A NON-ZERO OFFSET FOR THE POINT, SET UP THE D 3 X 3
!     MATRIX.
 
 370 IF (.NOT.offset) GO TO 380
 ke(10) =  0.0
 ke(11) =  ecpt(jofset+2)
 ke(12) = -ecpt(jofset+1)
 ke(13) = -ke(11)
 ke(14) =  0.0
 ke(15) =  ecpt(jofset)
 ke(16) = -ke(12)
 ke(17) = -ke(15)
 ke(18) =  0.0
 
!     FORM THE 3 X 3 PRODUCT H = G X D, I.E., KE(28) = KE(IG) X KE(10)
 
 CALL gmmats (ke(ig),3,3,0, ke(10),3,3,0, ke(28))
 
 
!     FORM THE W  MATRIX OR THE W  MATRIX IN KE(37) OR KE(73) DEPENDING
!               A                B
!     UPON WHICH POINT - A OR B - IS UNDER CONSIDERATION.  G WILL BE
!     STORED IN THE UPPER LEFT AND LOWER RIGHT CORNERS.  H, IF NON-ZERO,
!     WILL BE STORED IN THE UPPER RIGHT CORNER.
 
 
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
 IF (.NOT.offset) GO TO 390
 ke(iwbeg+40) = ke(28)
 ke(iwbeg+41) = ke(29)
 ke(iwbeg+42) = ke(30)
 ke(iwbeg+46) = ke(31)
 ke(iwbeg+47) = ke(32)
 ke(iwbeg+48) = ke(33)
 ke(iwbeg+52) = ke(34)
 ke(iwbeg+53) = ke(35)
 ke(iwbeg+54) = ke(36)
 
!                             E                    E
!     FORM THE PRODUCT  S =  K   * W   OR  S   = K    * W  , DEPENDING
!                        A    AA    A       B     AB     B
!     UPON WHICH POINT WE ARE WORKING WITH.
 
 390 CALL gmmats (kep(ikel),6,6,0, ke(iwbeg+37),6,6,0, sa(iab))
 
!     IF THE POINT UNDER CONSIDERATION IS POINT B WE ARE FINISHED.  IF
!     NOT, SET UP POINTS AND INDICATORS FOR WORKING WITH POINT B.
 
 IF (iwbeg == 36) GO TO 500
 basic  = bbasic
 jcsid  = jcsidb
 offset = bofset
 jofset = jofstb
 iwbeg  = 36
 ikel   = 37
 iab    = 37
 INDEX  = isilno(2)
 DO  i = 28,36
   ke(i)  = 0.0
 END DO
 GO TO 360
 
!     FILL REMAINDER OF OUTPUT BLOCK.
 
 500 jelid     = ielid
 jsilno(1) = isilno(1)
 jsilno(2) = isilno(2)
 out( 1) = a*e*alpha
 out( 2) = a*e/l
 out( 3) = a
 out( 4) = fj
 out( 5) = i1
 out( 6) = i2
 out( 7) = i12
 out( 8) = c1
 out( 9) = c2
 out(10) = d1
 out(11) = d2
 out(12) = f1
 out(13) = f2
 out(14) = g1
 out(15) = g2
 out(16) = t_sub_0
 out(17) = sigt
 out(18) = sigc
 out(19) = l
 RETURN
END SUBROUTINE sbar1
