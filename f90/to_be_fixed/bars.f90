SUBROUTINE bars
     
!     SINGLE PRECISION VERSION
 
!     THIS SUBROUTINE PROCESSES BAR  ELEMENT DATA TO PRODUCE STIFFNESS
!     AND MASS MATRICES. IF THE HEAT TRANSFER OPTION IS ON, CONDUCTIVITY
!     AND CAPACITY MATRICES ARE PRODUCED.
 
!     THIS ROUTINE WILL PRODUCE MASS MATRICES BY EITHER THE CONSISTENT
!     OR CONVENTIONAL MASS METHODS.
!     THE ECPT/EST ENTRIES FOR THE BAR (ELEMENT TYPE 34) ARE
 
!     ECPT( 1)  -  IELID          ELEMENT ID. NUMBER
!     ECPT( 2)  -  ISILNO(2)      * SCALAR INDEX NOS. OF THE GRID POINTS
!     ECPT( 3)  -    ...          *
!     ECPT( 4)  -  SMALLV(3)      $ REFERENCE VECTOR
!     ECPT( 5)  -    ...          $
!     ECPT( 6)  -    ...          $
!     ECPT( 7)  -  ICSSV          COOR. SYS. ID FOR SMALLV VECTOR
!     ECPT( 8)  -  IPINFL(2)      * PIN FLAGS
!     ECPT( 9)  -    ...          *
!     ECPT(10)  -  ZA(3)          $ OFFSET VECTOR FOR POINT A
!     ECPT(11)  -    ...          $
!     ECPT(12)  -    ...          $
!     ECPT(13)  -  ZB(3)          * OFFSET VECTOR FOR POINT B
!     ECPT(14)  -    ...          *
!     ECPT(15)  -    ...          *
!     ECPT(16)  -  IMATID         MATERIAL ID.
!     ECPT(17)  -  A              CROSS-SECTIONAL AREA
!     ECPT(18)  -  I1             $ AREA MOMENTS OF INERTIA
!     ECPT(19)  -  I2             $
!     ECPT(20)  -  FJ             TORSIONAL CONSTANT
!     ECPT(21)  -  NSM            NON-STRUCTURAL MASS
!     ECPT(22)  -  FE             FORCE ELEM DESCRIPTIONS (FORCE METHOD)
!     ECPT(23)  -  C1             * STRESS RECOVERY COEFFICIENTS
!     ECPT(24)  -  C2             *
!     ECPT(25)  -  D1             *
!     ECPT(26)  -  D2             *
!     ECPT(27)  -  F1             *
!     ECPT(28)  -  F2             *
!     ECPT(29)  -  G1             *
!     ECPT(30)  -  G2             *
!     ECPT(31)  -  K1             $ AREA FACTORS FOR SHEAR
!     ECPT(32)  -  K2             $
!     ECPT(33)  -  I12            AREA MOMENT OF INERTIA
!     ECPT(34)  -  MCSIDA         COOR. SYS. ID. FOR GRID POINT A
!     ECPT(35)  -  GPA(3)         * BASIC COORDINATES FOR GRID POINT A
!     ECPT(36)  -    ...          *
!     ECPT(37)  -    ...          *
!     ECPT(38)  -  MCSIDB         COOR. SYS. ID. FOR GRID POINT B
!     ECPT(39)  -  GPB(3)         $ BASIC COORDINATES FOR GRID POINT B
!     ECPT(40)  -    ...          $
!     ECPT(41)  -    ...          $
!     ECPT(42)  -  ELTEMP         AVG. ELEMENT TEMPERATURE
 
 LOGICAL :: basic,offset,nogo,aofset,bofset
 INTEGER :: dict(7),is12or(4),is21or(4),gsube,estid,iecpt(38)
 REAL :: k1,k2,i1,i2,i12,nsm,ke,kk,kep,m,mep,me,lr1,lr2,lb, l2b3,l2b6,limit
 DIMENSION       veci(3),vecj(3),veck(3),ecpt(42),ipin(10),ikk(4)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ksystm(100)
 COMMON /emgest/ ielid,isilno(2),smallv(3),icssv,ipinfl(2),za(3),  &
     zb(3),imatid,a,i1,i2,fj,nsm,fe,c1,c2,d1,d2,f1,f2,  &
     g1,g2,k1,k2,i12,mcsida,gpa(3),mcsidb,gpb(3),tempel
 COMMON /emgprm/ ixtra,jcore,ncore,dum(12),istif,imass,idamp,  &
     iprec,nogo,heat,icmbar,lcstm,lmat,lhmat
 COMMON /emgdic/ idumm, ldict,ngrids,elid,estid
 COMMON /emgtrx/ ke(144),kep(144),m(12,12),me(144),mep(144),  &
     kk(144),smalvn(6),ta(18),tb(9),vec(10), dela(6),delb(6)
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e,g,nu,rho,alpha,tsubo,gsube,sigt,sigc,sigs
 COMMON /hmtout/ fk
 EQUIVALENCE     (ksystm(2),ioutpt), (ksystm(56),iheat),  &
     (ecpt(1),iecpt(1),ielid), (ksystm(87),ksy87),  &
     (vec(1),veci(1)), (vec(4),vecj(1)), (vec(7),veck(1))
 DATA    ikk   / 1,7,73,79   /, epsi,epsi2 / 1.0E-18,1.0E-7 /
 DATA    is12or/ 1,37,109,73 /, is21or     / 73,109,37,1    /
 
 
 dict(1) = estid
 
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
 limit  = IABS(ksy87)*.01
 
!     NORMALIZE THE REFERENCE VECTOR WHICH LIES IN THE FIRST PRINCIPAL
!     AXIS PLANE  (FMMS - 36 P. 4)
 
 fl = SQRT(smallv(1)**2 + smallv(2)**2 + smallv(3)**2)
 IF (ABS(fl) < epsi) GO TO 7770
 DO  i = 1,3
   smalvn(i) = smallv(i)/fl
 END DO
 
!     DETERMINE IF POINT A AND B ARE IN BASIC COORDINATES OR NOT.
!     COMPUTE THE TRANSFORMATION MATRICES TA AND TB IF NECESSARY
 
 IF (icsida /= 0) CALL transs (ecpt(jcsida),ta)
 IF (icsidb /= 0) CALL transs (ecpt(jcsidb),tb)
 
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
 
 100 DO  i = 1,3
   jta = i + jcsida
   jtb = i + jcsidb
   veci(i) = ecpt(jta) - ecpt(jtb)
 END DO
 
!     SAVE IN A2B THE LENGTH OF BAR WITHOUT OFFSET, FROM GRID PT. A TO B
 
 a2b = SQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 
!     TRANSFORM THE OFFSET VECTORS IF NECESSARY
 
 IF (.NOT.aofset .AND. .NOT.bofset) GO TO 150
 
!     TRANSFORM THE OFFSET VECTOR FOR POINT A IF NECESSARY.
 
 idela = 1
 j = jofsta - 1
 DO  i = 1,3
   j = j + 1
   dela(i) = ecpt(j)
 END DO
 IF (icsida == 0) GO TO 120
 idela = 4
 CALL gmmats (ta,3,3,0,dela(1),3,1,0,dela(4))
 
!     TRANSFORM THE OFFSET VECTOR FOR POINT B IF NECESSARY
 
 120 idelb = 1
 j = jofstb - 1
 DO  i = 1,3
   j = j + 1
   delb(i) = ecpt(j)
 END DO
 IF (icsidb == 0) GO TO 140
 idelb = 4
 CALL gmmats (tb,3,3,0,delb(1),3,1,0, delb(4))
 
!     SINCE THERE WAS AT LEAST ONE NON-ZERO OFFSET VECTOR RECOMPUTE VECI
 
 140 DO  i = 1,3
   jta = i - 1 + idela
   jtb = i - 1 + idelb
   veci(i) = veci(i)+dela(jta) - delb(jtb)
 END DO
 
!     COMPUTE THE LENGTH OF THE BIG V (VECI) VECTOR AND NORMALIZE
 
 150 fl = 0.
 DO  i = 1,3
   veci(i) = - veci(i)
   fl = fl + veci(i)**2
 END DO
 fl = SQRT(fl)
 IF (ABS(fl) < epsi) GO TO 7770
 DO  i = 1,3
   veci(i) = veci(i)/fl
 END DO
 
!     NOW THAT LENGTH HAS BEEN COMPUTED, CHECK POSSIBLE OFFSET ERROR
!     ISSUE WARNING MESSAGE IF OFFSET EXCEEDS A2B BY 'LIMIT' PERCENT.
!     (DEFAULT IS 15 PERCENT, KSYSTM(87) WORD)
 
 IF (ABS(fl-a2b)/a2b <= limit) GO TO 170
 WRITE  (ioutpt,165) uwm,ielid
 165 FORMAT (a25,' - UNUSUALLY LARGE OFFSET IS DETECTED FOR CBAR ',  &
     'ELEMENT ID =',i8,' ***')
 IF (ksy87 <= 0) GO TO 170
 WRITE  (ioutpt,167) ksy87
 167 FORMAT (/5X,'(OFFSET BAR LENGTH EXCEEDS NON-OFFSET LENGTH BY',  &
     i4,' PERCENT, SET BY SYSTEM 87TH WORD)')
 ksy87 = -ksy87
 
!     BRANCH IF THIS IS A -HEAT- FORMULATION.
 
 170 IF (iheat == 1) GO TO 500
 
!     COMPUTE THE  SMALV0  VECTOR
 
 isv = 1
 IF (icssv == 0) GO TO 180
 isv = 4
 CALL gmmats (ta,3,3,0,smalvn(1),3,1,0,smalvn(4))
 
!     COMPUTE THE K VECTOR, VECK = VECI  X  SMALV0, AND NORMALIZE
 
 180 veck(1) = veci(2)*smalvn(isv+2) - veci(3)*smalvn(isv+1)
 veck(2) = veci(3)*smalvn(isv  ) - veci(1)*smalvn(isv+2)
 veck(3) = veci(1)*smalvn(isv+1) - veci(2)*smalvn(isv  )
 fll= SQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 IF (ABS(fll) < epsi2) GO TO 7770
 DO  i = 1,3
   veck(i) = veck(i)/fll
 END DO
 
!     COMPUTE THE J VECTOR, VECJ = VECK  X  VECI, AND NORMALIZE
 
 vecj(1) = veck(2)*veci(3) - veck(3)*veci(2)
 vecj(2) = veck(3)*veci(1) - veck(1)*veci(3)
 vecj(3) = veck(1)*veci(2) - veck(2)*veci(1)
 fll = SQRT(vecj(1)**2 + vecj(2)**2 + vecj(3)**2)
 IF (ABS(fll) < epsi2) GO TO 7770
 vecj(1) = vecj(1)/fll
 vecj(2) = vecj(2)/fll
 vecj(3) = vecj(3)/fll
 
!     SEARCH THE MATERIAL PROPERTIES TABLE FOR E,G AND THE DAMPING
!     CONSTANT.
 
 matidc = imatid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 
 IF (istif == 0) GO TO 600
 
!     IF ELASTICITY AND SHEAR MODULES BOTH ZERO, SKIP STIFFNESS
!     CALCULATION
 
 IF (e == 0. .AND. g == 0.) GO TO 600
 
!     SET UP INTERMEDIATE VARIABLES FOR ELEMENT STIFFNESS MATRIX
!     CALCULATION
 
 ASSIGN 305 TO k OR m
 205 bl    = fl
 blsq  = fl**2
 blcube= blsq*bl
 
!     COMPUTE SOME TERMS TO BE USED IN STIFFNESS MATRIX KE
 
 ei1  = e*i1
 ei2  = e*i2
 IF (k1 == 0.0 .OR. i12 /= 0.0) GO TO 210
 gak1 = g*a*k1
 r1   = (12.*ei1*gak1)/(gak1*blcube + 12.*bl*ei1)
 GO TO 220
 210 r1   = 12.*ei1/blcube
 220 IF (k2 == 0.0 .OR. i12 /= 0.0) GO TO 230
 gak2 = g*a*k2
 r2   = (12.*ei2*gak2)/(gak2*blcube + 12.*bl*ei2)
 GO TO 240
 230 r2   = 12.*ei2/blcube
 
 240 sk1 = .25*r1*blsq + ei1/bl
 sk2 = .25*r2*blsq + ei2/bl
 sk3 = .25*r1*blsq - ei1/bl
 sk4 = .25*r2*blsq - ei2/bl
 
 ael =  a*e/bl
 lr1 =  bl*r1/2.
 lr2 =  bl*r2/2.
 gjl =  g*fj/bl
 
!     CONSTRUCT  THE GENERAL 12X12 MATRIX FOR THE BAR ELEMENT
 
!                      **       **
!                      * K    K  *
!                      *  AA   AB*
!                K =   *  T      *
!                      * K    K  *
!                      *  AB   BB*
!                      **       **
 
 
!     FIRST SET THE COMPONENT CODE AND THE DOF
 
 icode = 63
 ndof  = 12
 nsq   = ndof**2
 
!     CONSTRUCT THE 12 X 12 MATRIX KE
 
 DO  i = 1,144
   ke(i)   = 0.
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
 IF (i12 == 0.) GO TO 303
 beta    =  12.*e*i12/blcube
 lb      =  bl*beta/2.0
 l2b3    =  blsq*beta/3.0
 l2b6    =  blsq*beta/6.0
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
 303 GO TO k OR m, (305,640,465,750)
 
!     DETERMINE IF THERE ARE NON-ZERO PIN FLAGS.
 
 305 ka = iecpt(jpina)
 kb = iecpt(jpinb)
 IF (ka == 0 .AND. kb == 0) GO TO 345
 
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
   IF (ke(ii) /= 0.) GO TO 320
   il = ipin(i)
   ii = ii - il
   DO  j = 1,12
     ii = ii + 1
     ke(ii) = 0.
     ke(il) = 0.
     il = il + 12
   END DO
   CYCLE
   320 DO  j = 1,12
     ji = 12*(j-1) + ipin(i)
     ij = 12*(ipin(i)-1) + j
     DO  ll = 1,12
       jll = 12*(j-1) + ll
       ill = 12*(ipin(i)-1) + ll
       kep(jll) = ke(jll) - (ke(ill)/ke(ii))*ke(ji)
     END DO
     kep(ij) = 0.
     kep(ji) = 0.
   END DO
   DO  k = 1,144
     ke(k) = kep(k)
   END DO
 END DO
 
!     DIVIDE KE INTO FOUR SUBMATRICES AND STORE IN OPEN CORE
 
!      E                   E                   E
!     K   = KK(1 TO 36)   K   = KK(37 TO 72)  K   = KK(73 TO 108)
!      AA                  AB                  BA
 
!      E
!     K   =  KK(109 TO 144)
!      BB
 
 
 345 j = 0
 DO  i = 1,72,12
   low = i
   lim = i + 5
   DO  k = low,lim
     j = j + 1
     kk(j    ) = ke(k   )
     kk(j+ 36) = ke(k+ 6)
     kk(j+ 72) = ke(k+72)
     kk(j+108) = ke(k+78)
   END DO
 END DO
 
 ASSIGN 465 TO k OR m
 
!     ZERO OUT THE ARRAY WHERE THE 3 X 3 MATRIX H AND THE W AND W  6 X 6
!     MATRICES WILL RESIDE.                                A      B
!      T
!     A   MATRIX NOW STORED IN KE
 
 358 DO  i = 1,9
   ke(i) = vec(i)
 END DO
 
 DO  i = 28,144
   ke(i) = 0.
 END DO
 
!     SET UP POINTERS
 
 basic  = icsida == 0
 jcsid  = jcsida
 offset = aofset
 jofset = jofsta
 DO  i = 1,2
   iwbeg = i*36
   
!     SET UP THE -G- MATRIX.  IG POINTS TO THE BEGINNING OF THE G MATRIX
!     G = AT X TI
   
   ig = 1
   IF (basic) GO TO 380
   CALL transs (ecpt(jcsid),ke(10))
   CALL gmmats (ke(1),3,3,0, ke(10),3,3,0, ke(19) )
   ig = 19
   
!     IF THERE IS A NON-ZERO OFFSET FOR THE POINT, SET UP THE D 3 X 3
!     MATRIX.
   
   380 IF (.NOT.offset) GO TO 385
   ke(10) =  0.
   ke(11) =  ecpt(jofset+2)
   ke(12) = -ecpt(jofset+1)
   ke(13) = -ke(11)
   ke(14) =  0.
   ke(15) =  ecpt(jofset)
   ke(16) = -ke(12)
   ke(17) = -ke(15)
   ke(18) =  0.
   
!     FORM THE 3 X 3 PRODUCT H = G X D, I.E., KE(28) = KE(IG) X KE(10)
   
   CALL gmmats (ke(ig),3,3,0, ke(10),3,3,0,ke(28))
   
   
!     FORM THE W SUBMATRICES IN KE(37) AND KE(73)
   
   
   385 ke(iwbeg+ 1) = ke(ig  )
   ke(iwbeg+ 2) = ke(ig+1)
   ke(iwbeg+ 3) = ke(ig+2)
   ke(iwbeg+ 7) = ke(ig+3)
   ke(iwbeg+ 8) = ke(ig+4)
   ke(iwbeg+ 9) = ke(ig+5)
   ke(iwbeg+13) = ke(ig+6)
   ke(iwbeg+14) = ke(ig+7)
   ke(iwbeg+15) = ke(ig+8)
   ke(iwbeg+22) = ke(ig  )
   ke(iwbeg+23) = ke(ig+1)
   ke(iwbeg+24) = ke(ig+2)
   ke(iwbeg+28) = ke(ig+3)
   ke(iwbeg+29) = ke(ig+4)
   ke(iwbeg+30) = ke(ig+5)
   ke(iwbeg+34) = ke(ig+6)
   ke(iwbeg+35) = ke(ig+7)
   ke(iwbeg+36) = ke(ig+8)
   IF (.NOT.offset) GO TO 390
   ke(iwbeg+ 4) = ke(28)
   ke(iwbeg+ 5) = ke(29)
   ke(iwbeg+ 6) = ke(30)
   ke(iwbeg+10) = ke(31)
   ke(iwbeg+11) = ke(32)
   ke(iwbeg+12) = ke(33)
   ke(iwbeg+16) = ke(34)
   ke(iwbeg+17) = ke(35)
   ke(iwbeg+18) = ke(36)
   390 basic  = icsidb == 0
   jcsid  = jcsidb
   offset = bofset
   jofset = jofstb
 END DO
 
!     CONVERT THE K PARTITIONS TO GLOBAL COORDINATES AND STORE IN KEP
 
 iaft = 37
 DO  i = 1,4
   ikx  = (i-1)*36 + 1
   ik   = ikx
   IF (i >= 3) ikx = (7-i-1)*36  + 1
   ifore = ((i-1)/2)*36 + 37
   CALL gmmats (ke(ifore),6,6,1, kk(ikx),6,6,0, ke(109))
   CALL gmmats (ke(109), 6,6,0, ke(iaft),6,6,0, kep(ik))
   iaft = 73
   IF (i == 3) iaft = 37
 END DO
 
!     REFORM THE K MATRIX (12X12) FROM THE FOUR SUBMATRICES (6X6) AND
!     ORDER  THE SUBMATRICES BY INCREASING SIL VALUE
 
 DO  ii = 1,4
   ix1 = ikk(ii)
   ix2 = ix1 + 60
   is  = is12or(ii)
   IF (isilno(1) > isilno(2)) is = is21or(ii)
   DO  i = ix1,ix2,12
     ip5 = i + 5
     DO  j = i,ip5
       ke(j) = kep(is)
       is = is + 1
     END DO
   END DO
 END DO
 
 GO TO k OR m, (305,640,465,750)
 
!     OUTPUT THE STIFFNESS MATRIX
 
 465 dict(2) = 1
 dict(3) = ndof
 dict(4) = icode
 dict(5) = gsube
 CALL emgout (ke(1),ke(1),nsq,1,dict,1,iprec)
 GO TO 600
 
!     THE MASS MATRIX IS GENERATED HERE.  IF THE PARAMETER ICMBAR IS
!     .LT. 0, CALL THE CONVENTIONAL MASS MATRIX GENERATION ROUTINE FOR
!     THE BAR.  OTHERWISE CALL THE ROUTINE TO GENERATE CONSISTENT MASS
!     MATRICES FOR THE BAR.
 
 600 const = (fl*(rho*a + nsm))/420.
 IF (imass == 0 .OR. const == 0.) RETURN
 IF (icmbar < 0) GO TO 800
 
!     CALCULATE THE CONSISTENT/CONVENTIONAL MASS MATRIX
 
!     CALL THE MAT ROUTINE TO FETCH SINGLE PRECISION MATERIAL PROPERTIES
 
 matidc = imatid
 matflg = 1
 eltemp = tempel
 CALL mat (iecpt(1))
 
 
!     COMPUTE TERMS OF THE ELEMENT MASS MATRIX
 
 bl22  = 22.*fl
 bl13  = 13.*fl
 blsq4 = 4.0*fl**2
 blsq3 = 3.0*fl**2
 
!     CONSTRUCT THE ELEMENT MASS MATRIX.
 
 DO  i = 1,12
   DO  j = 1,12
     m( i, j) = 0.
   END DO
 END DO
 m( 1, 1) = 175.
 m( 1, 7) = 35.
 m( 2, 2) = 156.
 m( 2, 6) = bl22
 m( 2, 8) = 54.
 m( 2,12) =-bl13
 m( 3, 3) = 156.
 m( 3, 5) =-bl22
 m( 3, 9) = 54.
 m( 3,11) = bl13
 m( 5, 5) = blsq4
 m( 5, 9) =-bl13
 m( 5,11) =-blsq3
 m( 6, 6) = blsq4
 m( 6, 8) = bl13
 m( 6,12) =-blsq3
 m( 7, 7) = 175.
 m( 8, 8) = 156.
 m( 8,12) =-bl22
 m( 9, 9) = 156.
 m( 9,11) = bl22
 m(11,11) = blsq4
 m(12,12) = blsq4
 
!     STORE THE UPPER TRIANGULAR PART OF THE MATRIX IN THE LOWER PART.
 
 DO  i = 1,10
   low = i + 1
   DO  j = low,12
     m(j,i) = m(i,j)
   END DO
 END DO
 
!     MULTIPLY BY CONSTANT AND STORE ROW-WISE IN THE ARRAY ME
 
 k = 0
 DO  i = 1,12
   DO  j = 1,12
     k = k + 1
     me(k) = const*m(i,j)
   END DO
 END DO
 
!     IF THERE ARE NO PIN FLAGS THERE IS NO NEED TO CALCULATE THE
!     ELEMENT STIFFNESS MATRIX
 
 ka = iecpt(jpina)
 kb = iecpt(jpinb)
 IF (ka == 0 .AND. kb == 0) GO TO 705
 
!     COMPUTE THE STIFFNESS MATRIX KE
 
 ASSIGN 640 TO k OR m
 GO TO 205
 
!     RETURN HERE AFTER COMPUTING THE STIFFNESS MATRIX
 
 
!     SET UP TNHE IPIN ARRAY
 
 640 DO  i = 1,5
   ipin(i  ) = MOD(ka,10)
   ipin(i+5) = MOD(kb,10) + 6
   IF (ipin(i+5) == 6) ipin(i+5) = 0
   ka = ka/10
   kb = kb/10
 END DO
 
!     ALTER THE ELEMENT MASS MATRIX DUE TO PIN FLAGS.  NOTE THAT THE
!     FOLLOWING CODE IS CONGRUENT AS IT WERE TO THE CODE IN SUBROUTINE
!     DBEAM IN THE DSMG1 MODULE.
 
 DO  j = 1,10
   IF (ipin(j) == 0) CYCLE
   jj = 12*(ipin(j)-1) + ipin(j)
   IF (ke(jj) == 0.) GO TO 680
   DO  i = 1,12
     ji = 12*(ipin(j)-1) + i
     ij = 12*(i-1) + ipin(j)
     DO  l1 = 1,12
       il = 12*(i-1) + l1
       lj = 12*(l1-1) + ipin(j)
       mep(il) = me(il) - ke(lj)*me(ji)/ke(jj) - ke(ji)*me(lj)/ke(jj)  &
           + ke(lj)*ke(ji)*me(jj)/ke(jj)**2
     END DO
   END DO
   DO  k = 1,144
     me(k) = mep(k)
   END DO
   
!     ZERO OUT THE IPIN(J) TH ROW AND COLUMN OF ME
   
   680 j1 = jj - ipin(j)
   j2 = ipin(j)
   DO  k = 1,12
     j1 = j1 + 1
     me(j1) = 0.
     me(j2) = 0.
     j2 = j2 + 12
   END DO
 END DO
 
!            E                  E                    E
!     STORE M  AT KK(1 TO 36), M  AT KK (37 TO 72), M  AT KK(73 TO 108)
!            AA                 AB                   BA
 
!            E
!     AND   M  AT KK(109 TO 144)
!            BB
 
 705 j = 0
 DO  i = 1,72,12
   low = i
   lim = low + 5
   DO  k = low,lim
     j = j + 1
     kk(j) = me(k)
     kk(j+ 36) = me(k+ 6)
     kk(j+ 72) = me(k+72)
     kk(j+108) = me(k+78)
   END DO
 END DO
 
!     CALCULATE THE TRANSFORMATION VECTORS
 
 ASSIGN 750 TO k OR m
 GO TO 358
 
!     OUTPUT THE CONSISTENT MASS MATRIX
 
 750 dict(2) = 1
 dict(3) = ndof
 dict(4) = icode
 dict(5) = 0
 CALL emgout (ke(1),ke(1),144,1,dict,2,iprec)
 RETURN
 
!     CALCULATE THE CONVENTIONAL MASS MATRIX HERE
 
!     GET RHO FROM MPT BY CALLING MAT
 
 800 matidc = imatid
 matflg = 4
 eltemp = tempel
 CALL mat (ecpt(1))
 DO  i = 1,72
   mep(i) = 0.
 END DO
 fm = .5*fl*(rho*a + nsm)
 
!     DETERMINE IF THE GRID POINT IS ASSOCIATED WITH A NON-ZERO OFFSET.
 
 jofset = 9
 DO  ii = 1,2
   ix = (ii-1)*36
   j  = jofset
   DO  i = 1,3
     j = j + 1
     IF (ecpt(j) /= 0.) GO TO 820
   END DO
   GO TO 840
   
!     FORM UPPER RIGHT CORNER OF THE MATRIX
   
   820 mep(ix+ 1) =  1.
   mep(ix+ 8) =  1.
   mep(ix+15) =  1.
   mep(ix+ 5) =  ecpt(jofset+3)
   mep(ix+ 6) = -ecpt(jofset+2)
   mep(ix+12) =  ecpt(jofset+1)
   mep(ix+10) = -mep(ix+ 5)
   mep(ix+16) = -mep(ix+ 6)
   mep(ix+17) = -mep(ix+12)
   mep(ix+20) = -mep(ix+ 5)
   mep(ix+21) = -mep(ix+ 6)
   mep(ix+25) = -mep(ix+10)
   mep(ix+27) = -mep(ix+12)
   mep(ix+31) = -mep(ix+16)
   mep(ix+32) = -mep(ix+17)
   mep(ix+22) =  ecpt(jofset+3)**2 + ecpt(jofset+2)**2
   mep(ix+29) =  ecpt(jofset+3)**2 + ecpt(jofset+1)**2
   mep(ix+36) =  ecpt(jofset+2)**2 + ecpt(jofset+1)**2
   mep(ix+23) = -ecpt(jofset+1)*ecpt(jofset+2)
   mep(ix+24) = -ecpt(jofset+1)*ecpt(jofset+3)
   mep(ix+30) = -ecpt(jofset+2)*ecpt(jofset+3)
   mep(ix+28) =  mep(ix+23)
   mep(ix+34) =  mep(ix+24)
   mep(ix+35) =  mep(ix+30)
   
!     MULTIPLY M BY THE CONSTANT FL
   
   DO  i = 1,36
     is = ix + i
     mep(is) = mep(is)*fm
   END DO
   GO TO 850
   
!     HERE WE HAVE A ZERO OFFSET VECTOR
   
   840 mep(ix+ 1) = fm
   mep(ix+ 8) = fm
   mep(ix+15) = fm
   850 jofset = 12
 END DO
 
!     INSERT M  AND M  SUBMATRICES INTO M ACCORDING TO INCREASING SIL
!             A      B
 
 DO  i = 1,144
   me(i) = 0.
 END DO
 
 IF (isilno(1)-isilno(2) > 0) THEN
   GO TO   880
 END IF
 870 ix1 = 1
 ix2 = 37
 GO TO 890
 880 ix1 = 37
 ix2 = 1
 890 CONTINUE
 DO  jj = 1,36
   mm = MOD(jj,6)
   IF (mm == 0) mm = 6
   i = ((jj-1)/6)*12  + mm
   j = i + 78
   me(i) = mep(ix1)
   me(j) = mep(ix2)
   ix1 = ix1 + 1
   ix2 = ix2 + 1
 END DO
 
!     OUTPUT THE CONVENTIONAL MASS MATRIX
 
 dict(2) = 1
 dict(3) = ndof
 dict(4) = icode
 dict(5) = 0
 
 CALL emgout ( me,me,144,1,dict,2,iprec)
 
 RETURN
 
!     HEAT FORMULATION CONTINUES HERE.  GET MATERIAL PROPERTY -K- FROM
!     HMAT
 
 500 matflg  = 1
 matidc  = iecpt(16)
 eltemp  = ecpt(42)
 dict(2) = 1
 dict(3) = 2
 dict(4) = 1
 dict(5) = 0
 IF (istif == 0) GO TO 540
 CALL hmat (ielid)
 
 kk(1) = fk*ecpt(17)/fl
 IF (kk(1) == 0.) GO TO 520
 kk(2) =-kk(1)
 kk(3) = kk(2)
 kk(4) = kk(1)
 CALL emgout (kk(1),kk(1),4,1,dict,1,iprec)
 
 520 matflg = 4
 
!     ERROR IN NEXT CARD FOR HEAT FORMULATION. REMOVED BY
!     G.CHAN/UNISYS, 1984
!     ALSO, CHANGE  GO TO 520  TO 540, 11-TH CARD ABOVE, AND
!     CALL EMGOUT BELOW AND WRITE TO THE 3-RD FILE INSTEAD OF THE 2-ND.
 
 
 CALL hmat (ielid)
 kk(1) = (fk*ecpt(17))*fl/2.
 IF (kk(1) == 0.) RETURN
 kk(2) = kk(1)
 dict(2) = 2
 CALL emgout (kk(1),kk(1),2,1,dict,3,iprec)
 540 RETURN
 
!     ERROR RETURNS
 
 7770 CONTINUE
 WRITE  (ioutpt,7775) ufm,ielid
 7775 FORMAT (a23,' 3176, BAR ELEMENT ID',i9,  &
     ' HAS ILLEGAL GEOMETRY OR CONNECTIONS.')
 nogo = .true.
 RETURN
END SUBROUTINE bars
