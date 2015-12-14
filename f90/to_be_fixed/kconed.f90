SUBROUTINE kconed
     
!     DOUBLE PRECISION CONEAX ROUTINE
 
!     FOUR KCONE VERSIONS
!     KCONES  FOR MACHINES WITH 60 OR 64 BIT WORD (e.g. CDC, CRAY).
!             S.P. COMPUTATION IS USED
!     KCONE2, SIMILAR TO KCONES, EXECPT CERTAIN CRITICAL AREAS ARE
!             COMPUTED IN D.P. FOR IMPROVED ACCURACY
!     KCONED  FOR MAHCINES WITH LESS THEN 60 BIT WORD, WITHOUT QUAD
!             PRECISION SOFTWARE SUPPORT (e.g. DEC3100).
!             D.P. COMPUTAION IS USED
!     KCONEQ, SIMILAR TO KCONED, EXECPT CERTAIN CRITICAL AREAS ARE
!             COMPUTED IN QUAD PREC. FOR IMPROVED ACCURACY
 
!     ORIGINALLY, THIS ROUTINE CALLS KCONEX AND KCONEY/Z. THESE THREE
!     SUPPORTING ROUTINES ARE NOW MOVED INTO KCONED (AND ALSO KCONES)
 
!     ECPT( 1) = ELEMENT ID             INTEGER        ECT
!     ECPT( 2) = SIL PT A               INTEGER        ECT
!     ECPT( 3) = SIL PT B               INTEGER        ECT
!     ECPT( 4) = MATID 1                INTEGER        EPT
!     ECPT( 5) = T   (MEMBRANE THICK)   REAL           EPT
!     ECPT( 6) = MATID 2                INTEGER        EPT
!     ECPT( 7) = I   (MOM.OF INERTIA)   REAL           EPT
!     ECPT( 8) = MATID 3                INTEGER        EPT
!     ECPT( 9) = TS  (SHEAR THICKNESS)  REAL           EPT
!     ECPT(10) = NON-STRUCTURAL-MASS    REAL           EPT
!     ECPT(11) = Z1                     REAL           EPT
!     ECPT(12) = Z2                     REAL           EPT
!     ECPT(13) = PHI  1                 REAL           EPT
!     ECPT(14) = PHI  2                 REAL           EPT
!     ECPT(15) = PHI  3                 REAL           EPT
!     ECPT(16) = PHI  4                 REAL           EPT
!     ECPT(17) = PHI  5                 REAL           EPT
!     ECPT(18) = PHI  6                 REAL           EPT
!     ECPT(19) = PHI  7                 REAL           EPT
!     ECPT(20) = PHI  8                 REAL           EPT
!     ECPT(21) = PHI  9                 REAL           EPT
!     ECPT(22) = PHI 10                 REAL           EPT
!     ECPT(23) = PHI 11                 REAL           EPT
!     ECPT(24) = PHI 12                 REAL           EPT
!     ECPT(25) = PHI 13                 REAL           EPT
!     ECPT(26) = PHI 14                 REAL           EPT
!     ECPT(27) = COORD. SYS. ID PT.1    INTEGER        BGPDT
!     ECPT(28) = RADIUS PT. 1           REAL           BGPDT
!     ECPT(29) = DISTANCE TO PT.1       REAL           BGPDT
!     ECPT(30) = NULL                   REAL           BGPDT
!     ECPT(31) = COORD. SYS. ID PT.2    INTEGER        BGPDT
!     ECPT(32) = RADIUS PT 2            REAL           BGPDT
!     ECPT(33) = DISTANCE TO PT. 2      REAL           BGPDT
!     ECPT(34) = NULL                   REAL           BGPDT
!     ECPT(35) = ELEMENT TEMPERATURE    REAL           GEOM3
 
 
 INTEGER :: nerror(2)     ,necpt(100)     ,na(7)  , oldpt1 ,oldpt2
 DOUBLE PRECISION :: i00    ,i01   ,i02    ,i03    ,i04    ,  &
     i10    ,i11   ,i12    ,i13    ,i14    ,  &
     i20    ,i21   ,i22    ,i23    ,i24    , i31   ,i32    ,i33    ,i34    ,  &
     i42    ,i43    ,i44    , i52    ,i53    ,i54    ,  &
     constd ,       i62    ,i63    ,i64
 DOUBLE PRECISION :: kqn(10,10)    ,kqx(10,10)     ,kqe(10,10)     ,  &
     kqy(10,10)    ,fac(7),h(120)  ,h11    ,h12    ,  &
     h13    ,h14   ,h15    ,h16    ,h17    ,h18    ,  &
     h19    ,h1ten ,determ ,pi     ,one    ,huq    ,  &
     integ  ,kij   ,nspopi ,hyq    ,temp60 ,hyqf   ,  &
     za     ,e11   ,d11    ,zb     ,e12    ,d12    ,  &
     a      ,e22   ,d22    ,b      ,e33    ,d33    ,  &
     SIGN   ,t     ,cp     ,ra     ,ts     ,sp     ,  &
     rb     ,n     ,cp2    ,rasq   ,n2     ,sp2    ,  &
     rbsq   ,sl    ,nsp    ,tn     ,l2     ,ncp    ,  &
     piovb  ,dl    ,spe12  ,td     ,temp   ,cpe12  ,  &
     n2e22  ,twod33,tnsp   ,n2e33  ,opi    ,oq     ,  &
     spe22  ,temp1 ,temp5  ,cpe22  ,temp2  ,temp6  ,  &
     sp2e22 ,temp3 ,temp7  ,cp2e22 ,temp4  ,sp2e33 , n2d33  ,sp2d22,spe33
 DOUBLE PRECISION :: sum    ,qq1   ,qq2    ,qq3    ,qq4
 
 COMMON /condad/  constd(5)
 COMMON /matin /  matid  ,inflag,eltemp ,stress ,sinth  ,costh
 COMMON /matout/  g11    ,g12   ,g13    ,g22    ,g23    ,g33    ,  &
     dum(5) ,gsube
 COMMON /sma1io/  dum1(10)      ,ifkgg  ,dum2   ,if4gg
 COMMON /sma1cl/  iopt4  ,k4ggsw,npvt   ,dumcl(7)       ,link(10),  &
     idetck ,dodet ,nogo
 COMMON /sma1et/  ecpt(100)
 COMMON /sma1dp/  integ(28)     ,kij(36),huq(100)       ,hyqf(10),  &
     hyq(10),temp60(60)    ,opi    ,za     ,e11     ,  &
     cp     ,spe22 ,zb     ,e12    ,sp     ,cpe22   ,  &
     a      ,e22   ,cp2    ,sp2e22 ,b      ,e33     ,  &
     sp2    ,cp2e22,SIGN   ,t      ,d11    ,temp1   ,  &
     ra     ,ts    ,d12    ,temp2  ,rb     ,n       ,  &
     d22    ,temp3 ,rasq   ,n2     ,d33    ,temp4   ,  &
     rbsq   ,sl    ,nsp    ,temp5  ,tn     ,l2      ,  &
     ncp    ,temp6 ,piovb  ,dl     ,spe12  ,temp7   ,  &
     td     ,temp  ,cpe12  ,oq     ,n2e22  ,twod33  , tnsp   ,n2e33  ,sp2e33,spe33
 EQUIVALENCE      (constd(1),pi  ), (ecpt(4),matid1),  &
     (ecpt(6),matid2), (ecpt(8),matid3), (ecpt(1),necpt(1))
 EQUIVALENCE      (g,g12), (kqn(1,1),kqe(1,1),kqx(1,1),kqy(1,1))
 EQUIVALENCE      (hyq(1),h11), (hyq(2),h12), (hyq(3),h13),  &
     (hyq(4),h14), (hyq(5),h15), (hyq(6),h16),  &
     (hyq(7),h17), (hyq(8),h18), (hyq(9),h19), (hyq(10),h1ten)
 EQUIVALENCE      (i00,integ( 1)), (i20,integ(11)),  &
     (i01,integ( 2)), (i21,integ(12)), (i02,integ( 3)), (i22,integ(13)),  &
     (i03,integ( 4)), (i23,integ(14)), (i04,integ( 5)), (i24,integ(15)),  &
     (i10,integ( 6)), (i31,integ(16)), (i11,integ( 7)), (i32,integ(17)),  &
     (i12,integ( 8)), (i33,integ(18)), (i13,integ( 9)), (i34,integ(19)),  &
     (i14,integ(10)), (i52,integ(23)), (i42,integ(20)), (i53,integ(24)),  &
     (i43,integ(21)), (i54,integ(25)), (i44,integ(22)), (i62,integ(26)),  &
     (i63,integ(27)), (i64,integ(28))
 DATA    oldpt1,  oldpt2 / 0, 0  /
 DATA    fac   /  1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,120.0D0,720.0D0 /
 DATA    na    /  1,1,1,2,3,3,3  /
 DATA    one   /  1.0D0  /
 
!     DOES PIVOT POINT EQUAL EITHER OF THE LAST TWO SILS
 
 IF (oldpt1 == necpt(2)) IF (oldpt2-necpt(3)) 10,110,10
 IF (oldpt2 == necpt(2)) IF (oldpt1-necpt(3)) 10,110,10
 10 CONTINUE
 
!     NO MATCH THUS DO ENTIRE COMPUTATION
 
 sinth = 0.0
 costh = 1.0
 nint  = necpt(1) - (necpt(1)/1000)*1000 - 1
 n     = nint
 ra    = ecpt(28)
 za    = ecpt(29)
 rb    = ecpt(32)
 zb    = ecpt(33)
 temp1 = rb - ra
 temp2 = zb - za
 sl    = DSQRT(temp1**2 + temp2**2)
 l2    = sl*sl
 IF (sl == 0.0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 20 nerror(1) = necpt(1)/1000
 nerror(2) = n + .3
 CALL mesage (30,39,nerror(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
 30 sp = temp1/sl
 cp = temp2/sl
 a  = ra
 b  = sp
 IF (b == 0.0) THEN
   GO TO    40
 ELSE
   GO TO    60
 END IF
 
!     GO TO 40 FOR B = 0
 
!                               1-N
!                         PI  RA     M+1
!     FOR B = 0,   I   = --------- SL    (FOR ALL M,N .GE. 0)
!                   M,N    M + 1
 
 40 idx = 0
 DO  i = 1,7
   nbegin = na(i)
   
   DO  j = nbegin,5
     
!     M = I - 1
!     N = J - 1
!     MPLUS1 THUS EQUALS I
     
     idx = idx + 1
     integ(idx) = (pi*sl**i)/(DBLE(FLOAT(i))*ra**(j-2))
   END DO
 END DO
 GO TO 100
 
!     ABOVE COMPLETES ALL INTEGRALS FOR B = 0
 
 60 CONTINUE
 
!     FOR B .NE. ZERO
 
!     IF AN OVERFLOW RESULTS BELOW POSSIBLY B IS NOT ZERO, BUT SMALL
 
!     OK BELOW IS FOR B NOT EQUAL TO ZERO
 
!     FIRST M = 0 CASE
 
!                             2-N     2-N
!                       PI (RB    - RA   )
!               I    = --------------------   (N NOT EQUAL TO 2)
!                0,N        (2-N) B
 
 
!     FOR N=2   I    = PI * (LOG RB  -  LOG RA) / B
!                0,2            E          E
 
 rasq  = ra*ra
 rbsq  = rb*rb
 piovb = pi/b
 
 integ(1) = 0.5D0*piovb*(rbsq - rasq)
 integ(2) = piovb*(rb - ra)
 integ(3) = piovb*DLOG(rb/ra)
 integ(4) =-piovb*(one/rb - one/ra)
 integ(5) =-0.5D0*piovb*(one/rbsq - one/rasq)
 
 idx  = 5
 DO  i = 1,6
   mplus1 = i + 1
   nbegin = na(mplus1)
   DO  j = nbegin,5
     
!     M = I
!     N = J - 1
     
!     WE ARE GETTING INTEGRAL(M,N)
!     M = POWER OF S
!     N = POWER OF R
     
!     EVALUATING AT R = RB,  THEN AT R = RA
     
!                                    K   MNK2
!                (M)FAC.     M   (-A) (R)
!     I  = (PI) (-------) ((SUM -------------------------) + (TERM-X))
!      MN          (M+1)    K=0  (M-K)FAC. (K)FAC. (MNK2)
!                 B        (FOR K.NE. MN2                   (FOR K=MN2)
     
!       WHERE    MNK2 = M-N-K+2
!                MN2  = M-N  +2
!             (X)FAC. = X!
!                             MN2
!                         (-A)    LOG(R)
!              TERM-X = --------------------
!                       (M-N+2)FAC. (N-2)FAC.
     
!     NOTE IN DATA STATEMENT THAT 0 FACTORIAL = FAC(1)
!                                 1 FACTORIAL = FAC(2)
!                                 2 FACTORIAL = FAC(3)    ETC.
     
     sum  = 0.0
     SIGN =-1.0D0
     DO  kk = 1,mplus1
       SIGN =-SIGN
       k    = kk - 1
       mn2  = i - j + 3
       qq1  = a
       qq2  = rb
       qq3  = ra
       IF (k == mn2) GO TO 70
       mnk2 = mn2 - k
       mk1  = mplus1 - k
       temp = mnk2
       
!     QQ4  = A**K*(RB**MNK2-RA**MNK2)/(FAC(MK1)*FAC(KK)*TEMP)
       
       qq1  = qq1**k
       qq2  = qq2**mnk2
       qq3  = qq3**mnk2
       qq2  = qq2 - qq3
       qq3  = fac(mk1)*fac(kk)*temp
       GO TO 75
       
!     QQ4 = A**MN2*DLOG(RB/RA)/(FAC(MN2+1)*FAC(J-2))
       
       70 qq1 = qq1**mn2
       qq3 = qq2/qq3
       qq2 = DLOG(qq3)
       qq3 = fac(mn2+1)*fac(j-2)
       75 qq4 = qq1*qq2/qq3
       sum = sum + SIGN*qq4
     END DO
     
     qq1 = pi*fac(mplus1)
     qq2 = b
     qq3 = qq2**mplus1
     qq4 = sum*qq1/qq3
     idx = idx + 1
     integ(idx) = DBLE(qq4)
   END DO
 END DO
 
 100 oldpt1 = necpt(2)
 oldpt2 = necpt(3)
 GO TO 140
 
!     WE HAVE A MATCH ON OLD SIL NUMBER 1
 
 110 IF (npvt-oldpt1 == 0) THEN
   GO TO   120
 ELSE
   GO TO   130
 END IF
 120 npivot = 1
 GO TO 410
 
!     WE HAVE A MATCH ON OLD SIL NUMBER 2
 
 130 npivot = 2
 GO TO 410
 
!     ZERO OUT THE KQN MATRIX
 
 140 DO  i = 1,10
   DO  j = 1,10
     kqn(i,j) = 0.0D0
   END DO
 END DO
 
!     IF MEMBRANE THICKNESS IS NOT ZERO FORM THE KQE MATRIX
 
 t = ecpt(5)
 IF (t == 0.0) THEN
   GO TO   200
 END IF
 160 ASSIGN 190 TO iretrn
 matid  = matid1
 170 inflag = 2
 180 eltemp = ecpt(35)
 CALL mat (ecpt(1))
 GO TO iretrn, (190,230,242)
 190 e11 = g11
 e12 = g12
 e22 = g22
 e33 = g33
 tn  = t * n
 cp2 = cp* cp
 sp2 = sp* sp
 n2  = n * n
 cp2e22= cp2* e22
 sp2e22= sp2* e22
 cpe22 = cp * e22
 spe22 = sp * e22
 cpe12 = cp * e12
 spe12 = sp * e12
 n2e33 = n2 * e33
 n2e22 = n2 * e22
 sp2e33= sp2* e33
 spe33 = sp * e33
 
! /// FURTHER REDUCTION IS NEEDED HERE ///
 
 kqe(1,1) = t*(n2e22 + sp2e33)*i02
 kqe(1,2) = t*(n2e22*i12 - spe33*i01 + sp2e33*i12)
 temp     = e22 + e33
 tnsp     = tn*sp
 kqe(1,3) = tnsp*temp*i02
 kqe(1,4) = tn*(e12*i01 + sp*temp*i12)
 temp     = tn*cp*e22
 kqe(1,5) = temp*i02
 kqe(1,6) = temp*i12
 kqe(1,7) = temp*i22
 kqe(1,8) = temp*i32
 temp4    = 2.d0*sp*i11
 kqe(2,2) = t *(n2e22*i22 + e33*(i00-temp4 + sp2*i22))
 kqe(2,3) = tn*(spe22*i12 - e33*i01 + spe33*i12)
 kqe(2,4) = tn*(e12*i11 + spe22*i22 - e33*i11 + spe33*i22)
 kqe(2,5) = kqe(1,6)
 kqe(2,6) = kqe(1,7)
 kqe(2,7) = kqe(1,8)
 kqe(2,8) = tn*cpe22 *i42
 kqe(3,3) = t *(sp2e22*i02 + n2e33*i02)
 kqe(3,4) = t *(spe12 *i01 + sp2e22*i12 + n2e33*i12)
 temp     = t *cp*spe22
 kqe(3,5) = temp*i02
 kqe(3,6) = temp*i12
 kqe(3,7) = temp*i22
 kqe(3,8) = temp*i32
 kqe(4,4) = t *(e11*i00 + temp4*e12 + sp2e22*i22 + n2e33*i22)
 temp     = sp*cpe22
 kqe(4,5) = t *(cpe12*i01 + temp*i12)
 kqe(4,6) = t *(cpe12*i11 + temp*i22)
 kqe(4,7) = t *(cpe12*i21 + temp*i32)
 kqe(4,8) = t *(cpe12*i31 + temp*i42)
 temp     = t *cp2e22
 kqe(5,5) = temp*i02
 kqe(5,6) = temp*i12
 kqe(5,7) = temp*i22
 kqe(5,8) = temp*i32
 kqe(6,6) = kqe(5,7)
 kqe(6,7) = kqe(5,8)
 kqe(6,8) = temp*i42
 kqe(7,7) = kqe(6,8)
 kqe(7,8) = temp*i52
 kqe(8,8) = temp*i62
 
 200 IF (ecpt(7) == 0.0) GO TO 270
 
!     NOW GET G MATERIAL MATRIX ID = MATID2
 
 matid = matid2
 ASSIGN 230 TO iretrn
 GO TO 170
 
!     NOW FORM D = I DOT G
 
 230 d11 = ecpt(7)*g11
 d12 = ecpt(7)*g12
 d22 = ecpt(7)*g22
 d33 = ecpt(7)*g33
 
!     IF SHEAR THICKNESS IS NOT ZERO FORM THE HYQ AND KQY MATRICES
 
 ts = ecpt(9)
 IF (ts == 0.0) THEN
   GO TO   265
 END IF
 240 CONTINUE
 
!     GET G FOR MATID3
 
 matid  = matid3
 inflag = 1
 ASSIGN 242 TO iretrn
 GO TO 180
 
 242 CONTINUE
 IF (g == 0.0) GO TO 261
 
!     FORMING 1.0/Q DIRECTLY
 
 opi = one / pi
 
! /// MAKE SURE ALL BASIC PRODUCTS ARE AT TOP BEFORE ANY SKIPS
 
 n2d33  = n2 *d33
 sp2d22 = sp2*d22
 oq     = sl*ts*DBLE(g)*(ra+rb)*0.5D0 + i02*(n2d33+sp2d22)*opi
 oq     = one/oq
 nsp    = n*sp
 ncp    = n*cp
 nspopi = nsp  *opi
 twod33 = 2.0D0*d33
 temp1  = d12*(one/rb - one/ra)
 temp2  = nspopi*(d22 + d33)
 temp3  = n  *nspopi*(twod33 + d22)
 temp4  = oq *0.5D0 *ncp*n*d33*opi
 temp5  = opi*(n2*twod33 + sp2*d22)
 temp6  = d12*n2*l2/rb
 temp7  = nspopi*cp*0.5D0
 hyq(1) = oq *(temp1*ncp - temp7*i03*(d33+2.0D0*d22))
 hyq(2) = oq *(ncp*sl/rb*d12 - temp7*i13*(3.0D0*d33+d22)  &
     + 1.5D0*ncp*opi*i02*d33)
 hyq(3) = temp4*i03
 hyq(4) = temp4*i13
 hyq(5) = oq*(temp1*n2 - temp3*i03)
 hyq(6) = oq*(d12*n2*sl/rb - temp3*i13 + temp5*i02)
 hyq(7) = oq*(2.0D0*d11*(ra-rb)+temp6+2.0D0*i12*temp5-temp3*i23)
 hyq(8) = oq*(-d11*6.d0*sl*rb+temp6*sl+3.d0*i22*temp5-temp3*i33)
 hyq(9) =-oq*temp2 * i02
 hyq(10)= oq*(n*sl*(d12+d33) - temp2*i12)
 
 temp = ts*DBLE(g)*i00
 DO  i = 1,10
   hyqf(i) = hyq(i)*temp
 END DO
 DO  i = 1,10
   DO  j = i,10
     kqy(i,j) = kqy(i,j) + hyq(i)*hyqf(j)
   END DO
 END DO
 
!     ADD IN TERMS PER EQUATION-90- PAGE -27- MS-28
 
 temp = ts*DBLE(g)
 kqy( 9,10) = kqy( 9,10) + temp*i10
 kqy(10,10) = kqy(10,10) + temp*i20
 kqy( 9, 9) = kqy( 9, 9) + temp*i00
 
!     END OF KQY COMPUTATION
 
 GO TO 265
 261 ts = 0.0D0
 265 CONTINUE
 
!     THE FOLLOWING CODES WERE MOVED HERE FROM KCONEX
 
!     KQX MATRIX FOR SHEAR THICKNESS CONSIDERATION
 
!     (THE FOLLOWING CODE WAS MACHINE GENERATED AND WILL NOT BE SIMPLI-
!     FIED FURTHER UNTIL FORMULATION VERIFICATION IS COMPLETED)
 
 kqx(1, 1) = kqx(1, 1) + cp*cp*i04*(+d22*n**2+2.25D0*d33*sp**2)
 kqx(1, 2) = kqx(1, 2) + cp*cp*(d33*sp*(+2.25D0*sp*i14-2.25D0*i03)  &
     + d22*n*n*i14)
 kqx(1, 3) = kqx(1, 3) + d33*cp*cp*sp*n*i04*(-7.5D-1)
 kqx(1, 4) = kqx(1, 4) + d33*cp*cp*sp*n*i14*(-7.5D-1)
 kqx(1, 5) = kqx(1, 5) + cp*n*i04*(+d22*n**2+3.0D0*d33*sp**2)
 kqx(1, 6) = kqx(1, 6) + cp*n*(sp*(d33*(+3.0D0*sp*i14-3.0D0*i03)  &
     - d22*i03)  + d22*n*n*i14)
 kqx(1, 7) = kqx(1, 7) + cp*n*(sp*(d33*(+3.0D0*sp*i24-6.0D0*i13)  &
     + d22*i13*(-2.0D0)) - 2.0D0*d12*i02 + d22*n**2*i24)
 kqx(1, 8) = kqx(1, 8) + cp*n*(sp*(d33*(+3.0D0*sp*i34-9.0D0*i23)  &
     + d22*i23*(-3.0D0)) - 6.0D0*d12*i12 + d22*n**2*i34)
 kqx(1, 9) = kqx(1, 9) + cp*i03*(+d22*n**2+1.5D0*d33*sp**2)
 kqx(1,10) = kqx(1,10) + cp*(d33*sp*(-1.5D0*i02+1.5D0*sp*i13) + d22*n*n*i13)
 kqx(2, 2) = kqx(2, 2) + cp*cp*(d33*(sp*(i13*(-4.5D0)  &
     + sp*i24*2.25D0) + i02*2.25D0) + d22*n*n*i24)
 kqx(2, 3) = kqx(2, 3) + d33*cp*cp*n*(-7.5D-1*sp*i14+7.5D-1*i03)
 kqx(2, 4) = kqx(2, 4) + d33*cp*cp*n*(-7.5D-1*sp*i24+7.5D-1*i13)
 kqx(2, 5) = kqx(2, 5) + cp*n*(d33*sp*(+3.0D0*sp*i14-3.0D0*i03) + d22*n*n*i14)
 kqx(2, 6) = kqx(2, 6) + cp*n*(d33*(sp*(i13*(-6.0D0)  &
     + sp*i24*3.0D0) + i02*3.0D0) + d22*(-sp*i13+n**2*i24))
 kqx(2, 7) = kqx(2, 7) + cp*n*(d33*(sp*(i23*(-9.0D0)  &
     + sp*i34*3.0D0) + i12*6.0D0)  &
     + d22*(-2.0D0*sp*i23 + n**2*i34) + d12*i12*(-2.0D0))
 kqx(2, 8) = kqx(2, 8) + cp*n*(d33*(sp*(i33*(-1.20D01)  &
     + sp*i44*3.0D0) + i22*9.0D0)  &
     + d22*(-3.0D0*sp*i33+n**2*i44) + d12*i22*(-6.0D0))
 kqx(2, 9) = kqx(2, 9) + cp*(d33*sp*(+1.5D0*sp*i13-1.5D0*i02) + d22*n*n*i13)
 kqx(2,10) = kqx(2,10) + cp*(d33*(sp*(i12*(-3.0D0)+sp*i23*1.5D0)  &
     + i01*1.5D0)+ d22*n*n*i23)
 kqx(3, 3) = kqx(3, 3) + d33*cp*cp*n*n*i04*2.5D-1
 kqx(3, 4) = kqx(3, 4) + d33*cp*cp*n*n*i14*2.5D-1
 kqx(3, 5) = kqx(3, 5) + d33*cp*sp*n*n*i04*(-1.0D0)
 kqx(3, 6) = kqx(3, 6) + d33*cp*n*n*(-sp*i14+i03)
 kqx(3, 7) = kqx(3, 7) + d33*cp*n*n*(-sp*i24+2.0D0*i13)
 kqx(3, 8) = kqx(3, 8) + d33*cp*n*n*(-sp*i34+3.0D0*i23)
 kqx(3, 9) = kqx(3, 9) + d33*cp*sp*n*i03*(-5.0D-1)
 kqx(3,10) = kqx(3,10) + d33*cp*n*(+5.0D-1*i02-5.0D-1*sp*i13)
 kqx(4, 4) = kqx(4, 4) + d33*cp*cp*n*n*i24*2.5D-1
 kqx(4, 5) = kqx(4, 5) + d33*cp*sp*n*n*i14*(-1.0D0)
 kqx(4, 6) = kqx(4, 6) + d33*cp*n*n*(-sp*i24+i13)
 kqx(4, 7) = kqx(4, 7) + d33*cp*n*n*(-sp*i34+2.0D0*i23)
 kqx(4, 8) = kqx(4, 8) + d33*cp*n*n*(-sp*i44+3.0D0*i33)
 kqx(4, 9) = kqx(4, 9) + d33*cp*sp*n*i13*(-5.0D-1)
 kqx(4,10) = kqx(4,10) + d33*cp*n*(+5.0D-1*i12-5.0D-1*sp*i23)
 kqx(5, 5) = kqx(5, 5) + n*n*i04*(+d22*n**2+4.0D0*d33*sp**2)
 kqx(5, 6) = kqx(5, 6) + n*n*(sp*(d33*(+4.0D0*sp*i14-4.0D0*i03)  &
     + d22*i03*(-1.0D0)) + d22*n*n*i14)
 kqx(5, 7) = kqx(5, 7) + n*n*(sp*(d33*(+4.0D0*sp*i24-8.0D0*i13)  &
     + d22*i13*(-2.0D0))-2.0D0*d12*i02 + d22*n**2*i24)
 kqx(5, 8) = kqx(5, 8) + n*n*(sp*(d33*(+4.0D0*sp*i34-1.20D01*i23)  &
     + d22*i23*(-3.0D0)) - 6.0D0*d12*i12 + d22*n**2*i34)
 kqx(5, 9) = kqx(5, 9) + n*i03*(+d22*n**2+2.0D0*d33*sp**2)
 kqx(5,10) = kqx(5,10) + n*(d33*sp*(-2.0D0*i02+2.0D0*sp*i13) + d22*n*n*i13)
 kqx(6, 6) = kqx(6, 6) + n*n*(sp*(i13*(d22*(-2.0D0)+d33*(-8.0D0))  &
     + d33*sp*i24*4.0D0) + d22*n**2*i24 + 4.0D0*d33*i02) + d22*sp*sp*i02
 kqx(6, 7) = kqx(6, 7) + n*n*(sp*(i23*(d22*(-3.0D0)+d33*(-1.20D01))  &
     + d33*sp*i34*4.0D0) + i12*(-2.0D0*d12+8.0D0*d33)  &
     + d22*n*n*i34) + sp*(+2.0D0*d12*i01+2.0D0*d22*sp*i12)
 kqx(6, 8) = kqx(6, 8) + n*n*(sp*(i33*(d22*(-4.0D0)+d33*(-1.6D01))  &
     + d33*sp*i44*4.0D0) + i22*(-6.0D0*d12+1.20D01*d33)  &
     + d22*n*n*i44)+sp*(+6.0D0*d12*i11+3.0D0*d22*sp*i22)
 kqx(6, 9) = kqx(6, 9) + n*(sp*(d33*(+2.0D0*sp*i13-2.0D0*i02)  &
     + d22*i02*(-1.0D0)) + d22*n*n*i13)
 kqx(6,10) = kqx(6,10) + n*(d33*(sp*(i12*(-4.0D0) + sp*i23*2.0D0)  &
     + i01*2.0D0)+ d22*(+n**2*i23-sp*i12))
 kqx(7, 7) = kqx(7, 7) + n*n*(sp*(i33*(d22*(-4.0D0)+d33*(-1.6D01))  &
     + d33*sp*i44*4.0D0) + i22*(d12*(-4.0D0) +d33*1.6D01)  &
     + d22*n*n*i44) + sp*(d12*i11*8.0D0+d22*sp*i22*4.0D0) + d11*i00*4.0D0
 kqx(7, 8) = kqx(7, 8) + n*n*(sp*(i43*(d22*(-5.0D0)+d33*(-2.0D01))  &
     + d33*sp*i54*4.0D0) + i32*(d12*(-8.0D0)+d33*2.40D01)  &
     + d22*n*n*i54) + sp*(d12*i21*1.80D01+d22*sp*i32*6.0D0) + d11*i10*1.20D01
 kqx(7, 9) = kqx(7, 9) + n*(sp*(d33*(+2.0D0*sp*i23-4.0D0*i12)  &
     + d22*i12*(-2.0D0)) - 2.0D0*d12*i01 + d22*n**2*i23)
 kqx(7,10) = kqx(7,10) + n*(d33*(sp*(i22*(-6.0D0)+sp*i33*2.0D0)  &
     + i11*4.0D0)+ d22*(+n**2*i33-2.0D0*sp*i22) + d12*i11*(-2.0D0))
 kqx(8, 8) = kqx(8, 8) + n*n*(sp*(i53*(d22*(-6.0D0)+d33*(-2.40D01))  &
     + d33*sp*i64*4.0D0) + i42*(d12*(-1.20D01) + d33*3.60D01)  &
     + d22*n*n*i64) + sp*(d12*i31*3.60D01+d22*sp*i42*9.0D0) + d11*i20*3.60D01
 kqx(8, 9) = kqx(8, 9) + n*(sp*(d33*(+2.0D0*sp*i33-6.0D0*i22)  &
     + d22*i22*(-3.0D0)) - 6.0D0*d12*i11 + d22*n**2*i33)
 kqx(8,10) = kqx(8,10) + n*(d33*(sp*(i32*(-8.0D0)+sp*i43*2.0D0)  &
     + i21*6.0D0)+ d22*(+n**2*i43-3.0D0*sp*i32) + d12*i21*(-6.0D0))
 kqx(9, 9) = kqx(9, 9) + i02*(+d22*n**2+d33*sp**2)
 kqx(9,10) = kqx(9,10) + d33*sp*(-i01+sp*i12) + d22*n*n*i12
 kqx(10,10)= kqx(10,10)+ d33*(sp*(i11*(-2.0D0)+ sp*i22)+i00) + d22*n*n*i22
 IF (ts == 0.0D0) GO TO 270
 
!     THE FOLLOWING CODES WERE MOVED HERE FROM KCONEY
 
 kqx(1, 1) = kqx(1, 1) + h11*(sp*(cp*n*i03*(d22*2.0D0+d33*3.0D0)  &
     + d22*sp*h11*i02) + d33*n*n*h11*i02)
 kqx(1, 2) = kqx(1, 2) + n*(cp*(sp*(d22*(+h12*i03+h11*i13)  &
     + d33*(+1.5D0*h12*i03+1.5D0*h11*i13))  &
     + d33*h11*i02*(-1.5D0))+d33*n*h11*h12*i02) + d22*sp*sp*h11*h12*i02
 kqx(1, 3) = kqx(1, 3) + n*(d33*(cp*i03*(+1.5D0*sp*h13  &
     - 5.0D-1*n*h11) + n*h11*h13*i02) + d22*cp*sp*h13*i03) + d22*sp*sp*h11*h13*i02
 kqx(1, 4) = kqx(1, 4) + n*(d33*(cp*(+1.5D0*sp*h14*i03  &
     - 5.0D-1*n*h11*i13)+n*h11*h14*i02) + d22*cp*sp*h14*i03)  &
     + d22*sp*sp*h11*h14*i02
 kqx(1, 5) = kqx(1, 5) + sp*(n*i03*(d22*(+cp*h15+n*h11)  &
     + d33*(+1.5D0*cp*h15+2.0D0*n*h11))  &
     + d22*sp*h11*h15*i02) + d33*n*n*h11*h15*i02
 kqx(1, 6) = kqx(1, 6) + sp*(d22*(h11*(sp*i02*(-1.0D0+h16)  &
     + n*n*i13)  + cp*n*h16*i03)+d33*n*(+1.5D0*cp*h16*i03  &
     + 2.0D0*n*h11*i13)) + d33*n*n*h11*i02*(-2.0D0+h16)
 kqx(1, 7) = kqx(1, 7) + sp*(h11*(d22*(sp*(-2.0D0*i12+h17*i02)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + cp*n*h17*i03*(+d22+1.5D0*d33)) + d33*n*n*h11*(-4.0D0*i12+h17*i02)
 kqx(1, 8) = kqx(1, 8) + sp*(h11*(d22*(sp*(-3.0D0*i22+h18*i02)  &
     + n*n*i33)  - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + cp*n*h18*i03*(+d22+1.5D0*d33)) + d33*n*n*h11*(-6.0D0*i22+h18*i02)
 kqx(1, 9) = kqx(1, 9) + sp*(n*(d22*(+cp*h19*i03+h11*i02)  &
     + d33*(+1.5D0*cp*h19*i03+h11*i02))  &
     + d22*sp*h11*h19*i02) + d33*n*n*h11*h19*i02
 kqx(1,10) = kqx(1,10) + n*(d33*(h11*(-i01+sp*i12+n*h1ten*i02)  &
     + cp*sp*h1ten*i03*1.5D0) + d22*sp*(+cp*h1ten*i03  &
     + h11*i12)) + d22*sp*sp*h11*h1ten*i02
 kqx(2, 2) = kqx(2, 2) + h12*(n*(cp*(d33*(sp*i13*3.d0+i02*(-3.d0))  &
     + d22*sp*i13*2.d0) + d33*n*h12*i02) + d22*sp*sp*h12*i02)
 kqx(2, 3) = kqx(2, 3) + n*(d33*(cp*(h13*(+1.5D0*sp*i13-1.5D0*i02)  &
     + n*h12*i03*(-5.0D-1)) + n*h12*h13*i02)  &
     + d22*cp*sp*h13*i13) + d22*sp*sp*h12*h13*i02
 kqx(2, 4) = kqx(2, 4) + n*(d33*(cp*(h14*(+1.5D0*sp*i13-1.5D0*i02)  &
     + n*h12*i13*(-5.0D-1)) + n*h12*h14*i02)  &
     + d22*cp*sp*h14*i13) + d22*sp*sp*h12*h14*i02
 kqx(2, 5) = kqx(2, 5) + n*(d33*(h15*(cp*(+1.5D0*sp*i13-1.5D0*i02)  &
     + n*h12*i02)+ sp*n*h12*i03*2.0D0) + d22*sp*(+cp*h15*i13  &
     + n*h12*i03)) + d22*sp*sp*h12*h15*i02
 kqx(2, 6) = kqx(2, 6) + n*(d33*(n*h12*(i02*(-2.0D0+h16)  &
     + sp*i13*2.0D0) + cp*h16*(+1.5D0*sp*i13-1.5D0*i02))  &
     + d22*sp*i13*(+cp*h16+n*h12)) + d22*sp*sp*h12*i02*(-1.0D0+h16)
 kqx(2, 7) = kqx(2, 7) + sp*(h12*(d22*(sp*(-2.0D0*i12+h17*i02)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + cp*n*h17*i13*(+d22+1.5D0*d33)) + d33*n*(n*h12*(-4.0D0*i12+h17*i02)  &
     + cp*h17*i02*(-1.5D0))
 kqx(2, 8) = kqx(2, 8) + sp*(h12*(d22*(sp*(-3.0D0*i22+h18*i02)  &
     + n*n*i33)  - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + cp*n*h18*i13*(+d22+1.5D0*d33)) + d33*n*(n*h12*(-6.0D0*i22+h18*i02)  &
     + cp*h18*i02*(-1.5D0))
 kqx(2, 9) = kqx(2, 9) + n*(d33*(h19*(cp*(+1.5D0*sp*i13-1.5D0*i02)  &
     + n*h12*i02)+ sp*h12*i02)+d22*sp*(+cp*h19*i13+h12*i02))  &
     + d22*sp*sp*h12*h19*i02
 kqx(2,10) = kqx(2,10) + n*(d33*(h12*(-i01+sp*i12+n*h1ten*i02)  &
     + cp*h1ten*(+1.5D0*sp*i13-1.5D0*i02)) + d22*sp*(+cp*h1ten*i13+h12*i12))  &
     + d22*sp*sp*h12*h1ten*i02
 kqx(3, 3) = kqx(3, 3) + h13*(d33*n*n*(cp*i03*(-1.0D0)+h13*i02)  &
     + d22*sp*sp*h13*i02)
 kqx(3, 4) = kqx(3, 4) + d33*n*n*(cp*(-5.0D-1*h14*i03-5.0D-1*h13  &
     * i13)+h13*h14*i02) + d22*sp*sp*h13*h14*i02
 kqx(3, 5) = kqx(3, 5) + n*n*(d33*(h13*(+2.0D0*sp*i03+h15*i02)  &
     + cp*h15*i03*(-5.0D-1)) + d22*sp*h13*i03) + d22*sp*sp*h13*h15*i02
 kqx(3, 6) = kqx(3, 6) + h13*(sp*(d22*(sp*i02*(-1.d0+h16)+n*n*i13)  &
     + d33*n*n*i13*2.0D0) + d33*n*n*i02*(-2.0D0+h16))  &
     + d33*cp*n*n*h16*i03*(-5.0D-1)
 kqx(3, 7) = kqx(3, 7) + h13*(sp*(d22*(sp*(-2.0D0*i12+h17*i02)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + d33*n*n*(-4.0D0*i12+h17*i02)) + d33*cp*n*n*h17*i03*(-5.0D-1)
 kqx(3, 8) = kqx(3, 8) + h13*(sp*(d22*(sp*(-3.0D0*i22+h18*i02)  &
     + n*n*i33)  - 6.0D0*d12*i11+2.0D0*d33*n**2*i33)  &
     + d33*n*n*(-6.0D0*i22+h18*i02)) + d33*cp*n*n*h18*i03*(-5.0D-1)
 kqx(3, 9) = kqx(3, 9) + n*(d33*(n*h19*(-5.0D-1*cp*i03+h13*i02)  &
     + sp*h13*i02) + d22*sp*h13*i02) + d22*sp*sp*h13*h19*i02
 kqx(3,10) = kqx(3,10) + n*(d33*(h13*(-i01+sp*i12+n*h1ten*i02)  &
     + cp*n*h1ten*i03*(-5.0D-1))+d22*sp*h13*i12) + d22*sp*sp*h13*h1ten*i02
 kqx(4, 4) = kqx(4, 4) + h14*(d33*n*n*(cp*i13*(-1.0D0)+h14*i02)  &
     + d22*sp*sp*h14*i02)
 kqx(4, 5) = kqx(4, 5) + n*n*(d33*(h14*(+2.0D0*sp*i03+h15*i02)  &
     + cp*h15*i13*(-5.0D-1)) + d22*sp*h14*i03) + d22*sp*sp*h14*h15*i02
 
!     THE FOLLOWING CODES, THRU 270, WERE MOVED HERE FROM KCONEZ
 
 kqx(4, 6) = kqx(4 ,6) + h14*(sp*(d22*(sp*i02*(-1.d0+h16)+n*n*i13)  &
     + d33*n*n*i13*2.0D0) + d33*n*n*i02*(-2.0D0+h16))  &
     + d33*cp*n*n*h16*i13*(-5.0D-1)
 kqx(4, 7) = kqx(4, 7) + h14*(sp*(d22*(sp*(-2.0D0*i12+h17*i02)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + d33*n*n*(-4.0D0*i12+h17*i02)) + d33*cp*n*n*h17*i13*(-5.0D-1)
 kqx(4, 8) = kqx(4, 8) + h14*(sp*(d22*(sp*(-3.0D0*i22+h18*i02)  &
     + n*n*i33)  - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + d33*n*n*(-6.0D0*i22+h18*i02)) + d33*cp*n*n*h18*i13*(-5.0D-1)
 kqx(4, 9) = kqx(4, 9) + n*(d33*(n*h19*(-5.0D-1*cp*i13+h14*i02)  &
     + sp*h14*i02)+d22*sp*h14*i02)+d22*sp*sp*h14*h19*i02
 kqx(4,10) = kqx(4,10) + n*(d33*(h14*(-i01+sp*i12+n*h1ten*i02)  &
     + cp*n*h1ten*i13*(-5.0D-1)) + d22*sp*h14*i12) + d22*sp*sp*h14*h1ten*i02
 kqx(5, 5) = kqx(5, 5) + h15*(sp*(n*n*i03*(d22*2.0D0+d33*4.0D0)  &
     + d22*sp*h15*i02) + d33*n*n*h15*i02)
 kqx(5, 6) = kqx(5, 6) + sp*(d22*(h15*(sp*i02*(-1.d0+h16)+n*n*i13)  &
     + n*n*h16*i03) + d33*n*n*(+2.0D0*h16*i03+2.0D0*h15*i13))  &
     + d33*n*n*h15*i02*(-2.0D0+h16)
 kqx(5, 7) = kqx(5, 7) + sp*(h15*(d22*(sp*(-2.0D0*i12+h17*i02)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + n*n*h17*i03*(+d22+2.0D0*d33)) + d33*n*n*h15*(-4.0D0*i12+h17*i02)
 kqx(5, 8) = kqx(5, 8) + sp*(h15*(d22*(sp*(-3.0D0*i22+h18*i02)  &
     + n*n*i33) - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + n*n*h18*i03*(+d22+2.0D0*d33)) + d33*n*n*h15*(-6.0D0*i22+h18*i02)
 kqx(5, 9) = kqx(5, 9) + sp*(n*(d22*(+n*h19*i03+h15*i02)  &
     + d33*(+2.0D0*n*h19*i03+h15*i02)) + d22*sp*h15*h19*i02) + d33*n*n*h15*h19*i02
 kqx(5,10) = kqx(5,10) + n*(d33*(h15*(-i01+sp*i12+n*h1ten*i02)  &
     + sp*n*h1ten*i03*2.d0) + d22*sp*(+n*h1ten*i03+h15*i12))  &
     + d22*sp*sp*h15*h1ten*i02
 kqx(6, 6) = kqx(6, 6) + h16*(sp*(d22*(sp*i02*(-2.0D0+h16)  &
     + n*n*i13*2.0D0) + d33*n*n*i13*4.0D0) + d33*n*n*i02*(-4.0D0+h16))
 kqx(6, 7) = kqx(6, 7) + sp*(d22*(sp*(h16*(-2.0D0*i12+h17*i02)  &
     + h17*i02*(-1.0D0)) + n*n*(+h17*i13+h16*i23))  &
     + d33*n*n*(+2.0D0*h17*i13 + 2.0D0*h16*i23)  &
     + d12*h16*i01*(-2.0D0))+d33*n*n*(h16*(-4.0D0*i12  &
     + h17*i02)  + h17*i02*(-2.0D0))
 kqx(6, 8) = kqx(6, 8) + sp*(d22*(sp*(h16*(-3.0D0*i22+h18*i02)  &
     + h18*i02*(-1.0D0)) + n*n*(+h18*i13+h16*i33))  &
     + d33*n*n*(+2.0D0*h18*i13 + 2.0D0*h16*i33)  &
     + d12*h16*i11*(-6.0D0)) + d33*n*n*(h16*(-6.0D0*i22  &
     + h18*i02)  + h18*i02*(-2.0D0))
 kqx(6, 9) = kqx(6, 9) + sp*(d22*(h19*(sp*i02*(-1.d0+h16)+n*n*i13)  &
     + n*h16*i02)+ d33*n*(+2.0D0*n*h19*i13+h16*i02))  &
     + d33*n*n*h19*i02*(-2.0D0+h16)
 kqx(6,10) = kqx(6,10) + n*(d33*(n*h1ten*(i02*(-2.0D0+h16)  &
     + sp*i13*2.0D0) + h16*(-i01+sp*i12)) + d22*sp*(+n*h1ten*i13+h16*i12))  &
     + d22*sp*sp*h1ten*i02*(-1.0D0+h16)
 kqx(7, 7) = kqx(7, 7) + h17*(sp*(d22*(sp*(i12*(-4.0D0)+h17*i02)  &
     + n*n*i23*2.0D0) + d12*i01*(-4.0D0)+d33*n*n*i23*4.0D0)  &
     + d33*n*n*(i12*(-8.0D0)+h17*i02))
 kqx(7, 8) = kqx(7, 8) + sp*(d22*(sp*(h17*(-3.0D0*i22+h18*i02)  &
     + h18*i12*(-2.0D0)) + n*n*(+h18*i23+h17*i33))  &
     + d12*(-6.0D0*h17*i11-2.0D0*h18*i01)  &
     + d33*n*n*(+2.0D0*h18*i23 + 2.0D0*h17*i33))  &
     + d33*n*n*(h17*(-6.0D0*i22+h18*i02) + h18*i12*(-4.0D0))
 kqx(7, 9) = kqx(7, 9) + sp*(h19*(d22*(sp*(+h17*i02-2.0D0*i12)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + n*h17*i02*(+d22+d33))+d33*n*n*h19*(-4.d0*i12+h17*i02)
 kqx(7,10) = kqx(7,10) + sp*(h1ten*(d22*(sp*(+h17*i02-2.0D0*i12)  &
     + n*n*i23)  - 2.0D0*d12*i01 + 2.0D0*d33*n**2*i23)  &
     + n*h17*i12*(+d22+d33))+d33*n*(n*h1ten*(-4.0D0*i12  &
     + h17*i02)  + h17*i01*(-1.0D0))
 kqx(8, 8) = kqx(8, 8) + h18*(sp*(d22*(sp*(i22*(-6.0D0)+h18*i02)  &
     + n*n*i33*2.0D0) + d12*i11*(-1.2D01)+d33*n*n*i33*4.0D0)  &
     + d33*n*n*(i22*(-1.2D01)+h18*i02))
 kqx(8, 9) = kqx(8, 9) + sp*(h19*(d22*(sp*(+h18*i02-3.0D0*i22)  &
     + n*n*i33)  - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + n*h18*i02*(+d22+d33))+d33*n*n*h19*(-6.d0*i22+h18*i02)
 kqx(8,10) = kqx(8,10) + sp*(h1ten*(d22*(sp*(+h18*i02-3.0D0*i22)  &
     + n*n*i33)  - 6.0D0*d12*i11 + 2.0D0*d33*n**2*i33)  &
     + n*h18*i12*(+d22+d33)) + d33*n*(n*h1ten*(-6.0D0*i22  &
     + h18*i02)  + h18*i01*(-1.0D0))
 kqx(9, 9) = kqx(9, 9) + h19*i02*(sp*(n*(d22*2.0D0+d33*2.0D0)  &
     + d22*sp*h19) + d33*n*n*h19)
 kqx(9,10) = kqx(9,10) + n*(d33*(h19*(-i01+sp*i12+n*h1ten*i02)  &
     + sp*h1ten*i02) + d22*sp*(+h1ten*i02+h19*i12)) + d22*sp*sp*h19*h1ten*i02
 kqx(10,10)= kqx(10,10)+ h1ten*(n*(d33*(sp*i12*2.0D0+i01*(-2.0D0)  &
     + n*h1ten*i02) + d22*sp*i12*2.0D0)+d22*sp*sp*h1ten*i02)
 
!     SET LOWER TRIANGLE EQUAL TO UPPER TRIANGLE OF KQN MATRIX
 
 270 DO  i = 1,10
   DO  j = i,10
     kqn(j,i) = kqn(i,j)
   END DO
 END DO
 
!     FILL HUQ PER PAGE 15 MS-28
 
 DO  i = 1,100
   huq(i) = 0.0D0
 END DO
 huq(  1) = one
 huq( 13) = one
 huq( 25) = one
 huq( 36) = one
 huq( 49) = one
 huq( 51) = one
 huq( 52) = sl
 huq( 63) = one
 huq( 64) = sl
 huq( 75) = one
 huq( 76) = sl
 huq( 77) = l2
 huq( 78) = huq(77)*sl
 huq( 86) = one
 huq( 87) = 2.0D0*sl
 huq( 88) = 3.0D0*huq(77)
 huq(100) = sl
 
 IF (ts == 0.0) THEN
   GO TO   320
 END IF
 300 huq( 41) = cp/ra
 huq( 45) = n /ra
 huq( 91) = cp/rb
 huq( 92) = huq(91)*sl
 huq( 95) = n/rb
 huq( 96) = huq(95)*sl
 huq( 97) = huq(95)*l2
 huq( 98) = huq(96)*l2
 huq( 99) = one
 
!     SUBTRACT FROM ROWS 4 AND 9 OF THE ABOVE MATRIX, THE HYQ MATRIX
 
 DO  i = 1,10
   huq(i+30) = huq(i+30) - hyq(i)
   huq(i+80) = huq(i+80) - hyq(i)
 END DO
 320 CONTINUE
 
!     NO NEED TO CALCULATE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY
 
 ising =-1
 CALL inverd (10,huq(1),10,dum,0,determ,ising,temp60(1))
!     CHECK SINGULARITY
 
 SELECT CASE ( ising )
   CASE (    1)
     GO TO 340
   CASE (    2)
     GO TO 330
 END SELECT
 330 CALL mesage (30,40,necpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
!     NOT SINGULAR, CONTINUE ON..
 
 340 CONTINUE
 IF (ts /= 0.0D0) GO TO 345
 huq( 85) = 0.0D0
 huq(100) = 0.0D0
 345 CONTINUE
 
!                                 T    N       T
!     NOW SOLVE FOR (K  ) = (E)(H  )(K  )(H )(E )    I = PIVOT A OR B
!                     IJ         I    Q    J         J = A,B
 
 
!                             T    N        T  T
!     WE WILL SOLVE FOR (E)(H  )(K  )((E)(H  ))
!                            A    Q        B
 
 
!                            T                      T
!     FIRST GET EHAT = (E)(H  ),  AND  EHBT = (E)(H  )
!                           A                      B
 
 
!     EHAT WILL BE STORED AT H(1)...H(60) AND EHBT AT H(61)...H(120)
 
!                0    SP   CP   0    0
!                1    0    0    0    0
!                0    CP  -SP   0    0
!     MATRIX E = 0    0    0    0    SP
!                0    0    0    1    0
!                0    0    0    0    CP
 
 inc1 = 0
 inc2 = 0
 350 DO  i = 1,10
   idx  = i + inc1
   iten = 10*i - 9 + inc2
   h(idx   ) = huq(iten+1)*sp + huq(iten+2)*cp
   h(idx+10) = huq(iten  )
   h(idx+20) = huq(iten+1)*cp - huq(iten+2)*sp
   h(idx+30) = huq(iten+4)*sp
   h(idx+40) = huq(iten+3)
   h(idx+50) = huq(iten+4)*cp
 END DO
 IF (inc1 == 0) THEN
   GO TO   370
 ELSE
   GO TO   380
 END IF
 370 inc1 = 60
 inc2 = 5
 GO TO 350
 380 CONTINUE
 
!     DETERMINE PIVOT POINT NUMBER
 
 IF (necpt(2) == npvt) GO TO 390
 IF (necpt(3) == npvt) GO TO 400
 CALL mesage (-30,34,necpt(1))
 390 npivot = 1
 GO TO 410
 400 npivot = 2
 GO TO 410
 
!     EHAT(1) IS AT H( 1)
!     EHBT(1) IS AT H(61)
 
 410 CALL gmmatd (h(60*npivot-59),6,10,0, kqn(1,1),10,10,0, temp60(1))
 
!     IF N = 0 DOUBLE RESULT FOR KIJ
 
 IF (n == 0) THEN
   GO TO   420
 ELSE
   GO TO   440
 END IF
 420 DO  i = 1,60
   temp60(i) = temp60(i)*2.0D0
 END DO
 
 440 DO  j = 1,2
   CALL gmmatd (temp60(1),6,10,0, h(60*j-59),6,10,1, kij(1))
   CALL sma1b (kij(1),necpt(j+1),-1,ifkgg,0.0D0)
   IF (iopt4 == 0) THEN
     GO TO   470
   END IF
   450 IF (gsube == 0.0) THEN
     GO TO   470
   END IF
   460 temp   = gsube
   k4ggsw = 1
   CALL sma1b (kij(1),necpt(j+1),-1,if4gg,temp)
 END DO
 
 RETURN
END SUBROUTINE kconed
