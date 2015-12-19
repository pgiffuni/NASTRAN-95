SUBROUTINE tker (x0,y0,z0,kr,br,sgr,cgr,sgs,cgs,t1,t2,m)
     
!     COMPUTES EITHER THE TOTAL KERNELS (IND=0) USED IN THE CALCULATION
!     OF A FINITE LENGTH DOUBLET LINE,  OR
!     THE INCREMENTAL OSCILLATORY KERNELS (IND=1) USED IN EVALUATING
!     THE INFLUENCE COEFFICIENT MATRIX ELEMENTS
 
 
 REAL, INTENT(IN)                         :: x0
 REAL, INTENT(IN)                         :: y0
 REAL, INTENT(IN)                         :: z0
 REAL, INTENT(IN)                         :: kr
 REAL, INTENT(IN)                         :: br
 REAL, INTENT(IN)                         :: sgr
 REAL, INTENT(IN)                         :: cgr
 REAL, INTENT(IN)                         :: sgs
 REAL, INTENT(IN)                         :: cgs
 REAL, INTENT(OUT)                        :: t1
 REAL, INTENT(OUT)                        :: t2
 REAL, INTENT(IN)                         :: m
 REAL :: i00r,i00i,j00r,j00i,i10r,i10i,i20r3,i20i3,i0ur,  &
     i0ui,j0ur,j0ui,i1ur,i1ui,i2ur3,i2ui3,k1,mu1,mu,k2,  &
     k10,k20,k1rt1,k1it1,k2rt2p,k2it2p,k10t1,k20t2p,kd1r, kd1i, kd2r,kd2i
 COMMON /dlm/ k10,k20,k1rt1,k1it1,k2rt2p,k2it2p,k10t1,k20t2p
 COMMON /kds/ ind,kd1r,kd1i,kd2r,kd2i
 
 eps    = 0.00001
 k10    = 0.0
 k20    = 0.0
 k1rt1  = 0.0
 k1it1  = 0.0
 k2rt2p = 0.0
 k2it2p = 0.0
 k10t1  = 0.0
 k20t2p = 0.0
 r1     = SQRT(y0*y0 + z0*z0)
 r1s    = r1
 IF (ABS(r1) > eps) GO TO 200
 IF (x0 < 0.0) THEN
   GO TO   905
 END IF
 120 c1     = kr*x0/br
 t1     = cgr*cgs + sgr*sgs
 k10    = 2.0
 k1rt1  = 2.0*t1*COS(c1)
 k1it1  =-2.0*t1*SIN(c1)
 k10t1  = 2.0*t1
 GO TO  905
 200 c1     = cgr
 c2     = sgr
 c3     = cgs
 c4     = sgs
 t2p    = (z0*z0*c1*c3 + y0*y0*c2*c4 - z0*y0*(c2*c3+c1*c4))
 t2     = (100.*t2p)/(br*br)
 IF (ABS(t2)-eps < 0.0) THEN
   GO TO   210
 ELSE
   GO TO   220
 END IF
 210 ichuz  = 1
 t1     = cgr*cgs + sgr*sgs
 t2     = 0.0
 GO TO 300
 220 t1     = cgr*cgs + sgr*sgs
 IF (ABS(t1)-eps < 0.0) THEN
   GO TO   230
 ELSE
   GO TO   240
 END IF
 230 ichuz  = 2
 t1     = 0.
 GO TO 300
 240 ichuz  = 3
 300 beta2  = (1.-m*m)
 bigr   = SQRT(x0*x0 + beta2*r1*r1)
 k1     = kr*r1/br
 mu1    = (m*bigr-x0)/(beta2*r1)
 mu     = ABS(mu1)
 k2     = k1*k1
 IF (mu1 < 0) THEN
   GO TO   310
 ELSE IF (mu1 == 0) THEN
   GO TO   320
 ELSE
   GO TO   330
 END IF
 310 ichuz  = ichuz + 3
 GO TO 330
 320 ichuz  = ichuz + 6
 
!     (N*C)**2  FOR  N = 1,11  AND C = .372 =
 
!       .138384      .553536     1.245456      2.214144
!      3.4596       4.981824     6.780816      8.856576
!     11.209104    13.8384      16.744464
 
!     (N*C)  FOR  N = 1,12  AND  14,16,18,20,22   =
 
!       .744        1.116        1.488         1.86      2.232
!      2.604        2.976        3.348         3.72      4.092
!      4.464        5.208        5.952         6.696     7.44
!      8.184
 
!     A(N)  FOR N = 1,11  =
 
!      .24186198   -2.7918027    24.991079    -111.59196
!      271.43549   -305.75288    -41.18363     545.98537
!     -644.78155    328.72755    -64.279511
 
 330 CONTINUE
 exarg = -0.372*mu
 
!     THE FOLLOWING TEST ON THE SIZE OF THE ARGUMENT TO  EXP  IS
!     NEEDED TO AVOID UNDERFLOW IN  SUBPROGRAM  EXP
 
 IF (exarg >= -180.0) GO TO 335
 e   = 0.0
 GO TO 337
 335 e   = EXP(exarg)
 337 CONTINUE
 c1  =  0.138384 + k2
 c2  =  0.553536 + k2
 c3  =  1.245456 + k2
 c4  =  2.214144 + k2
 c5  =  3.4596   + k2
 c6  =  4.981824 + k2
 c7  =  6.780816 + k2
 c8  =  8.856576 + k2
 c9  = 11.209104 + k2
 c10 = 13.8384   + k2
 c11 = 16.744464 + k2
 r1  = .24186198 / c1
 r2  =-2.7918027 / c2
 r3  = 24.991079 / c3
 r4  =-111.59196 / c4
 r5  = 271.43549 / c5
 r6  =-305.75288 / c6
 r7  =-41.18363  / c7
 r8  = 545.98537 / c8
 r9  =-644.78155 / c9
 r10 = 328.72755 / c10
 r11 =-64.279511 / c11
 IF (ichuz < 4) GO TO 340
 i00r = .372*(r1 + 2.*r2 + 3.*r3 + 4.*r4 + 5.*r5 + 6.*r6 + 7.*r7 +  &
     8.*r8 + 9.*r9 + 10.*r10 + 11.*r11)
 i00i =-k1*(r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11)
 340 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 420
   CASE (    2)
     GO TO 350
   CASE (    3)
     GO TO 350
   CASE (    4)
     GO TO 390
   CASE (    5)
     GO TO 350
   CASE (    6)
     GO TO 350
   CASE (    7)
     GO TO 380
   CASE (    8)
     GO TO 350
   CASE (    9)
     GO TO 350
 END SELECT
 350 q1  = r1/ c1
 q2  = r2/ c2
 q3  = r3/ c3
 q4  = r4/ c4
 q5  = r5/ c5
 q6  = r6/ c6
 q7  = r7/ c7
 q8  = r8/ c8
 q9  = r9/ c9
 q10 = r10/c10
 q11 = r11/c11
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 420
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 410
   CASE (    4)
     GO TO 390
   CASE (    5)
     GO TO 360
   CASE (    6)
     GO TO 360
   CASE (    7)
     GO TO 380
   CASE (    8)
     GO TO 360
   CASE (    9)
     GO TO 360
 END SELECT
 360 j00r  = q1*(.138384-k2)+q2*(.553536-k2)+q3*(1.245456-k2)+q4*  &
     (2.214144-k2)+q5*(3.4596-k2)+q6*(4.981824-k2)+q7*(6.780816  &
     -k2)+q8*(8.856576-k2)+q9*(11.209104-k2)+q10*(13.8384-k2)+ q11*(16.744464-k2)
 i20r3 = 2.+k1*i00i+k2*j00r
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 420
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 410
   CASE (    4)
     GO TO 390
   CASE (    5)
     GO TO 410
   CASE (    6)
     GO TO 390
   CASE (    7)
     GO TO 380
   CASE (    8)
     GO TO 370
   CASE (    9)
     GO TO 370
 END SELECT
 370 j00i  = -k1*(.744*q1+1.488*q2+2.232*q3+2.976*q4+3.72*q5+4.464*q6+  &
     5.208*q7+5.952*q8+6.696*q9+7.44*q10+8.184*q11)
 i20i3 = -k1*i00r+k2*j00i
 IF (ichuz == 8) GO TO 500
 380 i10i = -k1*i00r
 390 i10r = 1.+ k1*i00i
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 420
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 410
   CASE (    4)
     GO TO 420
   CASE (    5)
     GO TO 410
   CASE (    6)
     GO TO 410
   CASE (    7)
     GO TO 500
   CASE (    8)
     GO TO 500
   CASE (    9)
     GO TO 500
 END SELECT
 410 j0ur = e*(q1*(0.138384 - k2 + 0.372*mu*c1) +  &
     e*(q2*(0.553536 - k2 + 0.744*mu*c2) +  &
     e*(q3*(1.245456 - k2 + 1.116*mu*c3) +  &
     e*(q4*(2.214144 - k2 + 1.488*mu*c4) +  &
     e*(q5*(3.4596   - k2 + 1.860*mu*c5) +  &
     e*(q6*(4.981824 - k2 + 2.232*mu*c6) +  &
     e*(q7*(6.780816 - k2 + 2.604*mu*c7) +  &
     e*(q8*(8.856576 - k2 + 2.976*mu*c8) +  &
     e*(q9*(11.209104- k2 + 3.348*mu*c9) +  &
     e*(q10*(13.8384 - k2 + 3.72*mu*c10) +  &
     e*(q11*(16.744464-k2 + 4.092*mu*c11))))))))))))
 j0ui = -k1*(e*(q1*(0.744 + mu*c1) + e*(q2*(1.488 + mu*c2) +  &
     e*(q3*(2.232 + mu*c3) + e*(q4*(2.976 + mu*c4) +  &
     e*(q5*(3.720 + mu*c5) + e*(q6*(4.464 + mu*c6) +  &
     e*(q7*(5.208 + mu*c7) + e*(q8*(5.952 + mu*c8) +  &
     e*(q9*(6.696 + mu*c9) + e*(q10*(7.44 + mu*c10)+  &
     e*(q11*(8.184+ mu*c11)))))))))))))
 420 i0ur = .372*e*(r1+e*(2.*r2+e*(3.*r3+e*(4.*r4+e*(5.*r5+e*(6.*r6+  &
     e*(7.*r7+e*(8.*r8+e*(9.*r9+e*(10.*r10+e*11.*r11)))))) ))))
 i0ui = -k1*(e*(r1+e*(r2+e*(r3+e*(r4+e*(r5+e*(r6+e*(r7+e*(r8+e*(r9  &
     +e*(r10+e*r11)))))))))))
 r1   = r1s
 c6   = k1*mu
 c1   = SIN(c6)
 c2   = COS(c6)
 c3   = SQRT(1.+mu*mu)
 c4   = mu/c3
 c5   = c4/(1.+mu*mu)
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 430
   CASE (    2)
     GO TO 440
   CASE (    3)
     GO TO 430
   CASE (    4)
     GO TO 430
   CASE (    5)
     GO TO 440
   CASE (    6)
     GO TO 430
   CASE (    7)
     GO TO 500
   CASE (    8)
     GO TO 500
   CASE (    9)
     GO TO 500
 END SELECT
 430 i1ur = c2*(1.-c4+k1*i0ui) - c1*k1*i0ur
 i1ui =-c2*k1*i0ur - c1*(1.-c4+k1*i0ui)
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 500
   CASE (    2)
     GO TO 440
   CASE (    3)
     GO TO 440
   CASE (    4)
     GO TO 460
   CASE (    5)
     GO TO 440
   CASE (    6)
     GO TO 440
   CASE (    7)
     GO TO 500
   CASE (    8)
     GO TO 500
   CASE (    9)
     GO TO 500
 END SELECT
 440 i2ur3 = c2*(2.*(1.-c4)-c5+k1*i0ui+k2*j0ur)+c1*(c6*(1.-c4)-k1*i0ur  &
     + k2*j0ui)
 i2ui3 = c2*(c6*(1.-c4)-k1*i0ur+k2*j0ui)-c1*(2.*(1.-c4)-c5+k1*i0ui + k2*j0ur)
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 500
   CASE (    2)
     GO TO 500
   CASE (    3)
     GO TO 500
   CASE (    4)
     GO TO 460
   CASE (    5)
     GO TO 450
   CASE (    6)
     GO TO 450
   CASE (    7)
     GO TO 500
   CASE (    8)
     GO TO 500
   CASE (    9)
     GO TO 500
 END SELECT
 450 i2ur3 = 2.0*i20r3 - i2ur3
 IF (ichuz-6 == 0) THEN
   GO TO   460
 ELSE
   GO TO   500
 END IF
 460 car  = 2.*i10r - i1ur
 i1ur = car
 500 dk1r = 0.
 r1   = r1s
 dk1i = 0.
 dk2r = 0.
 dk2i = 0.
 c3   = k1*mu1
 c1   = COS(c3)
 c2   = SIN(c3)
 c3   = m*r1/bigr
 c4   = SQRT(1.+mu1*mu1)
 c5   = kr*x0/br
 c6   = COS(c5)
 c7   = SIN(c5)
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 530
   CASE (    2)
     GO TO 540
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 530
   CASE (    5)
     GO TO 540
   CASE (    6)
     GO TO 530
   CASE (    7)
     GO TO 510
   CASE (    8)
     GO TO 520
   CASE (    9)
     GO TO 510
 END SELECT
 510 i1ur = i10r
 i1ui = i10i
 IF (ichuz-7 == 0) THEN
   GO TO   530
 END IF
 520 i2ur3 = i20r3
 i2ui3 = i20i3
 IF (ichuz-8 == 0) THEN
   GO TO   540
 END IF
 530 ck1r = i1ur + c3*c1/c4
 ck1i = i1ui - c3*c2/c4
 k10  = 1.0 + x0/bigr
 dk1r = ck1r*c6 + ck1i*c7
 dk1i = ck1i*c6 - ck1r*c7
 SELECT CASE ( ichuz )
   CASE (    1)
     GO TO 900
   CASE (    2)
     GO TO 540
   CASE (    3)
     GO TO 540
   CASE (    4)
     GO TO 900
   CASE (    5)
     GO TO 540
   CASE (    6)
     GO TO 540
   CASE (    7)
     GO TO 900
   CASE (    8)
     GO TO 540
   CASE (    9)
     GO TO 540
 END SELECT
 540 c8   = (beta2*(r1/bigr)**2 + (2.+mu1*c3)/(c4*c4))*(-c3/c4)
 c9   = (k1*c3)*( c3/c4)
 ck2r = -i2ur3 + c8*c1 - c9*c2
 ck2i = -i2ui3 - c9*c1 - c8*c2
 k20  = -2.0 - x0*(2.0+beta2*(r1/bigr)**2)/bigr
 dk2r = ck2r*c6 + ck2i*c7
 dk2i = ck2i*c6 - ck2r*c7
 900 CONTINUE
 k1rt1  = t1 *dk1r
 k1it1  = t1 *dk1i
 k2rt2p = t2p*dk2r
 k2it2p = t2p*dk2i
 k10t1  = k10*t1
 k20t2p = k20*t2p
 905 CONTINUE
 kd1r = k1rt1  - k10t1*FLOAT(ind)
 kd1i = k1it1
 kd2r = k2rt2p - k20t2p*FLOAT(ind)
 kd2i = k2it2p
 
 RETURN
END SUBROUTINE tker
