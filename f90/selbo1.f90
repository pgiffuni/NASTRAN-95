SUBROUTINE selbo1
     
!     THIS ROUTINE IS PHASE 1 OF STRESS DATA RECOVERY FOR THE ELBOW
!     ELEMENT MUCH OF THE CODE WAS LIFTED FROM THE KELBOW SUBROUTINE
 
!     ECPT FOR THE ELBOW
 
!     ECPT( 1)  -  IELID          ELEMENT ID. NUMBER
!     ECPT( 2)  -  ISILNO(2)      * SCALAR INDEX NOS. OF THE GRID POINTS
!     ECPT( 3)  -    ...          *
!     ECPT( 4)  -  SMALLV(3)      $ REFERENCE VECTOR
!     ECPT( 5)  -    ...          $
!     ECPT( 6)  -    ...          $
!     ECPT( 7)  -  ICSSV          COOR. SYS. ID FOR SMALLV VECTOR
!     ECPT( 8)  -  IMATID         MATERIAL ID.
!     ECPT( 9)  -  A              CROSS-SECTIONAL AREA
!     ECPT(10)  -  I1             $ AREA MOMENTS OF INERTIA
!     ECPT(11)  -  I2             $
!     ECPT(12)  -  FJ             TORSIONAL CONSTANT
!     ECPT(13)  -  NSM            NON-STRUCTURAL MASS
!     ECPT(14)  -  FE             FORCE ELEM. DESCRIPTIONS, FORCE METHOD
!     ECPT(15)  -  R1             *STRESS RECOVERY COEFFICIENTS
!     ECPT(16)  -  T1             *  RI=RADIAL LOCATION
!     ECPT(17)  -  R2             *  TI=ANGULAR LOCATION
!     ECPT(18)  -  T2             *     OF STRESS RECOVERY POINTS
!     ECPT(19)  -  R3             *
!     ECPT(20)  -  T3             *
!     ECPT(21)  -  R4             *
!     ECPT(22)  -  T4             *
!     ECPT(23)  -  K1             $ AREA FACTOR FOR SHEAR
!     ECPT(24)  -  K2             $
!     ECPT(25)  -  C              STRESS INTENSIFICATION FACTOR
!     ECPT(26)  -  KX             * FLEXIBILITY CORRECTION FACTORS
!     ECPT(27)  -  KY             *
!     ECPT(28)  -  KZ             *
!     ECPT(29)  -  R              RADIUS OF CURVATURE
!     ECPT(30)  -  BETAR          ANGLE FROM GA TO GB
!     ECPT(31)  -  MCSIDA         COORD. SYS. ID. FOR GRID POINT A
!     ECPT(32)  -  GPA(3)         *BASIC COORD. FOR GRID POINT A
!     ECPT(33)  -   ...           *
!     ECPT(34)  -   ...           *
!     ECPT(35)  -  MCSIDB         COORD. SYS. ID. FOR GRID POINT B
!     ECPT(36)  -  GPB(3)         *BASIC COORD. FOR GRID POINT B
!     ECPT(37)  -   ...           *
!     ECPT(38)  -   ...           *
!     ECPT(39)  -  ELTEMP         AVG. ELEMENT TEMPERATURE
 
 
 LOGICAL :: abasic,bbasic,basic
 REAL :: l,i1,i2,k1,k2,ke,kep,nsm,hut( 6),kee(12,12), kx,ky,kz
 DIMENSION       veci(3),vecj(3),veck(3),ecpt(100),iecpt(100),  &
     ta(18),tb(9),smalv0(6),dp(20),f(6,6),s(12,12), h(6,6),df(6,6)
 COMMON /sdr2x5/ ielid,isilno(2),smallv(3),icssv,imatid,a,i1,i2,  &
     fj,nsm,fe,c1,c2,d1,d2,f1,f2,g1,g2,k1,k2,c,  &
     kx,ky,kz,r,betar,mcsida,gpa(3),mcsidb,gpb(3), tempel,dum3(61)
 COMMON /sdr2x5/ jelid,jsilno(2),sa(36),sb(36),out(21),therm(30)
 COMMON /sdr2x6/ ke(144),kep(144),dela(6),delb(6)
 COMMON /matin / matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/ e,g,nu,rho,alpha,t_sub_0,g_sub_e,sigt,sigc,sigs
 EQUIVALENCE     (ielid,ecpt(1),iecpt(1)), (ta(10),tb(1)),  &
     (kee(1,1),ke(1),s(1,1))
 DATA    dcr   / .017453292 /
 
 sid(x) = SIN(x*dcr)
 cod(x) = COS(x*dcr)
 dtr(x) = x*dcr
 
 x    = 1.0
 isop = -1
 
!     SET UP POINTERS TO COORD. SYSTEM IDS
 
 jcsida = 31
 jcsidb = 35
 icsida = iecpt(31)
 icsidb = iecpt(35)
 
!     DEFINE LOCATION OF END A, END B IN TERMS OF DP(1) THRU DP(6)
 
 dp(1) = ecpt(jcsida+1)
 dp(2) = ecpt(jcsida+2)
 dp(3) = ecpt(jcsida+3)
 dp(4) = ecpt(jcsidb+1)
 dp(5) = ecpt(jcsidb+2)
 dp(6) = ecpt(jcsidb+3)
 
!     DEFINE COMPONENTS OF VECTOR FROM END A TO CENTER OF CURVATURE,C
 
 dp(7) = ecpt(4)
 dp(8) = ecpt(5)
 dp(9) = ecpt(6)
 fld   = SQRT(dp(7)**2 + dp(8)**2 + dp(9)**2)
 dp(7) = dp(7)/fld
 dp(8) = dp(8)/fld
 dp(9) = dp(9)/fld
 
!     DETERMINE IF POINT A AND B ARE IN BASIC COORDINATES
 
 abasic =.true.
 bbasic =.true.
 IF (icsida /= 0) abasic =.false.
 IF (icsidb /= 0) bbasic =.false.
 
!     COMPUTE THE TRANSFORMATION MATRICES TA AND TB IF NECESSARY
 
 IF (abasic) GO TO 60
 CALL transs (ecpt(jcsida),ta)
 CALL gmmats (ta,3,3,0, dp(7),3,1,0, vecj)
 CALL gmmats (ta,3,3,0, dp(1),3,1,0, dp(14))
 dp(1) = dp(14)
 dp(2) = dp(15)
 dp(3) = dp(16)
 GO TO 61
 60 CONTINUE
 vecj(1) = dp(7)
 vecj(2) = dp(8)
 vecj(3) = dp(9)
 61 IF (bbasic) GO TO 62
 CALL transs (ecpt(jcsidb),tb)
 CALL gmmats (tb,3,3,0, dp(4),3,1,0, dp(14))
 dp(4) = dp(14)
 dp(5) = dp(15)
 dp(6) = dp(16)
 62 CONTINUE
 
!     CONSTRUCT VECTOR FROM A TO B
 
 smalv0(1) = dp(4) - dp(1)
 smalv0(2) = dp(5) - dp(2)
 smalv0(3) = dp(6) - dp(3)
 fll       = SQRT(smalv0(1)**2 + smalv0(2)**2 + smalv0(3)**2)
 smalv0(1) = smalv0(1)/fll
 smalv0(2) = smalv0(2)/fll
 smalv0(3) = smalv0(3)/fll
 
!     COMPUTE THE K VECTOR VECK = SMALV0 X VECJ
 
 veck(1) = smalv0(2)*vecj(3) - smalv0(3)*vecj(2)
 veck(2) = smalv0(3)*vecj(1) - smalv0(1)*vecj(3)
 veck(3) = smalv0(1)*vecj(2) - smalv0(2)*vecj(1)
 fll     = SQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 veck(1) = veck(1)/fll
 veck(2) = veck(2)/fll
 veck(3) = veck(3)/fll
 
!     COMPUTE THE I VECTOR  VECI = VECJ X VECK
 
 veci(1) = vecj(2)*veck(3) - vecj(3)*veck(2)
 veci(2) = vecj(3)*veck(1) - vecj(1)*veck(3)
 veci(3) = vecj(1)*veck(2) - vecj(2)*veck(1)
 fll     = SQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 veci(1) = veci(1)/fll
 veci(2) = veci(2)/fll
 veci(3) = veci(3)/fll
 
!     SEARCH THE MATERIAL PROPERTIES TABLE FOR E,G AND THE DAMPING
!     CONSTANT.
 
 matidc = imatid
 matflg = 1
 IF (isop == 3) matflg = 12
 eltemp = tempel
 CALL mat (iecpt(1))
 dampc = g_sub_e
 
!     SET UP INTERMEDIATE VARIABLES FOR ELEMENT STIFFNESS MATRIX
!     CALCULATION
 
 IF (kx < 1.0E-8) kx = 1.0
 IF (ky < 1.0E-8) ky = 1.0
 IF (kz < 1.0E-8) kz = 1.0
 fi1 = i1/kz
 fi2 = i2/ky
 fjk = fj/kx
 
 
!     THE FOLLOWING CODE WAS TAKEN FROM SAP4 BENDKS ROUTINE FOR A CURVED
!     PIPE ELEMENT
 
 
!     COMPUTE SECTION PROPERTY CONSTANTS
 
 t  = dtr(betar)
 ra = r/(a*e)
 rv1= k1*r/(2.*g*a)
 rv2= k2/k1*rv1
 rt = r/(g*fjk*2.)
 rb0= r/(e*fi2*2.)
 rb1= r/(e*fi1)
 r2 = r**2
 
!     COMPUTE COMMON TRIGONOMETRIC CONSTANTS
 
 st  = sid(betar)
 ct  = cod(betar)
 s2t = sid(2.0*betar)
 c2t = cod(2.0*betar)
 
!     FORM THE NODE FLEXIBILITY MATRIX AT NODE J REFERENCED TO THE
!     LOCAL (X,Y,Z) COORDINATE SYSTEM AT NODE I.
 
!     X - DIRECTION  IN-PLANE TANGENT TO THE BEND AT NODE I AND
!                    DIRECTED TOWARD NODE J
!     Y - DIRECTION  IN-PLANE AND DIRECTED RADIALLY INWARD TO THE
!                    CENTER OF CURVATURE
!     Z - DIRECTION  OUT OF PLANE AND ORTHOGONAL TO X AND Y
 
 DO  i = 1,6
   DO  k = i,6
     f(i,k) = 0.0
   END DO
 END DO
 
!     A X I A L
 
 f(1,1) = f(1,1) + 0.25*ra*(2.0*t + s2t)
 f(2,2) = f(2,2) + 0.25*ra*(2.0*t - s2t)
 
!     N O T E   (COEFFICIENT CHANGE)
 
 f(1,2) = f(1,2) + 0.50*ra*st**2
 
!     S H E A R
 
 f(1,1) = f(1,1) + 0.5*rv1*(2.0*t - s2t)
 f(2,2) = f(2,2) + 0.5*rv1*(2.0*t + s2t)
 f(3,3) = f(3,3) + 2.0*rv2*t
 
!     N O T E   (SIGN CHANGE)
 
 f(1,2) = f(1,2) - rv1*st**2
 
!     T O R S I O N
 
 f(3,3) = f(3,3) + 0.5*rt*r2*(6.0*t+s2t-8.0*st)
 f(4,4) = f(4,4) + 0.5*rt*   (2.0*t+s2t)
 f(5,5) = f(5,5) + 0.5*rt*   (2.0*t-s2t)
 f(3,4) = f(3,4) +     rt*r *(st-t*ct)
 f(3,5) = f(3,5) +     rt*r *(2.0-2.0*ct-t*st)
 f(4,5) = f(4,5) + 0.5*rt*   (1.0-c2t)
 
!     B E N D I N G
 
 f(1,1) = f(1,1) + 0.25*rb1*r2*(2.0*t*(2.0+c2t)-3.0*s2t)
 f(2,2) = f(2,2) + 0.25*rb1*r2*(2.0*t*(2.0-c2t)+3.0*s2t-8.0*st)
 f(3,3) = f(3,3) + 0.50*rb0*r2*(2.0*t-s2t)
 f(4,4) = f(4,4) + 0.50*rb0*   (2.0*t-s2t)
 f(5,5) = f(5,5) + 0.50*rb0*   (2.0*t+s2t)
 f(6,6) = f(6,6) +      rb1*t
 f(1,2) = f(1,2) + 0.25*rb1*r2*(1.0+3.0*c2t+2.0*t*s2t-4.0*ct)
 f(1,6) = f(1,6) -      rb1*r *(st-t*ct)
 f(2,6) = f(2,6) +      rb1*r *(t*st+ct-1.0)
 f(3,4) = f(3,4) +      rb0*r *(st-t*ct)
 f(3,5) = f(3,5) -      rb0*r *t*st
 f(4,5) = f(4,5) - 0.50*rb0*   (1.0-c2t)
 
 
!     FORM SYMMETRICAL UPPER PART OF FLEX MATRIX
 
 DO  i = 1,6
   DO  k = i,6
     df(k,i) = f(i,k)
     df(i,k) = df(k,i)
   END DO
 END DO
 
 
!     INVERT FLEX TO FORM STIFFNESS
 
 CALL invers (6,df,6,dum,0,determ,ising,h)
 IF (ising == 2) WRITE (6,4002) f
 IF (ising == 2) CALL mesage (-30,38,ecpt(1))
 4002 FORMAT (35H elbow stiffness matrix is singular, /,(5X,6E13.5))
 
 
!     SET UP THE FORCE TRANSFORMATION RELATING REACTIONS AT NODE I
!     ACTING ON THE MEMBER END DUE TO UNIT LOADS APPLIED TO THE MEMBER
!     END AT NODE J.
 
 DO  i = 1,6
   DO  k = 1,6
     h(i,k) = 0.0
   END DO
 END DO
 
 DO  k = 1,6
   h(k,k) =-1.0
 END DO
 
 h(4,3) =-(r*(1.0 - ct))
 h(5,3) = (r*st)
 h(6,1) =-h(4,3)
 h(6,2) =-h(5,3)
 
!     FORM THE UPPER TRIANGULAR PORTION OF THE LOCAL ELEMENT STIFFNESS
!     MATRIX FOR THE BEND
 
 DO  k = 1,6
   DO  i = k,6
     s(k+6,i+6) = df(k,i)
   END DO
 END DO
 
 DO  ir = 1,6
   DO  ic = 1,6
     s(ir,ic+6) = 0.0
     DO  in = 1,6
       s(ir,ic+6) = s(ir,ic+6) + h(ir,in)*df(in,ic)
     END DO
   END DO
 END DO
 
 DO  ir = 1,6
   DO  ic = ir,6
     s(ir,ic)  = 0.0
     DO  in = 1,6
       s(ir,ic) = s(ir,ic) + s(ir,in+6)*h(ic,in)
     END DO
   END DO
 END DO
 
!     REFLECT FOR SYMMETRY
 
 DO  i = 1,12
   DO  k = i,12
     s(k,i) = s(i,k)
   END DO
 END DO
 
!            E
!     STORE K   IN KEP(1) THRU KEP(36) AND
!            AA
 
!            E
!     STORE K   IN KEP(37) THRU KEP(72)
!            AB
 
 j = 0
 DO  i = 1,72,12
   low = i
   lim = low + 5
   DO  k = low,lim
     j = j + 1
     kep(j) = ke(k)
     kep(j+36) = ke(k+6)
   END DO
 END DO
 
!     COMPUTE THERMAL MATRIX
 
 l = dcr*ecpt(29)*ecpt(30)
 DO  i = 1,6
   hut(i) = 0.0
 END DO
 alphar = alpha*r
 hut(1) =-alphar*sid(betar)
 hut(2) =-alphar*(1.-cod(betar))
 hut(6) = 0.0
 CALL gmmats (kep(1),6,6,0, hut,6,1,0, therm(1))
 
!                                                             T
!     STORE VECI, VECJ, VECK IN KE(1) THRU KE(9) FORMING THE A  MATRIX.
 
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
 
 basic = abasic
 jcsid = jcsida
 iwbeg = 0
 ikel  = 1
 iab   = 1
 INDEX = isilno(1)
 
!     ZERO OUT THE ARRAY WHERE THE 3 X 3 MATRIX AND THE W  AND W  6 X 6
!     MATRICES WILL RESIDE.                              A      B
 
 DO  i = 28,108
   ke(i) = 0.0
 END DO
 
!     SET UP THE -G- MATRIX. IG POINTS TO THE BEGINNING OF THE G MATRIX.
!     G = AT X TI
 
 360 ig = 1
 IF (basic) GO TO 380
 CALL transs (ecpt(jcsid),ke(10))
 CALL gmmats (ke(1),3,3,0, ke(10),3,3,0, ke(19))
 ig = 19
 
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
 
!                              E                    E
!     FORM THE PRODUCT  S  =  K   X  W   OR  S   = K    X  W, DEPENDING
!                        A     AA     A       B     AB      B
!     UPON WHICH POINT WE ARE WORKING WITH.
 
 CALL gmmats (kep(ikel),6,6,0, ke(iwbeg+37),6,6,0, sa(iab))
 
!     IF THE POINT UNDER CONSIDERATION IS POINT B WE ARE FINISHED. IF
!     NOT, SET UP POINTS AND INDICATORS FOR WORKING WITH POINT B.
 
 IF (iwbeg == 36) GO TO 500
 basic = bbasic
 jcsid = jcsidb
 iwbeg = 36
 ikel  = 37
 iab   = 37
 INDEX = isilno(2)
 DO  i = 28,36
   ke(i) = 0.0
 END DO
 GO TO 360
 
!     FILL REMAINDER OF OUTPUT BLOCK.
 
 500 jelid   = ielid
 jsilno(1) = isilno(1)
 jsilno(2) = isilno(2)
 i12     = 0.
 out( 1) = a*e*alpha
 out( 2) = a*e/l
 out( 3) = a
 out( 4) = fj
 out( 5) = i1
 out( 6) = i2
 out( 7) = c
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
 out(20) = r
 out(21) = betar
 RETURN
END SUBROUTINE selbo1
