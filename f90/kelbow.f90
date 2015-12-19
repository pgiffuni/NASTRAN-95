SUBROUTINE kelbow
     
!     THIS ROUTINE COMPUTES THE TWO 6 X 6 MATRICES K(NPVT,NPVT) AND
!     K(NPVT,J) FOR A CURVED BAR ELEMENT HAVING END POINTS NUMBERED
!     NPVT AND J
 
!                   ECPT FOR THE ELBOW
 
!     ECPT( 1)  -  IELID         ELEMENT ID. NUMBER
!     ECPT( 2)  -  ISILNO(2)     * SCALAR INDEX NOS. OF THE GRID POINTS
!     ECPT( 3)  -    ...         *
!     ECPT( 4)  -  SMALLV(3)     $ REFERENCE VECTOR
!     ECPT( 5)  -    ...         $
!     ECPT( 6)  -    ...         $
!     ECPT( 7)  -  ICSSV         COOR. SYS. ID FOR SMALLV VECTOR
!     ECPT( 8)  -  IMATID        MATERIAL ID.
!     ECPT( 9)  -  A             CROSS-SECTIONAL AREA
!     ECPT(10)  -  I1            $ AREA MOMENTS OF INERTIA
!     ECPT(11)  -  I2            $
!     ECPT(12)  -  FJ            TORSIONAL CONSTANT
!     ECPT(13)  -  NSM           NON-STRUCTURAL MASS
!     ECPT(14)  -  FE            FORCE ELEM. DESCRIPTIONS, FORCE METHOD
!     ECPT(15)  -  R1            *STRESS RECOVERY COEFFICIENTS
!     ECPT(16)  -  T1            *  RI = RADIAL  LOCATION
!     ECPT(17)  -  R2            *  TI = ANGULAR LOCATION
!     ECPT(18)  -  T2            *       OF STRESS RECOVERY POINTS
!     ECPT(19)  -  R3            *
!     ECPT(20)  -  T3            *
!     ECPT(21)  -  R4            *
!     ECPT(22)  -  T4            *
!     ECPT(23)  -  K1            $  AREA FACTOR FOR SHEAR
!     ECPT(24)  -  K2            $
!     ECPT(25)  -  C             STRESS INTENSIFICATION FACTOR
!     ECPT(26)  -  KX            *  FLEXIBILITY CORRECTION FACTORS
!     ECPT(27)  -  KY            *
!     ECPT(28)  -  KZ            *
!     ECPT(29)  -  R             RADIUS OF CURVATURE
!     ECPT(30)  -  BETAR         ANGLE FROM GA TO GB
!     ECPT(31)  -  MCSIDA        COORD. SYS. ID. FOR GRID POINT A
!     ECPT(32)  -  GPA(3)        *  BASIC COORD. FOR GRID POINT A
!     ECPT(33)  -   ...          *
!     ECPT(34)  -   ...          *
!     ECPT(35)  -  MCSIDB        COORD. SYS. ID. FOR GRID POINT B
!     ECPT(36)  -  GPB(3)        *  BASIC COORD. FOR GRID POINT B
!     ECPT(37)  -   ...          *
!     ECPT(38)  -   ...          *
!     ECPT(39)  -  ELTEMP        AVG. ELEMENT TEMPERATURE
 
!     COMMENTS FROM G.CHAN/UNISYS   7/91
!     ABOUT K1 AND K2, THE AREA FACTORS FOR SHEAR
 
!     THE K1,K2 FOR BAR ARE 0. TO 1.0, AND ARE USED IN K1*G*A AND K2*G*A
!         THE K1,K2 ARE THEREFORE CORRECTION FACTORS FOR STIFFNESS
!     THE K1,K2 ARE USED IN ELBOW IN K1/G*A AND K2/G*A. AND THEREFORE
!         THE K1,K2 ARE COORECTION FACTORS FOR FLEXIBILITY. THE K1,K2
!         IN ELBOW ARE EQUIVALENT TO 1./K1 AND 1./K2 IN BAR ELEMENT.
!         THE PROPER VALUE FOR K1 AND K2 SHOULD BE INFINITY TO 1.0
 
!     IN 1992 COSMIC/NASTRAN, THE USE OF K1 AND K2 IN ELBOW AND BAR
!     ELMENTS ARE SYMCHRONIZED, WITH PROPER VALUES FROM 0. TO 1.0
!     THE K1 AND K2 ARE CHANGED TO 1./K1 AND 1./K2 IN ELBOW ELEMENT
!     SHEAR COMPUTATION. THAT IS, CORRECTION FACTORS FOR STIFFNESS IS
!     USED.
 
!     REFERENCE -  R.J. ROARK: FORMULAS FOR STRESS AND STRAIN,
!     SECTION 35, 'BEAMS FOR RELATIVELY GREAT DEPTH',
!     FOR BEAMS OF SAMLL SPAN/DEPTH RATIO
 
!     K = 1/F = 5/6 FOR RECTANGULAR SECTION
!             = 0.9 FOR SOLID CIRCULAR
!             = 0.5 FOR THIN-WALLED HOOLOW CIRCULAR SECTION
!             = 1.0 CAN BE USED FOR I-BEAM
 
 
 
 LOGICAL :: heat,abasic,bbasic,basic
 REAL :: k1,k2,i1,i2,nsm,kx,ky,kz
 DOUBLE PRECISION :: ta(18),tb(9),smalv0(6),dela,delb,ke,kep,veci,  &
     vecj,veck,fl,fll,df(6,6),determ,h(6,6),dp(16), s(12,12),dampc,kee(12,12)
 DIMENSION        veci(3),vecj(3),veck(3),ecpt(100),iecpt(100),  &
     iz(1),iwork(6,3),f(6,6)
 COMMON /sma1io/  ifcstm,ifmpt,ifdit,idum1,ifecpt,igecpt,ifgpct,  &
     iggpct,ifgei,iggei,ifkgg,igkgg,if4gg,ig4gg,  &
     ifgpst,iggpst,inrw,outrw,clsnrw,clsrw,neor, eor,mcbkgg(7),mcb4gg(7)
 COMMON /zzzzzz/  z(1)
 COMMON /sma1bk/  icstm,ncstm,igpct,ngpct,ipoint,npoint,i6x6k,  &
     n6x6k,i6x64,n6x64
 COMMON /sma1cl/  iopt4,k4ggsw,npvt,left,frowic,lrowic,nrowsc,  &
     tnrows,jmax,nlinks,link(10),idetck,dodet,nogo
 COMMON /sma1ht/  heat
 COMMON /sma1et/  ielid,isilno(2),smallv(3),icssv,imatid,a,i1,i2,  &
     fj,nsm,fe,c1,c2,d1,d2,f1,f2,g1,g2,k1,k2,c,kx,ky,  &
     kz,r,betar,mcsida,gpa(3),mcsidb,gpb(3),tempel
 COMMON /sma1dp/  ke(144),kep(144),dela(6),delb(6)
 COMMON /matin /  matidc,matflg,eltemp,stress,sinth,costh
 COMMON /matout/  e,g,nu,rho,alpha,tsubo,gsube,sigt,sigc,sigs
 COMMON /hmtout/  fk
 COMMON /system/  sysbuf,nout
 EQUIVALENCE      (ielid,ecpt(1),iecpt(1)),(iz(1),z(1)),  &
     (ta(10),tb(1)),(ecpt(71),dp(1)), (kee(1,1),ke(1),s(1,1))
 DATA    dcr   /  .017453292 /
 
 sid(x) = SIN(x*dcr)
 cod(x) = COS(x*dcr)
 dtr(x) = x*dcr
 
!     DETERMINE WHICH POINT IS THE PIVOT POINT.
 
 x    = 1.
 ipvt = 1
 IF (isilno(1) == npvt) GO TO 20
 ipvt = 2
 IF (isilno(2) /= npvt) CALL mesage (-30,34,iecpt(1))
 
!     SET UP POINTERS TO COORD. SYSTEM IDS
 
 20 jcsida = 31
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
 fld   = DSQRT(dp(7)**2 + dp(8)**2 + dp(9)**2)
 IF (fld <= 0.000) GO TO 1010
 dp(7) = dp(7)/fld
 dp(8) = dp(8)/fld
 dp(9) = dp(9)/fld
 
!     DETERMINE IF POINT A AND B ARE IN BASIC COORDINATES
 
 abasic =.true.
 bbasic =.true.
 IF (icsida /= 0) abasic =.false.
 IF (icsidb /= 0) bbasic =.false.
 
!     COMPUTE THE TRANSFORMATION MATRICES TA AND TB IF NECESSARY
 
 IF (abasic) GO TO 30
 CALL transd (ecpt(jcsida),ta)
 CALL gmmatd (ta,3,3,0, dp(7),3,1,0, vecj)
 CALL gmmatd (ta,3,3,0, dp(1),3,1,0, dp(14))
 dp(1) = dp(14)
 dp(2) = dp(15)
 dp(3) = dp(16)
 GO TO 35
 30 CONTINUE
 vecj(1) = dp(7)
 vecj(2) = dp(8)
 vecj(3) = dp(9)
 35 IF (bbasic) GO TO 40
 CALL transd (ecpt(jcsidb),tb)
 CALL gmmatd (tb,3,3,0, dp(4),3,1,0, dp(14))
 dp(4) = dp(14)
 dp(5) = dp(15)
 dp(6) = dp(16)
 40 CONTINUE
 
!     CALCULATE TRUE LENGTH OF ELBOW
 
 fl = DBLE(r*dtr(betar))
 IF (fl == 0.0D0) GO TO 1010
 
!     NOW THAT LENGTH HAS BEEN COMPUTED, BRANCH IF THIS IS A -HEAT-
!     FORMULATION.
 
 IF (heat) GO TO 2000
 
!     CONSTRUCT VECTOR FROM A TO B
 
 smalv0(1) = dp(4) - dp(1)
 smalv0(2) = dp(5) - dp(2)
 smalv0(3) = dp(6) - dp(3)
 fll = DSQRT(smalv0(1)**2 + smalv0(2)**2 + smalv0(3)**2)
 IF (fll == 0.0D0) GO TO 1010
 smalv0(1) = smalv0(1)/fll
 smalv0(2) = smalv0(2)/fll
 smalv0(3) = smalv0(3)/fll
 
!     COMPUTE THE K VECTOR VECK = SMALV0 X VECJ
 
 veck(1) = smalv0(2)*vecj(3) - smalv0(3)*vecj(2)
 veck(2) = smalv0(3)*vecj(1) - smalv0(1)*vecj(3)
 veck(3) = smalv0(1)*vecj(2) - smalv0(2)*vecj(1)
 fll = DSQRT(veck(1)**2 + veck(2)**2 + veck(3)**2)
 IF (fll == 0.0D0) GO TO 1010
 veck(1) = veck(1)/fll
 veck(2) = veck(2)/fll
 veck(3) = veck(3)/fll
 
!     COMPUTE THE I VECTOR  VECI = VECJ X VECK
 
 veci(1) = vecj(2)*veck(3) - vecj(3)*veck(2)
 veci(2) = vecj(3)*veck(1) - vecj(1)*veck(3)
 veci(3) = vecj(1)*veck(2) - vecj(2)*veck(1)
 fll = DSQRT(veci(1)**2 + veci(2)**2 + veci(3)**2)
 IF (fll == 0.0D0) GO TO 1010
 veci(1) = veci(1)/fll
 veci(2) = veci(2)/fll
 veci(3) = veci(3)/fll
 
!     SEARCH THE MATERIAL PROPERTIES TABLE FOR E,G AND THE DAMPING
!     CONSTANT.
 
 matidc = imatid
 matflg = 1
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
 
!     AREA FACTORS FOR SHEAR ARE FROM NEAR ZERO TO ONE
 
 IF (k1 < 1.0E-8) k1 = 1.0
 IF (k2 < 1.0E-8) k2 = 1.0
 IF (k1 > 1.0) k1 = 1.0/k1
 IF (k2 > 1.0) k2 = 1.0/k2
 
!     THE FOLLOWING CODE WAS TAKEN FROM SAP4 BENDKS ROUTINE
!     FOR A CURVED PIPE ELEMENT
 
!     COMPUTE SECTION PROPERTY CONSTANTS
 
 t   = dtr(betar)
 ra  = r/(a*e)
 rv1 = r/(2.*k1*g*a)
 rv2 = k1/k2*rv1
 rt  = r/(g*fjk*2.)
 rb0 = r/(e*fi2*2.)
 rb1 = r/(e*fi1)
 r2  = r**2
 
!     COMPUTE COMMON TRIGONOMETRIC CONSTANTS
 
 st  = sid(betar)
 ct  = cod(betar)
 s2t = sid(2.0*betar)
 c2t = cod(2.0*betar)
 
!     FORM THE NODE FLEXIBILITY MATRIX AT NODE J REFERENCED TO THE
!     LOCAL (X,Y,Z) COORDINATE SYSTEM AT NODE I.
 
!     X - DIRECTION...  IN-PLANE TANGENT TO THE BEND AT NODE I AND
!                       DIRECTED TOWARD NODE J
!     Y - DIRECTION...  IN-PLANE AND DIRECTED RADIALLY INWARD TO THE
!                       CENTER OF CURVATURE
!     Z - DIRECTION...  OUT OF PLANE AND ORTHOGONAL TO X AND Y
 
 DO  i = 1,6
   DO  k = i,6
     f(i,k)  = 0.0
   END DO
 END DO
 
!     A X I A L
 
 f(1,1) = f(1,1) + 0.25*ra*(2.0*t+s2t)
 f(2,2) = f(2,2) + 0.25*ra*(2.0*t-s2t)
 
!     N O T E   (COEFFICIENT CHANGE)
 
 f(1,2) = f(1,2) + 0.50*ra*st**2
 
!     S H E A R
 
 f(1,1) = f(1,1) + 0.5*rv1*(2.0*t-s2t)
 f(2,2) = f(2,2) + 0.5*rv1*(2.0*t+s2t)
 f(3,3) = f(3,3) + 2.0*rv2* t
 
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
     df(k,i) = DBLE(f(i,k))
     df(i,k) = df(k,i)
   END DO
 END DO
 
!     WRITE (6,4005) DF
 
!     INVERT FLEX TO FORM STIFFNESS
 
 CALL inverd (6,df,6,dum,0,determ,ising,iwork)
 IF (ising == 2) WRITE (6,4002) f
 IF (ising == 2) CALL mesage (-30,38,ecpt(1))
 4002 FORMAT (1X,34HELBOW stiffness matrix is singular, /,(5X,6E13.5))
 
 
!     SET UP THE FORCE TRANSFORMATION RELATING REACTIONS AT NODE I
!     ACTING ON THE MEMBER END DUE TO UNIT LOADS APPLIED TO THE MEMBER
!     END AT NODE J.
 
 DO  i = 1,6
   DO  k = 1,6
     h(i,k) = 0.0D0
   END DO
 END DO
 
 DO  k = 1,6
   h(k,k) =-1.0D0
 END DO
 
 h(4,3) =-DBLE(r*(1.0-ct))
 h(5,3) = DBLE(r*st)
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
     s(ir,ic+6) = 0.0D0
     DO  in = 1,6
       s(ir,ic+6) = s(ir,ic+6) + h(ir,in)*df(in,ic)
     END DO
   END DO
 END DO
 
 DO  ir = 1,6
   DO  ic = ir,6
     s(ir,ic)  = 0.0D0
     DO  in = 1,6
       s(ir,ic) = s(ir,ic) + s(ir,in+6)*h(ic,in)
     END DO
   END DO
 END DO
 
!     REFLECT FOR SYMMETRY
 
 DO  i = 1,12
   DO  k = i,12
     s(k,i)   = s(i,k)
   END DO
 END DO
 
 j = 0
 IF (ipvt == 2) GO TO 327
 ilow = 1
 ilim = 72
 GO TO 329
 327 ilow = 73
 ilim = 144
 329 DO  i = ilow,ilim,12
   low = i
   lim = low + 5
   DO  k = low,lim
     j = j + 1
     kep(j) = ke(k)
     kep(j+36) = ke(k+6)
   END DO
 END DO
 
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
 
!     ZERO OUT THE ARRAY WHERE THE 3 X 3 MATRIX H AND THE W  AND W
!     6 X 6  MATRICES WILL RESIDE.                         A      B
 
 DO  i = 28,108
   ke(i) = 0.0D0
 END DO
 ipass = 1
 iwbeg = 0
 
!     SET UP POINTERS
 
 IF (ipvt-1 == 0) THEN
   GO TO   360
 ELSE
   GO TO   365
 END IF
 360 basic  = abasic
 jcsid  = jcsida
 ikel   = 1
 INDEX  = isilno(1)
 GO TO 368
 365 basic  = bbasic
 jcsid  = jcsidb
 ikel   = 37
 INDEX  = isilno(2)
 
!     SET UP THE -G- MATRIX.  IG POINTS TO THE BEGINNING OF THE G
!     MATRIX. G = AT X TI
 
 368 ig = 1
 IF (basic) GO TO 380
 CALL transd (ecpt(jcsid),ke(10))
 CALL gmmatd (ke(1),3,3,0, ke(10),3,3,0, ke(19))
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
 
!                       T      E
!     FORM THE PRODUCT W   X  K   AND STORE IN KEP(73)
!                       NPVT
 
 CALL gmmatd (ke(37),6,6,1, kep(ikel),6,6,0, kep(73))
 
!     COMPUTE THE FINAL ANSWER AND STORE IN KEP(109)
 
 CALL gmmatd (kep(73),6,6,0, ke(iwbeg+37),6,6,0, kep(109))
 
!     INSERT THIS 6 X 6
 
 CALL sma1b (kep(109),INDEX,-1,ifkgg,0.0D0)
 IF (iopt4 == 0 .OR. gsube == 0.0) GO TO 400
 k4ggsw = 1
 CALL sma1b (kep(109),INDEX,-1,if4gg,dampc)
 
!     IF IPASS = 2, WE ARE DONE.  OTHERWISE COMPUTE THE OFF-DIAGONAL
!     6 X 6.
 
 400 IF (ipass == 2) GO TO 500
 iwbeg = 36
 ipass = 2
 DO  i = 28,36
   ke(i) = 0.0D0
 END DO
 IF (ipvt-1 == 0) THEN
   GO TO   365
 ELSE
   GO TO   360
 END IF
 500 RETURN
 
 1010 CALL mesage (30,26,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
 
!     HEAT FORMULATION CONTINUES HERE.  GET MATERIAL PROPERTY -K- FROM
!     HMAT
 
 2000 matflg = 1
 matidc = iecpt( 8)
 eltemp = ecpt(39)
 CALL hmat (ielid)
 
 fl = DBLE(fk)*DBLE(ecpt(9))/(dp(9)*dp(10)*DBLE(dcr))
 IF (npvt == iecpt(3)) fl = -fl
 DO  i = 1,2
   CALL sma1b (fl,iecpt(i+1),npvt,ifkgg,0.0D0)
   fl = -fl
 END DO
 RETURN
END SUBROUTINE kelbow
