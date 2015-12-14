SUBROUTINE fwmw (nd,NE,sgs,cgs,irb,a0,arb,xble,xbte,yb,zb,xs,  &
        ys,zs,nas,nasb,kr,beta2,cbar,avr,fwz,fwy)
     
!     CALCULATES THE EFFECT OF A DOUBLET PLUS ANY CONTRIBUTIONS DUE TO
!     IMAGES, SYMMETRY AND GROUND EFFECT ON BODY
 
 
 INTEGER, INTENT(IN)                      :: nd
 INTEGER, INTENT(IN)                      :: NE
 REAL, INTENT(IN)                         :: sgs
 REAL, INTENT(IN)                         :: cgs
 INTEGER, INTENT(IN OUT)                  :: irb
 REAL, INTENT(IN)                         :: a0
 REAL, INTENT(IN)                         :: arb(1)
 REAL, INTENT(IN OUT)                     :: xble
 REAL, INTENT(IN OUT)                     :: xbte
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN OUT)                     :: xs
 REAL, INTENT(IN)                         :: ys
 REAL, INTENT(IN)                         :: zs
 INTEGER, INTENT(IN)                      :: nas
 INTEGER, INTENT(IN)                      :: nasb(1)
 INTEGER, INTENT(IN OUT)                  :: kr
 REAL, INTENT(IN OUT)                     :: beta2
 REAL, INTENT(IN OUT)                     :: cbar
 REAL, INTENT(IN)                         :: avr(1)
 COMPLEX, INTENT(OUT)                     :: fwz
 COMPLEX, INTENT(OUT)                     :: fwy
 
 
 
!     ND        SYMMETRY FLAG
!     NE        GROUND EFFECTS FLAG
!     SGS       SINE   OF SENDING POINT DIHEDRAL ANGLE
!     CGS       COSINE OF SENDING POINT DIHEDRAL ANGLE
!     IRB       NUMBER OF THE RECEIVING BODY
!     A0        RADIUS OF THE BODY
!     ARB       ARRAY OF RATIOS OF BODY AXIS
!     XBLE      LEADING  EDGE COORDINATE OF SLENDER BODY ELEMENT
!     XBTE      TRAILING EDGE COORDINATE OF SLENDER BODY ELEMENT
!     YB        ARRAY CONTAINING THE Y-COORDINATES OF THE BODIES
!     ZB        ARRAY CONTAINING THE Y-COORDINATES OF THE BODIES
!     XS        1/4-CHORD  X-COORDINATE  OF SLENDER BODY ELEMENT
!     YS        Y-COORDINATE OF SENDING POINT
!     ZS        Z COORDINATE OF THE SENDING POINT
!     NAS       NUMBER OF ASSOCIATED BODIES
!     NASB      ARRAY CONTAINING THE ASSOCIATED BODY NOS.
!     KR        REDUCED FREQUENCY
!     BETA2     = 1 - MACH**2
!     CBAR      REFERENCE CHARD LENGTH
!     AVR       ARRAY OF BODY RADII
!     FWZ       OUTPUT Z-FORCE
!     FWY       OUTPUT Y FORCE
 
 fwz  = CMPLX(0.0,0.0)
 fwy  = CMPLX(0.0,0.0)
 
 dmmy = 0.0
 infl = 1
 
!     ARG-R  ARGUMENTS
 
 dyb  = yb(irb)
 dzb  = zb(irb)
 da   = a0
 deleps = 1.0
 c    = cgs
 s    =-sgs
 dy   = ys
 dz   = zs
 itype= 1
 k    = 1
 ASSIGN 100 TO iret1
 GO TO 2000
 100 sy   = 1.0
 sz   = 1.0
 sg   = sgs
 ASSIGN 200 TO iret1
 GO TO 5000
 200 CONTINUE
 
!     CHECK SYMMETRY FLAG. BRANCH IF EQUAL TO ZERO
 
 IF (nd == 0) GO TO 700
 
!     PORTION FOR SYMMETRIC CALCULATIONS
 
 deleps = nd
 c    = cgs
 s    = sgs
 dy   =-ys
 dz   = zs
 itype= 1
 k    = 2
 ASSIGN 300 TO iret1
 GO TO 2000
 300 CONTINUE
 sy   =-1.0
 sz   = 1.0
 sg   =-sgs
 ASSIGN 400 TO iret1
 GO TO 5000
 400 CONTINUE
 
!     CHECK GROUND EFFECTS FLAG. SKIP IF ZERO
 
 IF (NE == 0) GO TO 7000
 
!     PORTION FOR COMBINATION OF SYMMETRY AND GROUND EFFECTS
 
 itype = 1
 k     = 3
 deleps= nd*NE
 c     = cgs
 s     =-sgs
 dy    =-ys
 dz    =-zs
 ASSIGN 500 TO iret1
 GO TO 2000
 500 CONTINUE
 sy    =-1.0
 sg    = sgs
 sz    =-1.0
 ASSIGN 600 TO iret1
 GO TO 5000
 600 CONTINUE
 GO TO 800
 
!     SKIP GROUND EFFECTS CALCULATIONS IF FLAG IS ZERO
 
 700 IF (NE == 0) GO TO 7000
 
!     PORTION FOR GROUND EFFECTS ONLY
 
 800 CONTINUE
 deleps = NE
 dy   = ys
 dz   =-zs
 c    = cgs
 s    = sgs
 itype= 1
 k    = 4
 ASSIGN 900 TO iret1
 GO TO 2000
 900 CONTINUE
 sy   = 1.0
 sz   =-1.0
 sg   =-sgs
 ASSIGN 1000 TO iret1
 GO TO 5000
 1000 CONTINUE
 RETURN
 
!     CALCULATION OF EFFECTIVE FORCES
 
 2000 CONTINUE
 rho2 = (dy-dyb)**2 + (dz-dzb)**2
 rho  = SQRT(rho2)
 b    = avr(irb)*arb(irb)
 rhodb= rho/b
 f    = 1.0
 IF (rho <= b) GO TO 2020
 f    = rhodb/(arb(irb)*(rhodb-1.0)+1.0)
 2020  CONTINUE
 zbar = (dz-dzb)/(f*arb(irb)) + dzb
 CALL fzy2 (xs,xble,xbte,dy,zbar,dyb,dzb,da,beta2,cbar,kr,dfzzr,  &
     dfzzi,dfzyr,dfzyi,dfyzr,dfyzi,dfyyr,dfyyi)
 
 fwzr = c*dfzzr + s*dfzyr
 fwzi = c*dfzzi + s*dfzyi
 fwz  = fwz + deleps*CMPLX(fwzr,fwzi)
 fwyr = c*dfyzr + s*dfyyr
 fwyi = c*dfyzi + s*dfyyi
 fwy  = fwy + deleps*CMPLX(fwyr,fwyi)
 2060 SELECT CASE ( itype )
   CASE (    1)
     GO TO 3000
   CASE (    2)
     GO TO 6000
 END SELECT
 3000 GO TO iret1, (100,200,300,400,500,600,800,900,1000)
 
!     CALCULATION LOOP FOR ASSOCIATED BODIES
 
 5000 IF (nas <= 0) GO TO 3000
 i = 1
 itype = 2
 5100 ib = nasb(i)
 
!     CHECK TO SEE IF THE ASSOCIATED BODY IS THE RECEIVING BODY.
 
 IF (ib /= irb) GO TO 5800
 
!     IF IT IS DETERMINE IF THE SENDING POINT IS OUTSIDE OR INSIDE THE
!     BODY.
 
 SELECT CASE ( k )
   CASE (    1)
     GO TO 5600
   CASE (    2)
     GO TO 5500
   CASE (    3)
     GO TO 5400
   CASE (    4)
     GO TO 5300
 END SELECT
 5300 IF (dyb /= 0.0) GO TO 5800
 GO TO 5800
 5400 IF (dyb /= 0.0) GO TO 5800
 5500 IF (dzb /= 0.0) GO TO 5800
 5600 CONTINUE
 5800 CONTINUE
 eta   = sy*ys
 zeta  = sz*zs
 zbi   = sz*zb(ib)
 ybi   = sy*yb(ib)
 darib = arb(ib)
 daib  = avr(ib)
 CALL subi (daib,zbi,ybi,darib,eta,zeta,cgs,sg,dmmy,dmmy,dmmy,dy,  &
     dz,dmmy,dmmy,dmmy,dmmy,s,c,infl,ioutfl)
 IF (ioutfl == 0) THEN
   GO TO  2060
 ELSE
   GO TO  2000
 END IF
 6000 i = i + 1
 IF (i - nas > 0) THEN
   GO TO  3000
 ELSE
   GO TO  5100
 END IF
 7000 RETURN
END SUBROUTINE fwmw
