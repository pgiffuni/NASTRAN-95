SUBROUTINE dzy (x,y,z,sgr,cgr,xi1,xi2,eta,zeta,ar,ao,kr,cbar,  &
        beta,fmach,lsh,idzdy,dzdyr,dzdyi)
     
!     CALCULATION OF THE DZ AND DY MATRICES USED IN SLENDER BODY FLOW
 
!     X        X- COORDINATE OF THE RECEIVING POINT
!     Y        Y - COORDINATE OF THE RECEIVING POINT
!     Z        Z - COORDINATE OF THE RECEIVING POINT
!     SGR      SINE OF THE RECEIVING POINT DIHEDRAL ANGLE
!     CGR      COSINE OF RECEIVING POINT DIHEDRAL ANGLE
!     XI1
!     XI2
!     ETA
!     ZETA
!     AR       ASPECT RATIO OF THE SENDING BODY
!     A0       RADIUS OF THE SENDING BODY
!     KR       REDUCED FREQUENCY
!     CBAR     REFERENCE CHORD LENGTH
!     BETA     SQRT(1.0-M**2)
!     FMACH    MACH NUMBER
!     IDZDY    FLAG INDICATING WHETHER DZ OF DY IS TO BE
!              CALCULATED.  =0  DZ, OTHERWISE DY
!     DZDYR    REAL PART OF DZ OR DY
!     DZDYI    IMAGINARY PART OF DZ OR DY
 
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 REAL, INTENT(IN)                         :: z
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: xi1
 REAL, INTENT(IN)                         :: xi2
 REAL, INTENT(IN)                         :: eta
 REAL, INTENT(IN)                         :: zeta
 REAL, INTENT(IN)                         :: ar
 REAL, INTENT(IN)                         :: ao
 INTEGER, INTENT(IN OUT)                  :: kr
 REAL, INTENT(IN OUT)                     :: cbar
 REAL, INTENT(IN OUT)                     :: beta
 REAL, INTENT(IN OUT)                     :: fmach
 INTEGER, INTENT(IN OUT)                  :: lsh
 INTEGER, INTENT(IN OUT)                  :: idzdy
 REAL, INTENT(OUT)                        :: dzdyr
 REAL, INTENT(OUT)                        :: dzdyi
 REAL :: kd1pr, kd1pi, kd2pr, kd2pi, kd1mr, kd1mi, kd2mr, kd2mi
 DATA   pi16 / 50.265482 /
 
 
!     THE COMPLEX NUMBERS IN THIS ROUTINE ARE TREATED SEPERATLY AS
!     THE REAL PART,  NAME APPENDED BY AN  -R- ,  AND THE
!     IMAGINARY PART, NAME APPENDED BY AN  -I- .
 
 e   = ao*SQRT(ABS(1.0 - ar**2))/2.0
 x01 = x - xi1
 x02 = x - xi2
 
!     CHECK ON INPUT FLAG,  = 0  DZ ,  = 1  DY
 
 IF (idzdy == 1) GO TO 200
 
!     **     **
!     *  D Z  *
!     **     **
 
 sgs = 0.0
 cgs = 1.0
 IF (ar < 1.0)  GO TO 400
 
 z01 = z - (zeta+e)
 z02 = z - (zeta-e)
 y01 = y - eta
 y02 = y01
 GO TO 300
 
!     **     **
!     *  D Y  *
!     **     **
 
 200 sgs = -1.0
 IF (lsh == 1) sgs = 1.0
 cgs =  0.0
 IF (ar > 1.0) GO TO 400
 
 z01 = z - zeta
 z02 = z01
 y01 = y - (eta+e)
 y02 = y - (eta-e)
 
!     ****  DZ AR .GE. 1  ****
!     ****  DY AR .LE. 1  ****
 
 300 CONTINUE
 l   = 0
 z0  = z - zeta
 y0  = y - eta
 
 r1sqr = y01**2 + z01**2
 r2sqr = y02**2 + z02**2
 r1for = r1sqr**2
 r2for = r2sqr**2
 
 CALL flld (x01,x02,y01,z01,sgr,cgr,sgs,cgs,kr,cbar,fmach,e,l,  &
     kd1pr,kd1pi,kd2pr,kd2pi)
 
 IF (ar /= 1.0)  GO TO 320
 
!     IDENTICAL RESULTS FROM FLLD, THEREFORE SKIP SECOND CALL
 
 kd1mr = kd1pr
 kd1mi = kd1pi
 kd2mr = kd2pr
 kd2mi = kd2pi
 GO TO 360
 320 CONTINUE
 CALL flld (x01,x02,y02,z02,sgr,cgr,sgs,cgs,kr,cbar,fmach,e,l,  &
     kd1mr,kd1mi,kd2mr,kd2mi)
 360 CONTINUE
 dzdyr = 0.0
 dzdyi = 0.0
 IF (r1sqr <= 0.0001. OR .r2sqr <= 0.0001) GO TO 370
 
!     REAL
 
 dzdyr = ((kd1pr/r1sqr+kd1mr/r2sqr) + (kd2pr/r1for+kd2mr/r2for)) /pi16*(-1.0)
 
!     IMAGINARY
 
 dzdyi = ((kd1pi/r1sqr+kd1mi/r2sqr) + (kd2pi/r1for+kd2mi/r2for)) /pi16*(-1.0)
 370 CONTINUE
 
 RETURN
 
!     ****   DZ-AR .LT. 1   ****
!     ****   DY-AR .GT. 1   ****
 
 400 sl1 = 0.0
 tl1 = 0.0
 sl2 = 0.0
 tl2 = 0.0
 cl1 = 1.0
 cl2 = 1.0
 e   = 1.732051*e
 y0  = y - eta
 z0  = z - zeta
 
 CALL tvor (sl1,cl1,tl1,sl2,cl2,tl2,sgs,cgs,sgr,cgr,x01,x02,  &
     y0,z0,e,beta,cbar,fmach,kr,dzdyr,dzdyi)
 RETURN
END SUBROUTINE dzy
