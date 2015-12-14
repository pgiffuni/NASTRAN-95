SUBROUTINE incro (ax,ay,az,ax1,ay1,az1,ax2,ay2,az2,sgr,cgr,sgs,  &
        cgs,kr,fl,beta,sdelx,dely,delr,deli)
     
!     CALCULATES THE UNSTEADY PART OF THE INFLUENCE COEFFICIENT MATRIX
!     ELEMENTS USING  SUBROUTINES  KERNEL, IDF1  AND  IDF2
 
 
 REAL, INTENT(IN)                         :: ax
 REAL, INTENT(IN)                         :: ay
 REAL, INTENT(IN)                         :: az
 REAL, INTENT(IN)                         :: ax1
 REAL, INTENT(IN)                         :: ay1
 REAL, INTENT(IN)                         :: az1
 REAL, INTENT(IN)                         :: ax2
 REAL, INTENT(IN)                         :: ay2
 REAL, INTENT(IN)                         :: az2
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: sgs
 REAL, INTENT(IN)                         :: cgs
 REAL, INTENT(IN OUT)                     :: kr
 REAL, INTENT(IN)                         :: fl
 REAL, INTENT(IN)                         :: beta
 REAL, INTENT(IN)                         :: sdelx
 REAL, INTENT(IN)                         :: dely
 REAL, INTENT(OUT)                        :: delr
 REAL, INTENT(OUT)                        :: deli
 REAL :: k10,k20,k1rt1,k1it1,k2rt2p,k2it2p,k10t1,k20t2p, m
 COMMON /dlm/ k10,k20,k1rt1,k1it1,k2rt2p,k2it2p,k10t1,k20t2p
 COMMON /kds/ ind
 
!     DKRO = REAL PART OF THE PLANAR KERNEL  *  OUTBOARD POINT
!     DKIO = IMAGINARY PART OF THE PLANAR KERNEL  *  OUTBOARD POINT
!     XKRO = REAL PART OF THE NONPLANAR KERNEL  *  OUTBOARD POINT
!     XKIO = IMAGINARY PART OF THE NONPLANAR KERNEL  *  OUTBOARD POINT
!     DKRI = REAL PART OF THE PLANAR KERNEL  *   INBOARD POINT
!     DKII = IMAGINARY PART OF THE PLANAR KERNEL  *   INBOARD POINT
!     XKRI = REAL PART OF THE NONPLANAR KERNEL  *   INBOARD POINT
!     XKII = IMAGINARY PART OF THE NONPLANAR KERNEL  *   INBOARD POINT
 
 ind   = 1
 m     = SQRT(1.0 - beta**2)
 br    = fl/2.
 eps   = 0.00001
 pi    = 3.14159265
 xdelx = sdelx
 xdely = dely
 ee    = 0.5*xdely
 e2    = ee**2
 delr  = 0.0
 deli  = 0.0
 at1s  = 0.0
 at2s  = 0.0
 t1    = 0.0
 t2    = 0.0
 count = 0.
 x0    = ax
 y0    = ay
 z0    = az
 80 CONTINUE
 CALL tker (x0,y0,z0,kr,br,sgr,cgr,sgs,cgs,t1,t2,m)
 at1   = ABS(t1)
 at2   = ABS(t2)
 IF (at1 > at1s) at1s = at1
 IF (at2 > at2s) at2s = at2
 IF (count < 0.0) THEN
   GO TO   130
 ELSE IF (count == 0.0) THEN
   GO TO    90
 ELSE
   GO TO   150
 END IF
 90 dkrc  = k1rt1 - k10t1
 dkic  = k1it1
 xkrc  = k2rt2p - k20t2p
 xkic  = k2it2p
 at2   = ABS(t2)
 count = -1.
 x0    = ax1
 y0    = ay1
 z0    = az1
 GO TO 80
 130 dkri  = k1rt1 - k10t1
 dkii  = k1it1
 xkri  = k2rt2p - k20t2p
 xkii  = k2it2p
 count = 1.
 x0    = ax2
 y0    = ay2
 z0    = az2
 GO TO 80
 150 dkro  = k1rt1 - k10t1
 dkio  = k1it1
 xkro  = k2rt2p - k20t2p
 xkio  = k2it2p
 x0    = ax
 y0    = ay
 z0    = az
 zero  = 0.0
 xiijr = 0.
 xiiji = 0.
 diijr = 0.0
 diiji = 0.0
 xmult = xdelx/(8.0*pi)
 IF (y0 == zero .AND.  z0 == zero) GO TO 220
 IF (z0 == zero .AND. sgs == zero) GO TO 230
 eta01 = y0*cgs + z0*sgs
 zet01 =-y0*sgs + z0*cgs
 azet0 = ABS(zet01)
 IF (azet0 <= 0.0001) zet01 = 0.
 r1sqx = eta01**2 + zet01**2
 210 are   = (dkri - 2.*dkrc + dkro)/(2.0*e2)
 aim   = (dkii - 2.*dkic + dkio)/(2.0*e2)
 bre   = (dkro - dkri)/(2.0*ee)
 bim   = (dkio - dkii)/(2.0*ee)
 cre   =  dkrc
 cim   =  dkic
 GO TO 250
 220 eta01 = 0.0
 zet01 = 0.0
 r1sqx = 0.0
 GO TO  210
 230 eta01 = y0*cgs
 zet01 = 0.
 r1sqx = eta01**2
 GO TO  210
 250 CONTINUE
 IF (at1s == 0.0) GO TO 255
 CALL idf1 (ee,e2,eta01,zet01,are,aim,bre,bim,cre,cim,r1sqx,xiijr, xiiji)
 delr  = xmult*xiijr
 deli  = xmult*xiiji
 255 CONTINUE
 IF (at2s == 0.0) GO TO 260
 a2r   = (xkri - 2.0*xkrc + xkro)/(2.0*e2)
 a2i   = (xkii - 2.0*xkic + xkio)/(2.0*e2)
 b2r   = (xkro - xkri)/(2.0*ee)
 b2i   = (xkio - xkii)/(2.0*ee)
 c2r   =  xkrc
 c2i   =  xkic
 CALL idf2 (ee,e2,eta01,zet01,a2r,a2i,b2r,b2i,c2r,c2i,r1sqx,diijr, diiji)
 delr  = delr + xmult*diijr
 deli  = deli + xmult*diiji
 260 CONTINUE
 
 RETURN
END SUBROUTINE incro
