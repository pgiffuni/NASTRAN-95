SUBROUTINE t3gems (ierr,egpdt,iorder,gb,gs,lx,ly,edglen,shrflx,  &
        aic,jog,jok,k11,k22)
     
!     SINGLE PRECISION ROUTINE TO SET UP THE REQUIRED SHEAR-RELATED
!     TRANSFORMATION TO RELIEVE THE TRIA3 GEOMETRY BIAS IN BENDING.
 
!     INPUT :
!           EGPDT  - BGPDT DATA IN ELEMENT COORD. SYSTEM
!           IORDER - ARRAY OF ORDER INDICATORS FOR REARRANGED DATA
!           GB     - ARRAY OF BENDING MATERIAL PROPERTIES
!           GS     - ARRAY OF SHEAR   MATERIAL PROPERTIES
!           LX     - DIMENSION OF ELEMENT ALONG X-AXIS
!           LY     - DIMENSION OF ELEMENT ALONG Y-AXIS
!           EDGLEN - EDGE LENGTHS
!           SHRFLX - LOGICAL INDICATING THE PRESENCE OF SHEAR FLEX
!     OUTPUT:
!           IERR   - ERROR FLAG
!           AIC    - TRANSFORMATION TO RELIEVE GEOMETRY BIAS
!           JOG    - SHEAR   STIFFNESS FACTOR
!           JOK    - BENDING STIFFNESS FACTOR
!           K11    - BENDING STIFFNESS FACTOR
!           K22    - BENDING STIFFNESS FACTOR
 
 
!     [C]    - TRANSFORMATION TO YIELD GAMMAT ALONG THE ELEMENT SIDES.
 
!     [AA]   - TRANSFORMATION FROM GAMMA0 (AT THE ELEMENT CENTER) TO
!              GAMMAT (ALONG THE ELEMENT SIDES).
 
!                  -1
!     [AIC]  - [AA]  [C]
 
 
 
 INTEGER, INTENT(OUT)                     :: ierr
 REAL, INTENT(IN)                         :: egpdt(4,3)
 INTEGER, INTENT(IN)                      :: iorder(3)
 REAL, INTENT(IN)                         :: gb(9)
 REAL, INTENT(IN)                         :: gs(4)
 REAL, INTENT(IN)                         :: lx
 REAL, INTENT(IN)                         :: ly
 REAL, INTENT(IN)                         :: edglen(3)
 LOGICAL, INTENT(IN OUT)                  :: shrflx
 REAL, INTENT(OUT)                        :: aic(18)
 REAL, INTENT(OUT)                        :: jog
 REAL, INTENT(OUT)                        :: jok
 REAL, INTENT(OUT)                        :: k11
 REAL, INTENT(OUT)                        :: k22
 
 INTEGER :: INDEX(3,3)
 REAL :: xx(3),yy(3),aa(9),h1,h2,bdum(3), determ,cosa,sina,cosb,sinb,cosc,sinc
 
 
 ierr = 0
 DO  i = 1,3
   DO  j = 1,3
     jo = iorder(j)
     IF (i /= jo) CYCLE
     xx(i) = egpdt(2,j)
     yy(i) = egpdt(3,j)
   END DO
 END DO
 
 cosa = ((xx(2)-xx(1))/edglen(1))
 sina = ((yy(2)-yy(1))/edglen(1))
 cosb = ((xx(3)-xx(2))/edglen(2))
 sinb = ((yy(3)-yy(2))/edglen(2))
 cosc = ((xx(1)-xx(3))/edglen(3))
 sinc = ((yy(1)-yy(3))/edglen(3))
 
 aa(1) = sina
 aa(2) = cosa
 aa(3) = 1.0
 aa(4) = sinb
 aa(5) = cosb
 aa(6) = 1.0
 aa(7) = sinc
 aa(8) = cosc
 aa(9) = 1.0
 
 CALL invers (3,aa,3,bdum,0,determ,ising,INDEX)
 IF (ising /= 1) GO TO 30
 
 aic( 1) = aa(1)*sina
 aic( 2) = aa(1)*cosa
 aic( 3) = aa(2)*sinb
 aic( 4) = aa(2)*cosb
 aic( 5) = aa(3)*sinc
 aic( 6) = aa(3)*cosc
 aic( 7) = aa(4)*sina
 aic( 8) = aa(4)*cosa
 aic( 9) = aa(5)*sinb
 aic(10) = aa(5)*cosb
 aic(11) = aa(6)*sinc
 aic(12) = aa(6)*cosc
 aic(13) = aa(7)*sina
 aic(14) = aa(7)*cosa
 aic(15) = aa(8)*sinb
 aic(16) = aa(8)*cosb
 aic(17) = aa(9)*sinc
 aic(18) = aa(9)*cosc
 
!     CALCULATE THE BENDING STIFFNESS FACTORS
 
 h1  = ly
 h2  = lx
 k11 = 1.0/(h1*h1)*gb(5)
 k22 = 1.0/(h2*h2)*gb(1)
 
 jok = k11*k22
 IF (jok /= 0.0) jok = 1.0/jok
 jog = 0.0
 IF (shrflx) jog = gs(1)*gs(4) - gs(2)*gs(3)
 IF (jog /= 0.0) jog = 1.0/jog
 GO TO 40
 
 30 ierr = 1
 40 RETURN
END SUBROUTINE t3gems
