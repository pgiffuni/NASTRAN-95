SUBROUTINE tvor (sl1,cl1,tl1,sl2,cl2,tl2,sgs,cgs,sgr,cgr,x01,x02,  &
        y0,z0,e,beta,cbar,fmach,kr,bre,bim)
     
!     NORMALWASH AT A POINT (X,Y,Z) - OF A SURFACE DIHEDRAL -
!     DUE TO A TRAPEZOIDAL UNSTEADY VORTEX RING OF UNIT STRENGTH.
 
!     THIS SUBROUTINE CALLS - SNPDF, IDF1, IDF2, FLLD
 
!     SL1, CL1, TL1  SIN(LAMBDA-1), COS(LAMBDA-1), TAN(LAMBDA-1)
!     SL2, CL2, TL2  SIN(LAMBDA-2), .....
!     SGS, CGS       SIN(GAMMA-S),  ....
!     SGR, CGR       SIN(GAMMA-R),  ....
!     X01            X-XI1
!     X02            X-XI2
!     Y0             Y - ETA
!     Z0             Z - ZETA
!     E
!     BETA           SQRT(1-FMACH**2)
!     CV
!     BR
!     FMACH          MACH NO.
!     BRE            REAL PART OF B      (RETURNED)
!     BIM            IMAGINARY PART OF B (RETURNED)
 
 
 REAL, INTENT(IN)                         :: sl1
 REAL, INTENT(IN)                         :: cl1
 REAL, INTENT(IN)                         :: tl1
 REAL, INTENT(IN)                         :: sl2
 REAL, INTENT(IN)                         :: cl2
 REAL, INTENT(IN)                         :: tl2
 REAL, INTENT(IN)                         :: sgs
 REAL, INTENT(IN)                         :: cgs
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: x01
 REAL, INTENT(IN)                         :: x02
 REAL, INTENT(IN)                         :: y0
 REAL, INTENT(IN)                         :: z0
 REAL, INTENT(IN)                         :: e
 REAL, INTENT(IN OUT)                     :: beta
 REAL, INTENT(IN OUT)                     :: cbar
 REAL, INTENT(IN OUT)                     :: fmach
 REAL, INTENT(IN OUT)                     :: kr
 REAL, INTENT(OUT)                        :: bre
 REAL, INTENT(OUT)                        :: bim
 REAL :: kd1, kd2
 
!     VARIABLES DIMENSIONED (2), FIRST WORD IS THE REAL PART OF THE
!     VALUE AND THE SECOND IS THE IMAGINARY PART
 
 DIMENSION      dki(2), dkc(2), dko(2), kd1(2), kd2(2)
 
 DATA  pi48   / 150.79644720 /
 
!     CALCULATE  BS
 
 l  = 1
 cv = x01 - x02
 sl = sl1
 cl = cl1
 tl = tl1
 x0 = x01
 ee = e**2
 te = 2.0*e
 ASSIGN 50 TO isnp
 
!     CALL SNPDF
 
 GO TO 1000
 50 bs = dij
 sl = sl2
 cl = cl2
 tl = tl2
 x0 = x02
 ASSIGN 100 TO isnp
 
!     CALL SNPDF
 
 GO TO 1000
 100 bs = bs - dij
 
!     CALCULATE   DELTA-B
!     LIMITS FOR SMALL VALUES OF RADII
 
 eps = 0.25*ee
 ib  = 0
 fb  = 1.0
 fc  = 4.0
 
!     FIRST CALC.
!     DELTA-KD- 1I, 1C, AND 1O
 
 etl1 = e*tl1
 etl2 = e*tl2
 esgs = e*sgs
 ecgs = e*cgs
 
 dx01 = x01 + etl1
 dx02 = x02 + etl2
 dy0  = y0  + ecgs
 dz0  = z0  + esgs
 ASSIGN 200 TO iflld
 
!     CALCULATE  R-I  SQUARED AND CALL FLLD IF LARGE ENOUGH
 
 r2  = dy0**2 + dz0**2
 IF (r2 >= eps) GO TO 2000
 ib  = 1
 fc  = 6.0
 fb  = 0.0
 GO TO 230
 200 dki(1) = kd1(1)/r2 + kd2(1)/r4
 dki(2) = kd1(2)/r2 + kd2(2)/r4
 
!     KD1C AND KD2C
 
 230 dx01 = x01
 dx02 = x02
 dy0  = y0
 dz0  = z0
 ASSIGN 300 TO iflld
 
!     CALCULATE  R-C  SQUARED AND CALL FLLD IF LARGE ENOUGH
 
 r2 = dy0**2 + dz0**2
 IF (r2 >= eps) GO TO 2000
 fc = 0.0
 fb = 3.0
 GO TO 330
 300 dkc(1) = kd1(1)/r2 + kd2(1)/r4
 dkc(2) = kd1(2)/r2 + kd2(2)/r4
 
!     KD1O AND KD2O
!     SKIP IF  R-I IS TOO SMALL
 
 330 IF (ib /= 0) GO TO 430
 dx01 = x01 - etl1
 dx02 = x02 - etl2
 dy0  = y0  - ecgs
 dz0  = z0  - esgs
 ASSIGN 400 TO iflld
 
!     CALCULATE  R-O  SQUARED AND CALL FLLD IF LARGE ENOUGH
 
 r2 = dy0**2 + dz0**2
 IF (r2 >= eps) GO TO 2000
 fb = 0.0
 fc = 6.0
 ib = 1
 GO TO 430
 400 dko(1) = kd1(1)/r2 + kd2(1)/r4
 dko(2) = kd1(2)/r2 + kd2(2)/r4
 
 430 coef = 1.0/pi48
 bre  = bs/(te*cv) - coef*(fb*(dki(1) + dko(1)) + fc*dkc(1))
 bim  =            - coef*(fb*(dki(2) + dko(2)) + fc*dkc(2))
 RETURN
 
 1000 CALL snpdf (sl,cl,tl,sgs,cgs,sgr,cgr,x0,y0,z0,e,dij,beta,cv)
 GO TO isnp, (50,100)
 
 2000 CALL flld (dx01,dx02,dy0,dz0,sgr,cgr,sgs,cgs,kr,cbar,fmach,e,l,  &
     kd1(1),kd1(2),kd2(1),kd2(2))
 r4 = r2*r2
 GO TO iflld, (200,300,400)
 
END SUBROUTINE tvor
