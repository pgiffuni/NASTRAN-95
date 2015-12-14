SUBROUTINE flld (x01,x02,y0,z0,sgr,cgr,sgs,cgs,kr,cbar,fmach,e,  &
        l,kd1r,kd1i,kd2r,kd2i)
     
!     CALCULATION OF THE NUMERATOR OF A DOUBLET LINE OF FINITE LENGTH.
!     LIKE KERN, THERE ARE TWO OUTPUT COMPLEX VALUES REPRESENTED BY
!     FOUR REAL NUMBERS AND AN INPUT OPTION.
 
!     WRITTEN BY D. H. LARSON, STRUCTURAL MECHANICS MDAC 11/70
 
!     X01  -   X - XI1
!     X02  -   X - XI2
!     Y0   -   Y - ETA
!     Z0   -   Z - ZETA
!     SGR  -   SIN ( GAMMA-R)
!     CGR  -   COS ( GAMMA-R)
!     SGS  -   SIN ( GAMMA-S)
!     CGS  -   COS ( GAMMA-S)
!     KR   -   REDUCED FREQUENCY
!     BR   -   REFERENCE LENGTH
!     FMACH-   MACH NUMBER
!     E    -
!     L    -   OPTION FLAG USED IN TKER
!     KD1R -   REAL PART OF  KD1
!     KD1I -   IMAGINARY PART OF KD1
!     KD2R -   REAL PART OF  KD2
!     KD2I -   IMAGINARY PART OF KD2
 
 
 REAL, INTENT(IN)                         :: x01
 REAL, INTENT(IN)                         :: x02
 REAL, INTENT(IN OUT)                     :: y0
 REAL, INTENT(IN OUT)                     :: z0
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN OUT)                     :: sgs
 REAL, INTENT(IN OUT)                     :: cgs
 REAL, INTENT(IN)                         :: kr
 REAL, INTENT(IN)                         :: cbar
 REAL, INTENT(IN OUT)                     :: fmach
 REAL, INTENT(IN OUT)                     :: e
 INTEGER, INTENT(IN OUT)                  :: l
 REAL, INTENT(OUT)                        :: kd1r
 REAL, INTENT(OUT)                        :: kd1i
 REAL, INTENT(OUT)                        :: kd2r
 REAL, INTENT(OUT)                        :: kd2i
 REAL :: kk1r,kk1i,kk2r,kk2i, k10t1, k20t2p,k1rt1,k10,k2it2p,k20,k2rt2p,k1it1
 COMPLEX :: kd1,kd2,k1xi1,k1xi2,temp1,temp2,k2xi1,k2xi2
 COMMON  /kds/  ind,kk1r,kk1i,kk2r,kk2i
 COMMON  /dlm/  k10,k20,k1rt1,k1it1,k2rt2p,k2it2p,k10t1,k20t2p
 
!     X01 = X-XI1  AND  X02 = X-XI2, DELXI = XI2-XI1
 
 delxi = x01 - x02
 
!     FULL KERNEL FROM -TKER-
 
 ind  = 0
 kd1r = 0.0
 kd2r = 0.0
 t1   = kr*delxi/cbar
 br   = cbar/2.0
 st1  = SIN(t1)
 ct1  = COS(t1)
 i    = 1
 x0   = x01
 
 10 CALL tker (x0,y0,z0,kr,br,sgr,cgr,sgs,cgs,rt1,rt2,fmach)
 
 SELECT CASE ( i )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 40
 END SELECT
 30 k1xi1 = CMPLX(kk1r,kk1i)
 k2xi1 = CMPLX(kk2r,kk2i)
 IF (l == 0) GO TO 35
 kd1r  = kd1r - k10t1
 kd2r  = kd2r - k20t2p
 35 CONTINUE
 
!     NOW GO CALCULATE FOR XI = XI2
 
 x0 = x02
 i  = 2
 GO TO 10
 
 40 k1xi2 = CMPLX(kk1r,kk1i)
 k2xi2 = CMPLX(kk2r,kk2i)
 IF (l == 0) GO TO 50
 kd1r  = kd1r + k10t1
 kd2r  = kd2r + k20t2p
 50 CONTINUE
 
 temp1 = CMPLX(ct1, st1)
 temp2 = CMPLX(ct1,-st1)
 
!     DESIRED RESULTS (COMPLEX)
 
 kd1   = k1xi1*temp1 - k1xi2*temp2
 kd2   = k2xi1*temp1 - k2xi2*temp2
 
!     CONVERT TO REAL AND IMAGINARY PARTS
 
 kd1r  = REAL (kd1) + kd1r
 kd1i  = AIMAG(kd1)
 kd2r  = REAL (kd2) + kd2r
 kd2i  = AIMAG(kd2)
 RETURN
END SUBROUTINE flld
