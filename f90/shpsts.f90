SUBROUTINE shpsts (sigma,vonms,sigp)
     
!     TO CALCULATE PRINCIPAL STRESSES AND THEIR ANGLES FOR THE
!     ISOPARAMETRIC SHELL ELEMENTS
 
 
!     INPUT :
!           SIGMA  - ARRAY OF 3 STRESS COMPONENTS
!           VONMS  - LOGICAL FLAG INDICATING THE PRESENCE OF VON-MISES
!                    STRESS REQUEST
!     OUTPUT:
!           SIGP   - ARRAY OF PRINCIPAL STRESSES
 
 
 
 REAL, INTENT(IN)                         :: sigma(3)
 LOGICAL, INTENT(IN OUT)                  :: vonms
 REAL, INTENT(OUT)                        :: sigp(4)
 
!WKBNB 7/94 SPR94004
 LOGICAL :: ostrai
 COMMON / BLANK / app(2), sort2, idum(2), comps, skp(4), ostrai  &
     ,                sk2(39), midve
!WKBNE 7/94 SPR94004
 REAL :: sig,proj,taumax,eps,txy2
 DATA    eps / 1.0E-11 /
 
 
!     CALCULATE PRINCIPAL STRESSES
 
 sig  = 0.5*(sigma(1)+sigma(2))
 proj = 0.5*(sigma(1)-sigma(2))
 taumax = proj*proj + sigma(3)*sigma(3)
!WKBI 7/94 SPR94004
 IF ( ostrai ) taumax = proj*proj + sigma(3)*sigma(3)/4.
 IF (taumax /= 0.0) taumax = SQRT(taumax)
 IF (taumax <= eps) taumax = 0.0
 
!     CALCULATE THE PRINCIPAL ANGLE
 
 txy2 = sigma(3)*2.0
 proj = proj*2.0
 sigp(1) = 0.0
 IF (ABS(txy2) > eps .OR. ABS(proj) > eps)  &
     sigp(1) = 28.64788976*ATAN2(txy2,proj)
!                   28.64788976 = 90./PI
 
 sigp(2) = sig + taumax
 sigp(3) = sig - taumax
 sigp(4) = taumax
 
!     OUTPUT VON MISES YIELD STRESS IF REQUESTED
 
 IF (.NOT.vonms) RETURN
 sig = sigp(2)*sigp(2) + sigp(3)*sigp(3) - sigp(2)*sigp(3)
 IF (sig /= 0.0) sig = SQRT(sig)
 IF (sig <= eps) sig = 0.0
 sigp(4) = sig
 
 RETURN
END SUBROUTINE shpsts
