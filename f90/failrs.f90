SUBROUTINE failrs (fthr,ultstn,stresl,findex)
     
!     THIS ROUTINE COMPUTES THE FAILURE INDEX OF A LAYER IN A LAMINATED
!     COMPOSITE ELEMENT USING ONE OF THE FOLLOWING FIVE FAILURE THEORIES
!     CURRENTLY AVAILABLE
!        1.   HILL
!        2.   HOFFMAN
!        3.   TSAI-WU
!        4.   MAX STRESS
!        5.   MAX STRAIN
 
 
!     DEFINITIONS
 
!     XT = ULTIMATE UNIAXIAL TENSILE STRENGTH IN THE FIBER DIRECTION
!     XC = ULTIMATE UNIAXIAL COMPRESSIVE STRENGTH IN THE FIBER DIRECTION
!     YT = ULTIMATE UNIAXIAL TENSILE STRENGTH PERPENDICULAR TO THE FIBER
!          DIRECTION
!     YC = ULTIMATE UNIAXIAL COMPRESSIVE STRENGTH PERPENDICULAR TO THE
!          FIBER DIRECTION
!     S  = ULTIMATE PLANAR SHEAR STRENGTH UNDER PURE SHEAR LOADING
 
!     SIMILARILY FOR THE ULTIMATE STRAINS
 
 
 
 INTEGER, INTENT(IN OUT)                  :: fthr
 REAL, INTENT(IN)                         :: ultstn(6)
 REAL, INTENT(IN)                         :: stresl(3)
 REAL, INTENT(OUT)                        :: findex
 
 
 
 
!     CHECK FOR ZERO STRENGTH VALUES
 
 DO  i = 1,5
   IF (ultstn(i) == 0.0) GO TO 700
 END DO
 
!     ULTIMATE STRENGTH VALUES
 
 xt   = ultstn(1)
 xc   = ultstn(2)
 yt   = ultstn(3)
 yc   = ultstn(4)
 s    = ultstn(5)
 f12  = ultstn(6)
 
!     LAYER STRESSES
 
 sig1 = stresl(1)
 sig2 = stresl(2)
 tau12= stresl(3)
 
!     LAYER STRAINS
 
 eps1 = stresl(1)
 eps2 = stresl(2)
 gama = stresl(3)
 
 SELECT CASE ( fthr )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 300
   CASE (    4)
     GO TO 400
   CASE (    5)
     GO TO 500
 END SELECT
 
!     HILL FAILURE THEORY
!     -------------------
 
 100 x = xt
 IF (sig1 < 0.0) x = xc
 
 y = yt
 IF (sig2 < 0.0) y = yc
 
 xx = xt
 IF (sig1*sig2 < 0.0) xx = xc
 
 findex = (sig1*sig1)/(x*x)
 findex = findex + (sig2 * sig2)/(y * y)
 findex = findex - (sig1 * sig2)/(xx*xx)
 findex = findex + (tau12*tau12)/(s * s)
 GO TO 600
 
 
!     HOFFMAN FAILURE THEORY
!     ----------------------
 
 200 findex = (1.0/xt-1.0/xc)*sig1
 findex = findex + (1.0/yt-1.0/yc)*sig2
 findex = findex + (sig1 * sig1)/(xt*xc)
 findex = findex + (sig2 * sig2)/(yt*yc)
 findex = findex + (tau12*tau12)/(s * s)
 findex = findex + (sig1 * sig2)/(xt*xc)
 GO TO 600
 
 
!     TSAI-WU FAILURE THEORY
!     ----------------------
 
!     CHECK STABILITY CRITERIA FOR THE INTERACTION TERM F12
 
 300 IF (f12 == 0.0) GO TO 350
 
 crit = (1.0/(xt*xc))*(1.0/(yt*yc)) - f12*f12
 IF (crit > 0.0) GO TO 350
 
!     IF STABILITY CRITERIA IS VIOLATED THEN SET THE F12 THE INTERACTION
!     TERM TO ZERO
 
 f12 = 0.0
 
 
 350 findex = (1.0/xt-1.0/xc)*sig1
 findex = findex + (1.0/yt-1.0/yc)*sig2
 findex = findex + (sig1 * sig1)/(xt*xc)
 findex = findex + (sig2 * sig2)/(yt*yc)
 findex = findex + (tau12*tau12)/(s * s)
 IF (f12 == 0.0) GO TO 600
 findex = findex + 2.0*f12*sig1*sig2
 GO TO 600
 
 
!     MAX STRESS FAILURE THEORY
!     -------------------------
 
 400 fi1 = sig1/xt
 IF (sig1 < 0.0) fi1 = ABS(sig1/xc)
 
 fi2 = sig2/yt
 IF (sig2 < 0.0) fi2 = ABS(sig2/yc)
 
 fi12 = ABS(tau12)/s
 
 findex = fi1
 IF (fi2  > findex) findex = fi2
 IF (fi12 > findex) findex = fi12
 GO TO 600
 
 
!     MAX STRAIN FAILURE THEORY
!     -------------------------
 
 500 fi1 = eps1/xt
 IF (eps1 < 0.0) fi1 = ABS(eps1/xc)
 
 fi2 = eps2/yt
 IF (eps2 < 0.0) fi2 = ABS(eps2/yc)
 
 fi12 = ABS(gama)/s
 
 findex = fi1
 IF (fi2  > findex) findex = fi2
 IF (fi12 > findex) findex = fi12
 
 600 CONTINUE
 
 RETURN
 
 
!     NON-FATAL ERROR
 
 700 findex = 0.0
 RETURN
END SUBROUTINE failrs
