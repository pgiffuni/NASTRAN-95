SUBROUTINE failur (fthr,ultstn,stresl,findex)
     
!      THIS ROUTINE COMPUTES THE FAILURE INDEX OF A LAYER
!      IN A LAMINATED COMPOSITE ELEMENT USING ONE OF THE
!      FOLLOWING FIVE FAILURE THEORIES CURRENTLY AVAILABLE
 
!        1.   HILL
!        2.   HOFFMAN
!        3.   TSAI-WU
!        4.   MAX STRESS
!        5.   MAX STRAIN
 
!        DEFINITIONS
!        -----------
!        XT = ULTIMATE UNIAXIAL TENSILE STRENGTH IN THE FIBER
!             DIRECTION
!        XC = ULTIMATE UNIAXIAL COMPRESSIVE STRENGTH IN THE
!             FIBER DIRECTION
!        YT = ULTIMATE UNIAXIAL TENSILE STRENGTH PERPENDICULAR TO
!             THE FIBER DIRECTION
!        YC = ULTIMATE UNIAXIAL COMPRESSIVE STRENGTH PERPENDICULAR
!             TO THE FIBER DIRECTION
!        S  = ULTIMATE PLANAR SHEAR STRENGTH UNDER PURE SHEAR
!             LOADING
!        SIMILARILY FOR THE ULTIMATE STRAINS
 
 
 INTEGER, INTENT(IN OUT)                  :: fthr
 REAL, INTENT(IN)                         :: ultstn(6)
 REAL, INTENT(IN)                         :: stresl(3)
 REAL, INTENT(OUT)                        :: findex
 
 
 COMMON /system/ ibuf,nout
 
!**** CHECK FOR ZERO STRENGTH VALUES
 
 DO  i = 1,5
   IF (ultstn(i) == 0.0) GO TO 90
 END DO
 
!**** ULTIMATE STRENGTH VALUES
 
 xt      = ultstn(1)
 xc      = ultstn(2)
 yt      = ultstn(3)
 yc      = ultstn(4)
 s       = ultstn(5)
 f12     = ultstn(6)
 
!**** LAYER STRESSES
 
 sig1    = stresl(1)
 sig2    = stresl(2)
 tau12   = stresl(3)
 
!**** LAYER STRAINS
 
 eps1    = stresl(1)
 eps2    = stresl(2)
 gama    = stresl(3)
 
 
 SELECT CASE ( fthr )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 60
   CASE (    5)
     GO TO 70
 END SELECT
 
!     H I L L   F A I L U R E  T H E O R Y
!     ====================================
 
 20 x = xt
 IF (sig1 < 0.0) x = xc
 
 y = yt
 IF (sig2 < 0.0) y = yc
 
 xx = xt
 IF ((sig1*sig2) < 0.0) xx = xc
 
 findex = ( sig1*sig1 )/( x*x )
 findex = findex + ( sig2*sig2 )/( y*y )
 findex = findex - ( sig1*sig2 )/( xx*xx )
 findex = findex + ( tau12*tau12 )/( s*s )
 GO TO 80
 
 
!     H O F F M A N  F A I L U R E  T H E O R Y
!     =========================================
 
 30 findex = ( 1.0/xt - 1.0/xc )*sig1
 findex = findex + ( 1.0/yt - 1.0/yc )*sig2
 findex = findex + ( sig1*sig1 )/( xt*xc )
 findex = findex + ( sig2*sig2 )/( yt*yc )
 findex = findex + ( tau12*tau12 )/( s*s )
 findex = findex - ( sig1*sig2 )/( xt*xc )
 GO TO 80
 
 
!     T S A I-W U  F A I L U R E  T H E O R Y
!     =======================================
 
!**** CHECK STABILITY CRITERIA FOR THE INTERACTION TERM F12
 40 IF (f12 == 0.0) GO TO 50
 
 crit = ( 1.0/(xt*xc) )*( 1.0/(yt*yc) ) - ( f12*f12 )
 IF (crit > 0.0) GO TO 50
 
!**** IF STABILITY CRITERIA IS VIOLATED THEN SET THE
!     F12 THE INTERACTION TERM TO ZERO
 
 f12 = 0.0
 
 
 50 findex = ( 1.0/xt - 1.0/xc )*sig1
 findex = findex + ( 1.0/yt - 1.0/yc )*sig2
 findex = findex + ( sig1*sig1 )/( xt*xc )
 findex = findex + ( sig2*sig2 )/( yt*yc )
 findex = findex + ( tau12*tau12 )/( s*s )
 IF (f12 == 0.0) GO TO 80
 findex = findex + ( 2.0*f12*sig1*sig2 )
 GO TO 80
 
 
!     M A X  S T R E S S  F A I L U R E  T H E O R Y
!     ==============================================
 
 60 fi1 = sig1/xt
 IF (sig1 < 0.0) fi1 = sig1/xc
 
 fi2 = sig2/yt
 IF (sig2 < 0.0) fi2 = sig2/yc
 
 fi12 = ABS(tau12)/s
 
 findex = fi1
 IF (fi2  > findex) findex = fi2
 IF (fi12 > findex) findex = fi12
 GO TO 80
 
 
!     M A X  S T R A I N  F A I L U R E  T H E O R Y
!     ==============================================
 
 70 fi1 = eps1/xt
 IF (eps1 < 0.0) fi1 = eps1/xc
 
 fi2 = eps2/yt
 IF (eps2 < 0.0) fi2 = eps2/yc
 
 fi12 = ABS(gama)/s
 
 findex = fi1
 IF (fi2  > findex) findex = fi2
 IF (fi12 > findex) findex = fi12
 
 80 CONTINUE
 
 RETURN
 
 
!     NON-FATAL ERROR
 
 
 90 findex = 0.0
 RETURN
END SUBROUTINE failur
