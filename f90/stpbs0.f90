SUBROUTINE stpbs0(x,ncode,bj0,by0)
!     SUBROUTINE  BES0.  J AND Y BESSEL FUNCTIONS OF ORDER ZERO
!     E. ALBANO, ORGN 3721, EXT 1022, OCT. 1967
!     COMPUTES J0(X)  IF X IS GREATER THAN -3.
!     COMPUTES Y0(X)  IF (X IS GREATER THAN E AND NCODE = 1 ),
!           WHERE
 
 REAL, INTENT(IN)                         :: x
 INTEGER, INTENT(IN OUT)                  :: ncode
 REAL, INTENT(OUT)                        :: bj0
 REAL, INTENT(OUT)                        :: by0
 INTEGER :: NAME(2)
 DATA NAME /4HSTPB,4HS0  /
 
 e=0.00001
!                 REF. US DEPT OF COMMERCE HANDBOOK (AMS 55)  PG. 369
 a=ABS(x)
 IF(a-3. > 0.0) THEN
   GO TO   100
 END IF
 10 z=x*x/9.
 bj0=1.+z*(-2.2499997+z*(1.2656208+z*(-0.3163866+z*(0.0444479  &
     +z*(-0.0039444+z* 0.00021)))))
 IF(ncode-1 == 0) THEN
   GO TO    20
 END IF
 15 RETURN
 20 IF(x-e < 0.0) THEN
   GO TO   200
 END IF
 25 by0=0.63661977*bj0*(ALOG(x)-.69314718)+.36746691+z*(0.60559366+z*  &
     (-0.74350384+z*(0.25300117+z*(-0.04261214+z*(0.00427916 -0.00024846*z)))))
 RETURN
 100 IF(x    > 0.0) THEN
   GO TO   110
 ELSE
   GO TO   250
 END IF
 110 u=1./SQRT(x)
 z=3./x
 w=0.79788456+z*(-0.00000077+z*(-0.0055274+z*(-0.00009512+z*  &
     (0.00137237+z*(-0.00072805+0.00014476*z)))))
 t=x-0.78539816+z*(-0.04166397+z*(-0.00003954+z*(0.00262573+z*  &
     (-0.00054125+z*(-0.00029333+0.00013558*z)))))
 uw=u*w
 bj0=uw*COS(t)
 IF(ncode-1 == 0) THEN
   GO TO   120
 ELSE
   GO TO    15
 END IF
 120 by0=uw*SIN(t)
 1000 RETURN
 200 CONTINUE
 250 CONTINUE
 CALL mesage(-7,0,NAME)
 GO TO 1000
END SUBROUTINE stpbs0
