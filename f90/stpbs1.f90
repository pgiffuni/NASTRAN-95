SUBROUTINE stpbs1(x,ncode,bj1,by1)
!     SUBROUTINE BES1        J AND Y BESSEL FUNCTIONS OF FIRST ORDER
!     E. ALBANO, ORGN 3721, EXT 1022, OCT 1967
!     COMPUTES J1(X) IF X IS GREATER THAN -3.
!     COMPUTES Y1(X) IF (X IS GREATER THAN E AND NCODE = 1 ),
!           WHERE
 
 REAL, INTENT(IN)                         :: x
 INTEGER, INTENT(IN OUT)                  :: ncode
 REAL, INTENT(OUT)                        :: bj1
 REAL, INTENT(OUT)                        :: by1
 INTEGER :: NAME(2)
 DATA NAME /4HSTPB,4HS1  /
 
 e=0.00001
!                 REF. US DEPT OF COMMERCE HANBOOK (AMS 58)  PG. 370
 a=ABS(x)
 IF(a-3. > 0.0) THEN
   GO TO   100
 END IF
 10 z=x*x/9.
 bj1=x*(0.5+z*(-0.56249985+z*(0.21093573+z*(-0.03954289+z*  &
     (0.00443319+z*(-0.00031761+0.00001109*z))))))
 IF(ncode-1 == 0) THEN
   GO TO    20
 END IF
 15 RETURN
 20 IF(x-e < 0.0) THEN
   GO TO   200
 END IF
 25 by1=0.63661977*bj1*(ALOG(x)-.69314718)+(-0.6366198 +z*  &
     (0.2212091+z*(2.1682709 +z*(-1.3164827+z*(0.3123951+z*  &
     (-0.0400976+0.0027873*z))))))/x
 RETURN
 100 IF(x > 0.0) THEN
   GO TO   110
 ELSE
   GO TO   250
 END IF
 110 u=1./SQRT(x)
 z=3./x
 w=0.79788456+z*(0.00000156+z*(0.01659667+z*(0.00017105+z*  &
     (-0.00249511+z*(0.00113653-0.00020033*z)))))
 t=x-2.35619449+z*(0.12499612+z*(0.00005650+z*(-0.00637879+z*  &
     (0.00074348+z*(0.00079824-0.00029166*z)))))
 uw=u*w
 bj1=uw*COS(t)
 IF(ncode-1 == 0) THEN
   GO TO   120
 ELSE
   GO TO    15
 END IF
 120 by1=uw*SIN(t)
 1000 RETURN
 200 CONTINUE
 250 CONTINUE
 CALL mesage(-7,0,NAME)
 GO TO 1000
END SUBROUTINE stpbs1
