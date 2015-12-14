SUBROUTINE perpec (x,stereo)
     
 
 REAL, INTENT(IN OUT)                     :: x(3,1)
 INTEGER, INTENT(IN OUT)                  :: stereo
 INTEGER :: fvp,prject,gp
 REAL :: MIN,MAX
 DOUBLE PRECISION :: diam,r
 COMMON /BLANK /  skp(5),ngpset
 COMMON /xxparm/  skpplt(6),penpap(30),scalx1,objmod,scalx2(3),  &
     view(15),fvp,r0,s0l,s0r,t0,d0,d02,d03,prject,s0s
 COMMON /rstxxx/  cstm(3,3),MIN(3),MAX(3),d(3),aver(3)
 DATA    rdist /  29. /
 
!                               I====================I
!                         T     I                    I
!                         1     I     PROJECTION     I
!                         1     I                    I
!                         1     I       PLANE        I
!                         1     I                    I
!                         1     I====================I
!                         1    /                    /
!                         1   /                    /
!                         1  /   * OBSERVER       /
!                         1 /    1               /D0
!                         1/     1              /
!                         +--------------------/-----S
!                        /       1  /         /
!                       /      T01 /R0
!                      /         1/         /
!                     /----------+- - - - -/
!                    /    S0
!                   R
 
 IF (prject == 1) GO TO 140
 IF (fvp    == 0) GO TO 110
 IF (prject == 3) GO TO 100
 
!     PERSPECTIVE PROJECTION...FIND VANTAGE POINT
 
 r    = d(1)**2 + d(2)**2 + d(3)**2
 diam = DSQRT(r)
 r0   = 2.*diam + aver(1)
 s0l  = aver(2)
 t0   = diam + aver(3)
 d0   = 1.5*diam
 GO TO 110
 
!     STEREO PROJECTION...FIND VANTAGE POINT
 
 100 r0  = rdist + aver(1)*objmod
 s0l = aver(2)*objmod - s0s/2.
 s0r = aver(2)*objmod + s0s/2.
 t0  = aver(3)*objmod
 d0  = d03
 GO TO 140
 
 110 scal = 1.
 IF (prject == 3) scal = objmod
 slr = s0l
 IF (stereo /= 0) slr = s0r
 DO  gp = 1,ngpset
   r = d0/(r0-scal*x(1,gp))
   s = slr + r*(scal*x(2,gp)-slr)
   t = t0  + r*(scal*x(3,gp)-t0 )
   x(2,gp) = s
   x(3,gp) = t
   IF (prject == 3) CYCLE
   MIN(2) = AMIN1(MIN(2),s)
   MIN(3) = AMIN1(MIN(3),t)
   MAX(2) = AMAX1(MAX(2),s)
   MAX(3) = AMAX1(MAX(3),t)
 END DO
 IF (prject == 3) GO TO 140
 
!     FIND MINIMA + MAXIMA DIFFERENCES + AVERAGES
 
 DO  i = 2,3
   d(i) = MAX(i) - MIN(i)
   aver(i) = (MAX(i)+MIN(i))/2.
 END DO
 
 140 RETURN
END SUBROUTINE perpec
