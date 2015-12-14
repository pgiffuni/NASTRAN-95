SUBROUTINE proces (x)
     
 
 REAL, INTENT(IN OUT)                     :: x(3,1)
 INTEGER :: axes,axis,gp,prject
 REAL :: xmin(3),xmax(3),MIN,MAX
 DOUBLE PRECISION :: sum,v(3)
 COMMON /condas/  consts(5)
 COMMON /BLANK /  skpcom(5),ngpset
 COMMON /xxparm/  pbufsz,ploter(5),penpap(30),scale(5),axes(6),  &
     alpha,beta,gamma,beta13,beta2,view(4), vanput(8),prject
 COMMON /rstxxx/  cstm(3,3),MIN(3),MAX(3),d(3),aver(3), axis(3),SIGN(3)
 COMMON /drwaxs/  g(3,3)
 EQUIVALENCE      (consts(3),rad)
 
!     INITIALIZATION.
 
 DO  i = 1,3
   axis(i) = IABS(axes(i))
   SIGN(i) = 1.
   IF (axes(i) < 0) SIGN(i) = -1.
   MIN(i)  = +1.e+20
   MAX(i)  = -1.e+20
   IF (prject  /= 3) CYCLE
   xmin(i) = +1.e+20
   xmax(i) = -1.e+20
 END DO
 
!     CALCULATE THE CO-ORDINATE SYSTEM ROTATION MATRIX.
 
 IF (beta > -1.e+10) GO TO 20
 IF (prject /= 2) beta = beta13
 IF (prject == 2) beta = beta2
 20 sina = SIN(alpha/rad)
 sinb = SIN(beta /rad)
 sing = SIN(gamma/rad)
 cosa = COS(alpha/rad)
 cosb = COS(beta /rad)
 cosg = COS(gamma/rad)
 
 cstm(1,1) = cosb*cosg
 cstm(2,1) = cosa*sing + sina*sinb*cosg
 cstm(3,1) = sina*sing - cosa*sinb*cosg
 cstm(1,2) =-cosb*sing
 cstm(2,2) = cosa*cosg - sina*sinb*sing
 cstm(3,2) = sina*cosg + cosa*sinb*sing
 cstm(1,3) = sinb
 cstm(2,3) =-sina*cosb
 cstm(3,3) = cosa*cosb
 
!     SWITCH AXES + ROTATE THE GRID POINT CO-ORDINATES.
 
 DO  gp = 1,ngpset
   DO   i = 1,3
     j    = axis(i)
     v(i) = SIGN(i)*x(j,gp)
     IF (prject /= 3) CYCLE
     val  = v(i)
     xmin(i) = AMIN1(xmin(i),val)
     xmax(i) = AMAX1(xmax(i),val)
   END DO
   DO  j = 1,3
     sum = 0.d0
     DO  i = 1,3
       sum = sum + cstm(j,i)*v(i)
     END DO
     val = sum
     x(j,gp) = val
     MIN(j)  = AMIN1(MIN(j),val)
     MAX(j)  = AMAX1(MAX(j),val)
   END DO
 END DO
 
!     CALCULATE THE MINIMA-MAXIMA DIFFERENCES + AVERAGES.
 
 DO  i = 1,3
   IF (prject /= 3) d(i) =  MAX(i) -  MIN(i)
   IF (prject == 3) d(i) = xmax(i) - xmin(i)
   aver(i) = (MAX(i)+MIN(i))/2.
 END DO
 
!     CREATE A X-Y-Z UNIT COORDINATES IN /DRWAXS/ FOR VIEW PLOTTING
 
 DO  i = 1,9
   g(i,1) = 0.0
 END DO
 g(1,1) = 1.0
 g(2,2) = 1.0
 g(3,3) = 1.0
 
 DO  gp = 1,3
   DO   i = 1,3
     j    = axis(i)
     v(i) = SIGN(i)*g(j,gp)
   END DO
   DO  j = 1,3
     sum = 0.d0
     DO  i = 1,3
       sum = sum + cstm(j,i)*v(i)
     END DO
     g(j,gp) = sum
   END DO
 END DO
 
 RETURN
END SUBROUTINE proces
