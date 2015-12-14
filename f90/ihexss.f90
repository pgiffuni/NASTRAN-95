SUBROUTINE ihexss (TYPE,shp,dshp,jacob,detj,eid,xi,eta,zeta,bxyz)
     
!     SINGLE PRECISION VERSION
 
!     ISOPARAMETRIC UTILITY ROUTINE.  THIS ROUTINE WILL COMPUTE
!     VALUES OF THE SHAPE FUNCTIONS, THEIR DERIVATIVES WITH RESPECT TO
!     XI,ETA, AND ZETA, THE JACOBIAN MATRIX INVERSE, AND ITS DETERMINANT
 
!                       TYPE = 1       IHEX1
!                       TYPE = 2       IHEX2
!                       TYPE = 3       IHEX3
 
!     SHP    = VALUES OF SHAPE FUNCTIONS
!     DSHP   = DERIVATIVES OF SHAPE FUNCTIONS W.R.T. XI, ETA, ZETA
!     JACOB  = JACOBIAN MATRIX INVERSE
!     DETJ   = DETERMINANT OF JACOBIAN MATRIX
!     XI, ETA, ZETA = ELEMENT COORDINATES AT WHICH THESE COMPUTATIONS
!                     TAKE PLACE
!     BXYZ   = BASIC SYSTEM COORDINATES FOR GRID POINTS
 
!     LOCAL VARIABLES
!     X,Y,Z  = CONSTANTS FOR EACH SHAPE FUNCTION
!     NGP    = NUMBER OF SHAPE FUNCTIONS, ALSO NUMBER OF GRID POINTS
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 REAL, INTENT(OUT)                        :: shp(8)
 REAL, INTENT(OUT)                        :: dshp(3,8)
 REAL, INTENT(OUT)                        :: jacob(3,3)
 REAL, INTENT(OUT)                        :: detj
 INTEGER, INTENT(IN OUT)                  :: eid
 REAL, INTENT(IN)                         :: xi
 REAL, INTENT(IN)                         :: eta
 REAL, INTENT(IN)                         :: zeta
 REAL, INTENT(IN)                         :: bxyz(3,8)
 INTEGER :: op
 REAL :: work(3,3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf     ,op
 
 ngp = 12*TYPE - 4
 y   =-1.0
 z   =-1.0
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 310
 END SELECT
 
!     LINEAR ELEMENT IHEX1
 
 100 DO  j = 1,2
   IF (j == 2) z = 1.0
   x = -1.0
   y = -1.0
   DO  i = 1,4
     IF (i == 3) y = 1.0
     IF (i == 2) x = 1.0
     IF (i == 4) x =-1.0
     k    = i + (j-1)*4
     qxi  = 1.0 + xi*x
     qeta = 1.0 + eta*y
     qzeta= 1.0 + zeta*z
     shp(k) = qxi*qeta*qzeta/8.0
     dshp(1,k) = x*qeta*qzeta/8.0
     dshp(2,k) = y*qxi*qzeta/8.0
     dshp(3,k) = z*qxi*qeta/8.0
   END DO
 END DO
 GO TO 430
 
!     QUADRATIC ELEMENT IHEX2
 
 200 d = 1.0
 x = 0.0
 DO  i = 1,20
!            1   2   3   4   5   6   7   8   9   10
   SELECT CASE ( i )
     CASE (    1)
       GO TO 220
     CASE (    2)
       GO TO 210
     CASE (    3)
       GO TO 210
     CASE (    4)
       GO TO 230
     CASE (    5)
       GO TO 230
     CASE (    6)
       GO TO 220
     CASE (    7)
       GO TO 220
     CASE (    8)
       GO TO 240
     CASE (    9)
       GO TO 250
     CASE (   10)
       GO TO 210
     CASE (   11)
       GO TO 230
     CASE (   12)
       GO TO 220
     CASE (   13)
       GO TO 250
     CASE (   14)
       GO TO 210
     CASE (   15)
       GO TO 210
     CASE (   16)
       GO TO 230
     CASE (   17)
       GO TO 230
     CASE (   18)
       GO TO 220
     CASE (   19)
       GO TO 220
     CASE (   20)
       GO TO 240
   END SELECT
   210 x = x+d
   GO TO 260
   220 x = x - d
   GO TO 260
   230 y = y + d
   GO TO 260
   240 y = y - d
   GO TO 260
   250 z = z + 1.0
   y =-1.0
   d = 3.0-d
   260 IF (x == 0.0) GO TO 270
   IF (y == 0.0) GO TO 280
   IF (z == 0.0) GO TO 290
   
!     CORNER POINT
   
   qxi   = 1.0 + x*xi
   qeta  = 1.0 + y*eta
   qzeta = 1.0 + z*zeta
   qxyz  = x*xi + y*eta+z*zeta
   shp(i)= qxi*qeta*qzeta*(qxyz-2.0)/8.0
   dshp(1,i) = x*qeta*qzeta*(x*xi+qxyz-1.0)/8.0
   dshp(2,i) = y*qxi*qzeta*(y*eta+qxyz-1.0)/8.0
   dshp(3,i) = z*qxi*qeta*(z*zeta+qxyz-1.0)/8.0
   CYCLE
   
!     MID-EDGE POINT, X = 0.0
   
   270 qxi   = 1.0 - xi**2
   qeta  = 1.0 + y*eta
   qzeta = 1.0 + z*zeta
   shp(i)= qxi*qeta*qzeta/4.0
   dshp(1,i) =-xi*qeta*qzeta/2.0
   dshp(2,i) = qxi*qzeta*y/4.0
   dshp(3,i) = qxi*qeta*z/4.0
   CYCLE
   
!     MID-EDGE POINT, Y = 0.0
   
   280 qxi   = 1.0 + x*xi
   qeta  = 1.0 - eta**2
   qzeta = 1.0 + z*zeta
   shp(i)= qeta*qxi*qzeta/4.0
   dshp(1,i) = qeta*qzeta*x/4.0
   dshp(2,i) =-eta*qzeta*qxi/2.0
   dshp(3,i) = qeta*qxi*z/4.0
   CYCLE
   
!     MID-EDGE POINT, Z = 0.0
   
   290 qxi   = 1.0 + x*xi
   qeta  = 1.0 + y*eta
   qzeta = 1.0 - zeta**2
   shp(i)= qzeta*qxi*qeta/4.0
   dshp(1,i) = qzeta*qeta*x/4.0
   dshp(2,i) = qzeta*qxi*y/4.0
   dshp(3,i) =-zeta*qxi*qeta/2.0
 END DO
 GO TO 430
 
!     CUBIC ELEMENT IHEX3
 
 310 d = 2.0/3.0
 x =-1.0/3.0
 DO  i = 1,32
!            1   2   3   4   5   6   7   8   9   10
   GO TO (320,330,330,330,340,340,340,320,320,320,  &
       350,350,360,330,340,320,360,330,340,320,  &
       360,330,330,330,340,340,340,320,320,320, 350,350), i
   320 x = x - d
   GO TO 370
   330 x = x + d
   GO TO 370
   340 y = y + d
   GO TO 370
   350 y = y - d
   GO TO 370
   360 y =-1.0
   z = z + 2.0/3.0
   IF (z > -1.0) d = 2.0
   IF (z >  0.4) d = 2.0/3.0
   370 IF (ABS(x) < 0.4) GO TO 390
   IF (ABS(y) < 0.4) GO TO 400
   IF (ABS(z) < 0.4) GO TO 410
   
!     CORNER POINT
   
   qxi   = 1.0 + x*xi
   qeta  = 1.0 + y*eta
   qzeta = 1.0 + z*zeta
   qxyz  = xi**2 + eta**2 + zeta**2 - 19.0/9.0
   shp(i)= 9.0*qxi*qeta*qzeta*qxyz/64.0
   dshp(1,i) = 9.0*qeta*qzeta*(x*(2.0*xi**2+qxyz)+2.0*xi)/64.0
   dshp(2,i) = 9.0*qxi*qzeta*(y*(2.0*eta**2+qxyz)+2.0*eta)/64.0
   dshp(3,i) = 9.0*qxi*qeta*(z*(2.0*zeta**2+qxyz)+2.0*zeta)/64.0
   CYCLE
   
!     MID-EDGE POINT, X = + OR - 1/3
   
   390 qxi   = 9.0*(1.0-xi**2)*(1.0+9.0*x*xi)/64.0
   qeta  = 1.0 + y*eta
   qzeta = 1.0 + z*zeta
   qxyz  = 9.0*(-2.0*xi+9.0*x-27.0*xi*x*xi)/64.0
   shp(i)= qxi*qeta*qzeta
   dshp(1,i) = qeta*qzeta*qxyz
   dshp(2,i) = qxi*qzeta*y
   dshp(3,i) = qxi*qeta*z
   CYCLE
   
!     MID-EDGE POINT Y =+ OR - 1/3
   
   400 qxi   = 1.0 + x*xi
   qeta  = 9.0*(1.0-eta**2)*(1.0+9.0*eta*y)/64.0
   qzeta = 1.0 + z*zeta
   qxyz  = 9.0*(-2.0*eta+9.0*y-27.0*eta*y*eta)/64.0
   shp(i)= qeta*qxi*qzeta
   dshp(1,i) = qeta*qzeta*x
   dshp(2,i) = qxi*qzeta*qxyz
   dshp(3,i) = qeta*qxi*z
   CYCLE
   
!     MID-EDGE POINTS Z =+ OR - 1/3
   
   410 qxi   = 1.0 + x*xi
   qeta  = 1.0 + y*eta
   qzeta = 9.0*(1.0-zeta**2)*(1.0+9.0*z*zeta)/64.0
   qxyz  = 9.0*(-2.0*zeta+9.0*z-27.0*z*zeta**2)/64.0
   shp(i)= qzeta*qxi*qeta
   dshp(1,i) = qzeta*qeta*x
   dshp(2,i) = qzeta*qxi*y
   dshp(3,i) = qxi*qeta*qxyz
 END DO
 
!     COMPUTE JACOBIAN MATRIX
 
 430 DO  i = 1,3
   DO  j = 1,3
     jacob(i,j) = 0.0
     DO  k = 1,ngp
       jacob(i,j) = jacob(i,j)+dshp(i,k)*bxyz(j,k)
     END DO
   END DO
 END DO
 
!     COMPUTE INVERSE AND DETERMINANT OF JACOBIAN MATRIX
 
 work(1,1) = jacob(2,2)*jacob(3,3) - jacob(2,3)*jacob(3,2)
 work(2,1) = jacob(2,3)*jacob(3,1) - jacob(2,1)*jacob(3,3)
 work(3,1) = jacob(2,1)*jacob(3,2) - jacob(2,2)*jacob(3,1)
 work(1,2) = jacob(1,3)*jacob(3,2) - jacob(1,2)*jacob(3,3)
 work(2,2) = jacob(1,1)*jacob(3,3) - jacob(1,3)*jacob(3,1)
 work(3,2) = jacob(1,2)*jacob(3,1) - jacob(1,1)*jacob(3,2)
 work(1,3) = jacob(1,2)*jacob(2,3) - jacob(1,3)*jacob(2,2)
 work(2,3) = jacob(1,3)*jacob(2,1) - jacob(1,1)*jacob(2,3)
 work(3,3) = jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1)
 detj = 0.0
 DO  i = 1,3
   detj = detj + jacob(i,2)*work(2,i)
 END DO
 IF (detj == 0.0) GO TO 470
 DO  i = 1,3
   DO  j = 1,3
     jacob(i,j) = work(i,j)/detj
   END DO
 END DO
 RETURN
 
!     JACOBIAN MATRIX WAS SINGULAR.
 
 470 WRITE  (op,480) ufm,eid
 480 FORMAT (a23,' 3306, SINGULAR JACOBIAN MATRIX FOR ISOPARAMETRIC ',  &
     'ELEMENT NO.',i9)
 RETURN
END SUBROUTINE ihexss
