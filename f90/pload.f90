SUBROUTINE pload
     
 INTEGER :: NAME(2),gridp,pont
 DIMENSION       gridp(5),igpco(4,4),gpco1(3),gpco2(3),gpco3(3),  &
     pont(4),iord(4),vect(3),vect1(3),vect2(3), ploads(3,4),gpco4(3),vect3(3)
 COMMON /loadx / lcore,slt,bgpdt,OLD,cstm,nn(11),nobld
 COMMON /system/ ksys(87),ksys88
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (pmag,gridp(1)),  &
     (igpco(2,1),gpco1(1)),(igpco(2,2),gpco2(1)),  &
     (igpco(2,3),gpco3(1)),(igpco(2,4),gpco4(1))
 DATA    NAME  / 4HPLOA,4HD   /, pi / 3.141592654 /
 
 
 DO  i = 1,3
   ploads(i,4) = 0.0
 END DO
 CALL READ (*150,*150,slt,gridp(1),5,0,flag)
 pont(1) = gridp(2)
 pont(2) = gridp(3)
 pont(3) = gridp(4)
 pont(4) = gridp(5)
 n1 = 4
 IF (gridp(5) == 0) n1 = 3
 CALL permut (pont(1),iord(1),n1,OLD)
 DO  i = 1,n1
   l = iord(i)
   CALL fndpnt (igpco(1,l),pont(l))
 END DO
 IF (n1 == 4) GO TO 160
 
!     THREE  POINTS
 
 DO  i = 1,3
   vect3(i) = gpco1(i) - gpco2(i)
   vect2(i) = gpco3(i) - gpco1(i)
   vect1(i) = gpco2(i) - gpco3(i)
 END DO
 CALL cross (vect3(1),vect1(1),vect(1))
 
 DO  i = 1,3
   DO  j = 1,3
     ploads(j,i) = -vect(j)
   END DO
 END DO
 
 IF (ksys88 == 1) GO TO 50
 
!     KSYS88 = 0, PRESSURE LOAD IS DISTRIBUTED EVENLY (ONE-THIRD) TO
!     EACH OF THE 3 GRID POINTS. TRIANGULAR GEOMETRY IS NOT CONSIDERED.
 
 pmag    = pmag/6
 vect(1) = pmag
 vect(2) = pmag
 vect(3) = pmag
 GO TO 80
 
!     IMPLEMENTED BY G.CHAN/UNISYS   3/1990
!     KSYS88 = 1, PRESSURE LOAD IS DISTRIBUTED PROPORTIONALLY TO THE
!     THREE ANGLE SIZES.
!     E.G. A 45-90-45 DEGREE TRIANGLE ELEMENT WILL HAVE TWICE THE LOAD
!     AT THE 90 DEGREE ANGLE TO THAT OF THE 45 DEGREE ANGLE.
!     RECTANGULAR ELEMENT (4 POINTS) IS NOT AFFECTED
 
!     GET AREA(2X), SIDES (VI) AND ANGLES (AI) OF THE TRIANGLE
 
 50 CONTINUE
 area = SQRT(vect (1)**2 + vect (2)**2 + vect (3)**2)
 v1   = SQRT(vect1(1)**2 + vect1(2)**2 + vect1(3)**2)
 v2   = SQRT(vect2(1)**2 + vect2(2)**2 + vect2(3)**2)
 v3   = SQRT(vect3(1)**2 + vect3(2)**2 + vect3(3)**2)
 
!     CHOOSE AN ANGLE, WHICH IS NOT THE LARGEST, TO START COMPUTING
!     THE THREE ANGLES
 
 IF (v2 > v1 .AND. v2 > v3) GO TO 60
 sin2 = area/(v3*v1)
 sin1 = v1*sin2/v2
 sin3 = v3*sin2/v2
 a2   = ASIN(sin2)
 IF (sin1 >= 0.0) a1 = ASIN(sin1)
 IF (sin3 >= 0.0) a3 = ASIN(sin3)
 IF (v1 > v3) a1 = pi - a2 - a3
 IF (v3 > v1) a3 = pi - a2 - a1
 GO TO 70
 
 60 sin3 = area/(v2*v1)
 sin2 = v2*sin3/v3
 sin1 = v1*sin3/v3
 a3   = ASIN(sin3)
 IF (sin2 >= 0.0) a2 = ASIN(sin2)
 IF (sin1 >= 0.0) a1 = ASIN(sin1)
 IF (v1 > v2) a1 = pi - a3 - a2
 IF (v2 > v1) a2 = pi - a3 - a1
 70 pmag    = 0.5*pmag/pi
 vect(1) = pmag*a1
 vect(2) = pmag*a2
 vect(3) = pmag*a3
 
!     TRANSFORM TO GLOBAL AND ADD CONTRIBUTIONS
 
 80 DO  i = 1,n1
   DO   j = 1,3
     IF (n1 == 4) ploads(j,i) = -ploads(j,i)*pmag
     IF (n1 == 3) ploads(j,i) = -ploads(j,i)*vect(i)
   END DO
   IF (igpco(1,i) /= 0) CALL basglb (ploads(1,i),ploads(1,i),  &
       igpco(2,i),igpco(1,i))
   CALL fndsil (pont(i))
   DO  j = 1,3
     in = pont(i) + j - 1
     core(in) = ploads(j,i) + core(in)
   END DO
 END DO
 140 RETURN
 
 150 CALL mesage (-1,slt,NAME)
 GO TO 140
 
!     FOUR  POINTS
 
 
!     TRIANGLE  1,2,3
 
 160 DO  i = 1,3
   vect1(i) = gpco1(i) - gpco2(i)
   vect2(i) = gpco3(i) - gpco2(i)
 END DO
 CALL cross (vect1(1),vect2(1),vect(1))
 DO  i = 1,3
   DO  j = 1,3
     ploads(j,i) = vect(j)
   END DO
 END DO
 
!     TRIANGLE  2,3,4
 
 DO  i  =1,3
   vect1(i) = gpco2(i) - gpco3(i)
   vect2(i) = gpco4(i) - gpco3(i)
 END DO
 CALL cross (vect1(1),vect2(1),vect(1))
 DO  i = 2,4
   DO  j = 1,3
     ploads(j,i) = ploads(j,i) + vect(j)
   END DO
 END DO
 
!     TRIANGLE  3,1,4
 
 DO  i = 1,3
   vect1(i) = gpco4(i) - gpco1(i)
   vect2(i) = gpco3(i) - gpco1(i)
 END DO
 CALL cross (vect1(1),vect2(1),vect(1))
 DO  i = 1,4
   IF (i == 2) CYCLE
   DO  j = 1,3
     ploads(j,i) = ploads(j,i)+vect(j)
   END DO
 END DO
 
!     TRIANGLE (4,1,2)
 
 DO  i = 1,3
   vect1(i) = gpco4(i) - gpco1(i)
   vect2(i) = gpco2(i) - gpco1(i)
 END DO
 CALL cross (vect1(1),vect2(1),vect(1))
 DO  i = 1,4
   IF (i == 3) CYCLE
   DO  j = 1,3
     ploads(j,i) = ploads(j,i) + vect(j)
   END DO
 END DO
 pmag = pmag/12.0
 GO TO 80
END SUBROUTINE pload
