SUBROUTINE deltkl (del,r,z,kode)
     
!     EVAULATE -
!        DELT(K,L) = SURFACE-INTEGRAL((R**K)*(Z**L)) DR*DZ
!        WHERE  DR*DZ IS EITHER A TRIANGLE OR A TRAPEZOID.
 
!     USAGE -
!        CALL  DELTKL (DEL,R,Z,KODE)
!        WHERE  DEL   =  DOUBLE PRECISION ARRAY OF 15 LOCATIONS.
!                        CONTAINING THE RESULTS.
!        WHERE  R     =  DOUBLE PRECISION ARRAY OF  4 LOCATIONS.
!                        CONTAINING THE R-COORDINATES OF THE ELEMENT
!        WHERE  Z     =  DOUBLE PRECISION ARRAY OF  4 LOCATIONS.
!                        CONTAINING THE Z-COORDINATES OF THE ELEMENT
!               KODE  =  0  FOR TRIANGULAR  ELEMENT
!               KODE  =  1  FOR TRAPEZOIDAL ELEMENT
 
!     PROCEDURE -
!        INFORMATION IS COMPUTED AND STORED AS FOLLOWS.
 
!           COMPUTED  FOR        ELEMENT         STORED
!        TRIANGLE   TRAPEZOID    DELT(K,L)      DEL(LOC)
!        ================================================
!           X          X              0,0            1
!           X          X              1,0            2
!           X          X              0,1            3
!           X          X             -1,0            4
!           X          X             -1,1            5
!           X          X             -1,2            6
!                      X              1,1            7
!                      X              1,2            8
!                      X              2,1            9
!                      X              2,0           10
!                      X              0,2           11
!                      X              3,0           12
!                      X              3,1           13
!                      X              3,2           14
!                      X              2,2           15
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: del(15)
 DOUBLE PRECISION, INTENT(IN)             :: r(4)
 DOUBLE PRECISION, INTENT(IN)             :: z(4)
 INTEGER, INTENT(IN OUT)                  :: kode
 INTEGER :: goback, n,i,m,l
 DOUBLE PRECISION :: rm,ri,zm,zi,ln,a,b,si,sm,  &
     r1,r2,r3,r4,r5,ab,a2,b2,a3,b3,aab,abb,r6,r7
 
!     ZERO ARRAY (ONLY THAT PORTION USING)
 
 n = 15
 DO  l = 1,n
   del(l) = 0.0D+0
 END DO
 
!     HERE FOR LINE 1-2
 
 i = 1
 m = 2
 ASSIGN  23 TO GO back
 GO TO 50
 
!     HERE FOR LINE 2-3
 
 23 CONTINUE
 i = 2
 m = 3
 ASSIGN 34 TO GO back
 GO TO 50
 
!     HERE FOR LINE 31 (TRIANGLE),  LINE 3-4 (TRAP)
 
 34 CONTINUE
 i = 3
 IF (kode > 0) GO TO 35
 m = 1
 ASSIGN 90 TO GO back
 GO TO 50
 35 m = 4
 ASSIGN 41 TO GO back
 GO TO 50
 41 i = 4
 m = 1
 ASSIGN 90 TO GO back
 
!     BEGIN LOCAL SUBROUTINE (DEL-KL-I,M)
 
 50 rm = r(m)
 ri = r(i)
 r1 = rm - ri
 IF (DABS(r1) < 1.0D-07) GO TO 80
 
!     THIS LINE IS NOT PARALLEL TO Z-AXIS
 
 zm = z(m)
 zi = z(i)
 IF (zi == 0.0D+0 .AND. zm == 0.0D+0) GO TO 80
 
!     SPECIAL CASE, ZM = ZI = 0   THUS ALL  A,B = 0  AND
!     ALL DEL TERMS  = 0 .   THUS SKIP AND SAVE CPU.
 
 a   = (rm*zi - ri*zm)/r1
 b   = (zm - zi)/r1
 ln  = DLOG(rm/ri)
 si  = ri * ri
 sm  = rm * rm
 r2  = sm - si
 si  = si * ri
 sm  = sm * rm
 r3  = sm - si
 si  = si * ri
 sm  = sm * rm
 r4  = sm - si
 si  = si * ri
 sm  = sm * rm
 r5  = sm - si
 a2  = a  * a
 a3  = a  * a2
 b2  = b  * b
 b3  = b  * b2
 ab  = a  * b
 aab = a  * ab
 abb = b  * ab
 del( 1) = a*r1 + b*r2/2.0D+0 + del(1)
 del( 2) = a*r2/  2.0D+0 + b*r3 /3.0D+0 + del(2)
 del( 3) = a2*r1/ 2.0D+0 + ab*r2/2.0D+0 + b2*r3/6.0D+0 + del(3)
 del( 4) = a*ln + b*r1   + del(4)
 del( 5) = a2*ln/ 2.0D+0 + ab*r1  + b2*r2 /4.0D+0 + del(5)
 del( 6) = a3*ln/ 3.0D+0 + aab*r1 + abb*r2/2.0D+0 + b3*r3/9.0D+0 + del(6)
 del( 7) = a2*r2/ 4.0D+0 + ab*r3/3.0D+0  + b2*r4 /8.0D+0 + del(7)
 del( 8) = a3*r2/ 6.0D+0 + aab*r3/3.0D+0 + abb*r4/4.0D+0  &
     + b3*r5/15.0D+0 + del(8)
 del( 9) = a2*r3/ 6.0D+0 + ab*r4/4.0D+0 + b2*r5/10.0D+0 + del(9)
 del(10) = a *r3/ 3.0D+0 + b*r4 /4.0D+0 + del(10)
 del(12) = a *r4/ 4.0D+0 + b*r5 /5.0D+0 + del(12)
 IF (kode < 1) GO TO 80
 si      = si*ri
 sm      = sm*rm
 r6      = sm - si
 r7      = (sm*rm - si*ri)
 del(11) = a3*r1/ 3.0D+0 + aab*r2/2.0D+0 + abb*r3/3.0D+0  &
     + b3*r4/12.0D+0 + del(11)
 del(13) = a2*r4/ 8.0D+0 + ab*r5 /5.0D+0 + b2*r6/12.0D+0 + del(13)
 del(14) = a3*r4/12.0D+0 + aab*r5/5.0D+0 + abb*r6/6.0D+0  &
     + b3*r7/21.0D+0 + del(14)
 del(15) = a3*r3/ 9.0D+0 + aab*r4/4.0D+0 + abb*r5/5.0D+0  &
     + b3*r6/18.0D+0 + del(15)
 80 GO TO GO back, (23,34,41,90)
 
!     THE ABSOLUTE VALUE IS CHOSEN SO THAT NODES INPUT MAY BE ORDERED
!     CW OR CCW.   RESULTS ARE SAME FOR A GIVEN ELEMENT.
 
 90 DO  l = 1,n
   del(l) = DABS(del(l))
 END DO
 RETURN
END SUBROUTINE deltkl
