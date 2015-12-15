SUBROUTINE delkls (del,r,z,kode)
     
!     SINGLE PRECISION VERSION USE DELKLS (DEL,R,Z,KODE)
 
!     PURPOSE
!        EVAULATE THE FOLLOWING FUNCTION
!          DELT(K,L) = SURFACE-INTEGRAL((R**K)*(Z**L)) DR*DZ
!        WHERE  DR*DZ  IS EITHER A TRIANGLE OR A TRAPEZOID.
 
!     USAGE
!        WHERE  DEL   =  DOUBLE PRECISION ARRAY OF 15 LOCATIONS.
!                        CONTAINING THE RESULTS.
!        WHERE  R     =  DOUBLE PRECISION ARRAY OF  4 LOCATIONS.
!                        CONTAINING THE R-COORDINATES OF THE ELEM.
!        WHERE  Z     =  DOUBLE PRECISION ARRAY OF  4 LOCATIONS.
!                        CONTAINING THE Z-COORDINATES OF THE ELEM.
!               KODE  =  0  FOR TRIANGULAR  ELEMENT
!               KODE  =  1  FOR TRAPEZOIDAL ELEMENT
 
!     PROCEDURE
!        INFORMATION IS COMPUTED AND STORED AS FOLLOWS.
!           COMPUTED  FOR        ELEMENT         STORED
!        TRIANGLE   TRAPEZOID   DELT(K,L)       DEL(LOC)
!        ************************************************
!           X          X              0,0            01
!           X          X              1,0            02
!           X          X              0,1            03
!           X          X             -1,0            04
!           X          X             -1,1            05
!           X          X             -1,2            06
!                      X              1,1            07
!                      X              1,2            08
!                      X              2,1            09
!                      X              2,0            10
!                      X              0,2            11
!                      X              3,0            12
!                      X              3,1            13
!                      X              3,2            14
!                      X              2,2            15
 
 
 REAL, INTENT(OUT)                        :: del(15)
 REAL, INTENT(IN)                         :: r(4)
 REAL, INTENT(IN)                         :: z(4)
 INTEGER, INTENT(IN OUT)                  :: kode
 INTEGER :: go_back
 REAL :: ln
 
!     ZERO ARRAY (ONLY THAT PORTION USING)
 
 n      = 15
 DO  l = 1,n
   del(l) = 0.
 END DO
 
!     HERE FOR  LINE 1-2
 
 i      = 1
 m      = 2
 ASSIGN 23 TO go_back
 GO TO 500
 
!     HERE  FOR  LINE 2-3
 
 23 CONTINUE
 i      = 2
 m      = 3
 ASSIGN 3134 TO go_back
 GO TO 500
 
!     HERE  FOR  LINE 31 (TRIANGLE),  LINE 3-4 (TRAP)
 
 3134 CONTINUE
 i      = 3
 IF (kode > 0) GO TO 35
 m      = 1
 ASSIGN 90 TO go_back
 GO TO 500
 35 m      = 4
 ASSIGN 41 TO go_back
 GO TO 500
 41 i      = 4
 m      = 1
 ASSIGN 90 TO go_back
 
!     BEGIN LOCAL SUBROUTINE  (DEL-KL-I,M)
 
 500 rm     = r(m)
 ri     = r(i)
 r1     = rm - ri
 IF (ABS(r1) < 1.e-5) GO TO 599
 
!     THIS LINE IS NOT PARALLEL TO Z-AXIS
 
 zm     = z(m)
 zi     = z(i)
 IF (zi == 0. .AND. zm == 0.) GO TO 599
 
!     SPECIAL CASE   ZM=ZI=0   THUS ALL  A,B = 0  AND
!     ALL DEL TERMS  = 0 .     THUS SKIP AND SAVE CPU.
 
 a      = (rm*zi - ri*zm)/r1
 b      = (zm - zi)/r1
 ln     = ALOG(rm/ri)
 si     = ri * ri
 sm     = rm * rm
 r2     = sm - si
 si     = si * ri
 sm     = sm * rm
 r3     = sm - si
 si     = si * ri
 sm     = sm * rm
 r4     = sm - si
 si     = si * ri
 sm     = sm * rm
 r5     = sm - si
 a2     = a  * a
 a3     = a  * a2
 b2     = b  * b
 b3     = b  * b2
 ab     = a  * b
 aab    = a  * ab
 abb    = b  * ab
 del(1) = a*r1 + b*r2/2. + del(1)
 del(2) = a*r2/2.  +  b*r3/3.  + del(2)
 del(3) = a2*r1/2. + ab*r2 /2. + b2*r3/6.  + del(3)
 del(4) = a *ln    +  b*r1     + del(04)
 del(5) = a2*ln/2. + ab*r1     + b2*r2/4.  + del(5)
 del(6) = a3*ln/3. + aab*r1    + abb*r2/2. + b3*r3/9.  + del(6)
 del(7) = a2*r2/4. + ab*r3/3.  + b2*r4/8.  + del(7)
 del(8) = a3*r2/6. + aab*r3/3. + abb*r4/4. + b3*r5/15. + del(8)
 del(9) = a2*r3/6. + ab*r4/4.  + b2*r5/10. + del(9)
 del(10)= a*r3/3.  +  b*r4/4.  + del(10)
 del(12)= a*r4/4.  + b*r5/5.   + del(12)
 IF (kode < 1) GO TO 599
 si     = si * ri
 sm     = sm * rm
 r6     = sm - si
 r7     = (sm*rm - si*ri)
 del(11)= a3*r1/3. + aab*r2/2. + abb*r3/3. + b3*r4/12. + del(11)
 del(13)= a2*r4/8. + ab*r5/5.  + b2*r6/12. + del(13)
 del(14)= a3*r4/12.+ aab*r5/5. + abb*r6/6. + b3*r7/21. + del(14)
 del(15)= a3*r3/9. + aab*r4/4. + abb*r5/5. + b3*r6/18. + del(15)
 599 GO TO go_back, (23,3134,41,90)
 
!     THE ABSOLUTE VALUE IS CHOSEN SO THAT NODES INPUT MAY BE ORDERED
!     CW OR CCW. RESULTS ARE SAME FOR A GIVEN ELEMENT.
 
 90 DO  l = 1,n
   del(l) = ABS(del(l))
 END DO
 
 RETURN
END SUBROUTINE delkls
