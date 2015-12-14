SUBROUTINE dpse4
     
!     PRESSURE STIFFNESS CALCULATIONS FOR A QUADRILATERAL MEMBRANE
!     ELEMENT, WHICH HAS 4 GRID POINTS.
!     THREE 6X6 STIFFNESS MATRICES FOR THE PIVOT POINT ARE INSERTED.
 
!     DOUBLE PRECISION VERSION
 
!     WRITTEN BY E. R. CHRISTENSEN/SVERDRUP,  9/91,  VERSION 1.1
!     INSTALLED IN NASTRAN AS ELEMENT DPSE4 BY G.CHAN/UNISYS, 2/92
 
!     REFERENCE - E. CHRISTENEN: 'ADVACED SOLID ROCKET MOTOR (ASRM)
!                 MATH MODELS - PRESSURE STIFFNESS EFFECTS ANALYSIS',
!                 NASA TD 612-001-02, AUGUST 1991
 
!     LIMITATION -
!     (1) ALL GRID POINTS USED BY ANY OF THE CPSE2/3/4 ELEMENTS MUST BE
!         IN BASIC COORDINATE SYSTEM!!!
!     (2) CONSTANT PRESSURE APPLIED OVER AN ENCLOSED VOLUMN ENCOMPASSED
!         BY THE CPSE2/3/4 ELEMENTRS
!     (3) PRESSURE ACTS NORMALLY TO THE CPSE2/3/4 SURFACES
 
!     SEE NASTRAN DEMONSTRATION PROBLEM -  T13022A
 
 DOUBLE PRECISION :: gamma,kij,dp,c,SIGN
 DIMENSION        necpt(6)
!     COMMON /SYSTEM/  IBUF,NOUT
 COMMON /ds1aaa/  npvt,icstm,ncstm
 COMMON /ds1aet/  ecpt(26),dum2(2),dum12(12)
 COMMON /ds1adp/  gamma,kij(36),dp(26),c(12),SIGN(3),nk(3),ik(3)
 EQUIVALENCE      (necpt(1),ecpt(1))
 
!     ECPT FOR THE PRESSURE STIFFNESS CPES4 ELEMENT
 
!     ECPT( 1) = ELEMENT ID
!     ECPT( 2) = SIL FOR GRID POINT A OR 1
!     ECPT( 3) = SIL FOR GRID POINT B OR 2
!     ECPT( 4) = SIL FOR GRID POINT C OR 3
!     ECPT( 5) = SIL FOR GRID POINT C OR 4
!     ECPT( 6) = PRESSURE
!     ECPT( 7) = NOT USED
!     ECPT( 8) = NOT USED
!     ECPT( 9) = NOT USED
!     ECPT(10) = COORD. SYSTEM ID 1
!     ECPT(11) = X1
!     ECPT(12) = Y1
!     ECPT(13) = Z1
!     ECPT(14) = COORD. SYSTEM ID 2
!     ECPT(15) = X2
!     ECPT(16) = Y2
!     ECPT(17) = Z2
!     ECPT(18) = COORD. SYSTEM ID 3
!     ECPT(19) = X3
!     ECPT(20) = Y3
!     ECPT(21) = Z3
!     ECPT(22) = COORD. SYSTEM ID 4
!     ECPT(23) = X4
!     ECPT(24) = Y4
!     ECPT(25) = Z4
!     ECPT(26) = ELEMENT TEMPERATURE
!     ECPT(27) THRU ECPT(40) = DUM2 AND DUM12, NOT USED IN THIS ROUTINE
 
!     STORE ECPT IN DOUBLE PRECISION
 
 dp(6) = ecpt(6)
 k = 10
 DO  i = 1,4
   DO  j = 1,3
     k = k + 1
     dp(k) = ecpt(k)
   END DO
   k = k + 1
 END DO
 
!     CALCULATE THE FOUR VECTORS GAB, GAC, GAD, AND GBD USED IN
!     COMPUTING THE PRESSURE STIFFNESS MATRIC
 
!     GAB = RA + RB - RC - RD
!     GAC = RB - RD
!     GAD =-RA + RB + RC - RD
!     GBD =-RA + RC
 
!     GAB STORED IN C( 1), C( 2), C( 3)
!     GAC STORED IN C( 4), C( 5), C( 6)
!     GAD STORED IN C( 7), C( 8), C( 9)
!     GBD STORED IN C(10), C(11), C(12)
 
 c(1) = dp(11) + dp(15) - dp(19) - dp(23)
 c(2) = dp(12) + dp(16) - dp(20) - dp(24)
 c(3) = dp(13) + dp(17) - dp(21) - dp(25)
 
 c(4) = dp(15) - dp(23)
 c(5) = dp(16) - dp(24)
 c(6) = dp(17) - dp(25)
 
 c(7) =-dp(11) + dp(15) + dp(19) - dp(23)
 c(8) =-dp(12) + dp(16) + dp(20) - dp(24)
 c(9) =-dp(13) + dp(17) + dp(21) - dp(25)
 
 c(10)=-dp(11) + dp(19)
 c(11)=-dp(12) + dp(20)
 c(12)=-dp(13) + dp(21)
 
 DO  i = 1,4
   IF (necpt(i+1) /= npvt) CYCLE
   npivot = i
   GO TO 40
 END DO
 RETURN
 
!     GENERATE THE THREE BY THREE PARTITIONS IN GLOBAL COORDINATES HERE
 
!     SET COUNTERS ACCORDING TO WHICH GRID POINT IS THE PIVOT
 
 40 IF (npivot == 4) GO TO 80
 IF (npivot-2 < 0) THEN
   GO TO    50
 ELSE IF (npivot-2 == 0) THEN
   GO TO    60
 ELSE
   GO TO    70
 END IF
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KAB, KAC, KAD
 
 50 nk(1) = 2
 nk(2) = 3
 nk(3) = 4
 ik(1) = 1
 ik(2) = 4
 ik(3) = 7
 SIGN(1) = 1.0D0
 SIGN(2) = 1.0D0
 SIGN(3) = 1.0D0
 GO TO 90
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KBA, KBC, KBD
!     NOTE THAT KBA = -KAB
 
 60 nk(1) = 1
 nk(2) = 3
 nk(3) = 4
 ik(1) = 1
 ik(2) = 7
 ik(3) = 10
 SIGN(1) =-1.0D0
 SIGN(2) = 1.0D0
 SIGN(3) = 1.0D0
 GO TO 90
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KCA, KCB, KCD
!     NOTE THAT KCA = -KAC, KCB = -KBC
 
 70 nk(1) = 1
 nk(2) = 2
 nk(3) = 4
 ik(1) = 4
 ik(2) = 7
 ik(3) = 1
 SIGN(1) =-1.0D0
 SIGN(2) =-1.0D0
 SIGN(3) =-1.0D0
 GO TO 90
 
 80 nk(1) = 1
 nk(2) = 2
 nk(3) = 3
 ik(1) = 7
 ik(2) = 10
 ik(3) = 1
 SIGN(1) =-1.0D0
 SIGN(2) =-1.0D0
 SIGN(3) = 1.0D0
 
 90 gamma =-dp(6)/12.0D0
 DO  i = 1,3
   DO  j = 1,36
     kij(j) = 0.0D0
   END DO
   k1 = ik(i)
   k2 = k1 + 1
   k3 = k1 + 2
   sg = gamma*SIGN(i)
   kij( 2) =-c(k3)*sg
   kij( 3) = c(k2)*sg
   kij( 7) = c(k3)*sg
   kij( 9) =-c(k1)*sg
   kij(13) =-c(k2)*sg
   kij(14) = c(k1)*sg
   
!     ASSEMBLE INTO THE GLOBAL STIFFNESS MATRIX
   
   ias = nk(i)
   CALL ds1b (kij(1),necpt(ias+1))
 END DO
 
 RETURN
END SUBROUTINE dpse4
