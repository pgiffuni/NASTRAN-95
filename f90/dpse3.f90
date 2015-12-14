SUBROUTINE dpse3
     
!     PRESSURE STIFFNESS CALCULATIONS FOR A TRIANGULAR MEMBRANE
!     ELEMENT (3 GRID POINTS).
!     THREE 6X6 STIFFNESS MATRICES FOR THE PIVOT POINT ARE INSERTED.
 
!     DOUBLE PRECISION VERSION
 
!     WRITTEN BY E. R. CHRISTENSEN/SVERDRUP, 9/91, VERSION 1.1
!     INSTALLED IN NASTRAN AS ELEMENT DPSE3 BY G.CHAN/UNISYS, 2/92
 
!     REFERENCE - E. CHRISTENEN: 'ADVACED SOLID ROCKET MOTOR (ASRM)
!                 MATH MODELS - PRESSURE STIFFNESS EFFECTS ANALYSIS',
!                 NASA TD 612-001-02, AUGUST 1991
 
!     LIMITATION -
!     (1) ALL GRID POINTS USED BY ANY IF THE CPSE2/3/4 ELEMENTS MUST BE
!         IN BASIC COORDINATE SYSTEM!!!
!     (2) CONSTANT PRESSURE APPLIED OVER AN ENCLOSED VOLUMN ENCOMPASSED
!         BY THE CPSE2/3/4 ELEMENTRS
!     (3) PRESSURE ACTS NORMALLY TO THE CPSE2/3/4 SURFACES
 
!     SEE NASTRAN DEMONSTRATION PROBLEM -  T13022A
 
 DOUBLE PRECISION :: gamma,kij,sgn,sgn1,sgn2,dp,c
 DIMENSION        necpt(5)
!     COMMON /SYSTEM/  IBUF,NOUT
 COMMON /ds1aaa/  npvt,icstm,ncstm
 COMMON /ds1aet/  ecpt(21),dum2(2),dum9(9)
 COMMON /ds1adp/  gamma,kij(36),sgn,sgn1,sgn2,dp(21),c(9)
 EQUIVALENCE      (necpt(1),ecpt(1))
 
!     ECPT FOR THE PRESSURE STIFFNESS CPES3 ELEMENT
 
!     ECPT( 1) = ELEMENT ID
!     ECPT( 2) = SIL FOR GRID POINT A OR 1
!     ECPT( 3) = SIL FOR GRID POINT B OR 2
!     ECPT( 4) = SIL FOR GRID POINT C OR 3
!     ECPT( 5) = PRESSURE
!     ECPT( 6) = NOT USED
!     ECPT( 7) = NOT USED
!     ECPT( 8) = NOT USED
!     ECPT( 9) = COORD. SYSTEM ID 1
!     ECPT(10) = X1
!     ECPT(11) = Y1
!     ECPT(12) = Z1
!     ECPT(13) = COORD. SYSTEM ID 2
!     ECPT(14) = X2
!     ECPT(15) = Y2
!     ECPT(16) = Z2
!     ECPT(17) = COORD. SYSTEM ID 3
!     ECPT(18) = X3
!     ECPT(19) = Y3
!     ECPT(20) = Z3
!     ECPT(21) = ELEMENT TEMPERATURE
!     ECPT(22) THRU (32) = DUM2 AND DUM9, NOT USED IN THIS ROUTINE
 
!     STORE ECPT IN DOUBLE PRECISION
 
 dp(5) = ecpt(5)
 k = 9
 DO  i = 1,3
   DO  j = 1,3
     k = k + 1
     dp(k) = ecpt(k)
   END DO
   k = k + 1
 END DO
 
!     CALCULATE THE THREE VECTORS R1, R2 AND R2 USED IN COMPUTING
!     THE PRESSURE STIFFNESS MATRICES:
 
!     R1 = RA - 2*RC + RB
!     R2 = 2*RB - RA - RC
!     R3 = RB - 2*RA + RC
 
!     R1 STORED IN C(1), C(2), C(3)
!     R2 STORED IN C(4), C(5), C(6)
!     R3 STORED IN C(7), C(8), C(9)
 
 c(1) = dp(10) - 2.0D0*dp(18) + dp(14)
 c(2) = dp(11) - 2.0D0*dp(19) + dp(15)
 c(3) = dp(12) - 2.0D0*dp(20) + dp(16)
 
 c(4) = 2.0D0*dp(14) - dp(10) - dp(18)
 c(5) = 2.0D0*dp(15) - dp(11) - dp(19)
 c(6) = 2.0D0*dp(16) - dp(12) - dp(20)
 
 c(7) = dp(14) - 2.0D0*dp(10) + dp(18)
 c(8) = dp(15) - 2.0D0*dp(11) + dp(19)
 c(9) = dp(16) - 2.0D0*dp(12) + dp(20)
 
 DO  i = 1,3
   IF (necpt(i+1) /= npvt) CYCLE
   npivot = i
   GO TO 40
 END DO
 RETURN
 
!     GENERATE THE THREE BY THREE PARTITIONS IN GLOBAL COORDINATES HERE
 
!     SET COUNTERS ACCORDING TO WHICH GRID POINT IS THE PIVOT
 
 40 IF (npivot-2 < 0) THEN
   GO TO    50
 ELSE IF (npivot-2 == 0) THEN
   GO TO    60
 ELSE
   GO TO    70
 END IF
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KAB, KAC
 
 50 ni = 2
 nj = 3
 nk = 1
 k1 = 1
 k2 = 4
 sgn1 = 1.0D0
 sgn2 = 1.0D0
 GO TO 80
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KBA, KBC
!     NOTE THAT KBA = -KAB
 
 60 ni = 1
 nj = 3
 nk = 2
 k1 = 1
 k2 = 7
 sgn1 =-1.0D0
 sgn2 = 1.0D0
 GO TO 80
 
!     SET COUNTERS AND POINTERS FOR CALCULATING KCA, KCB
!     NOTE THAT KCA = -KAC, KCB = -KBC
 
 70 ni = 1
 nj = 2
 nk = 1
 k1 = 4
 k2 = 7
 sgn1 =-1.0D0
 sgn2 =-1.0D0
 
 80 gamma =-dp(5)/12.0D0
 sgn   = sgn1*gamma
 k = k1
 DO  i = ni,nj,nk
   DO   j = 1,36
     kij(j) = 0.0D0
   END DO
   kk2 = k + 1
   kk3 = k + 2
   kij( 2) =-c(kk3)*sgn
   kij( 3) = c(kk2)*sgn
   kij( 7) = c(kk3)*sgn
   kij( 9) =-c(k  )*sgn
   kij(13) =-c(kk2)*sgn
   kij(14) = c(k  )*sgn
   CALL ds1b (kij(1),necpt(i+1))
   sgn = sgn2*gamma
   k   = k2
 END DO
 
 RETURN
END SUBROUTINE dpse3
