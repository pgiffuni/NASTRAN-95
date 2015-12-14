SUBROUTINE strqd1 ( ntype )
     
!**************** PHASE I  STRESS DATA RECOVERY ************************
! **********************************************************************
 
!     9/12/67         E C P T     L I S T I N G
!                    ***************************
! ECPT  TRMEM   QDMEM   TRPLT   QDPLT   TRIA1   QUAD1   TRIA2   QUAD2
! **********************************************************************
!   1   EL.ID   EL.ID   EL.ID   EL.ID   EL.ID   EL.ID   EL.ID   EL.ID
!   2   GRID A  GRID A  GRID A  GRID A  GRID A  GRID A  GRID A  GRID A
!   3   GRID B  GRID B  GRID B  GRID B  GRID B  GRID B  GRID B  GRID B
!   4   GRID C  GRID C  GRID C  GRID C  GRID C  GRID C  GRID C  GRID C
!   5   THETA   GRID D  THETA   GRID D  THETA   GRID D  THETA   GRID D
!   6   MATID   THETA   MATID1  THETA   MATID1  THETA   MAT ID  THETA
!   7   T       MAT ID  I       MATID1  T1      MATID1  T       MAT ID
!   8   NS MASS T       MATID2  I       MATID2  T1      NS MASS T
!   9   CSID 1  NS MASS T2      MATID2  I       MATID2  CSID 1  NS MASS
!  10   X1      CSID 1  NS MASS T2      MATID3  I       X1      CSID 1
!  11   Y1      X1      Z1      NS MASS T2      MATID3  Y1      X1
!  12   Z1      Y1      Z2      Z1      NS MASS T2      Z1      Y1
!  13   CSID 2  Z1      CSID 1  Z2      Z1      NS MASS CSID 2  Z1
!  14   X2      CSID 2  X1      CSID 1  Z2      Z1      X2      CSID 2
!  15   Y2      X2      Y1      X1      CSID 1  Z2      Y2      X2
!  16   Z2      Y2      Z1      Y1      X1      CSID 1  Z2      Y2
!  17   CSID 3  Z2      CSID 2  Z1      Y1      X1      CSID 3  Z2
!  18   X3      CSID 3  X2      CSID 2  Z1      Y1      X3      CSID 3
!  19   Y3      X3      Y2      X2      CSID 2  Z1      Y3      X3
!  20   Z3      Y3      Z2      Y2      X2      CSID 2  Z3      Y3
!  21   TEMP    Z3      CSID 3  Z2      Y2      X2      TEMP    Z3
!  22           CSID 4  X3      CSID 3  Z2      Y2              CSID 4
!  23           X4      Y3      X3      CSID 3  Z2              X4
!  24           Y4      Z3      Y3      X3      CSID 3          Y4
!  25           Z4      TEMP    Z3      Y3      X3              Z4
!  26           TEMP            CSID 4  Z3      Y3              TEMP
!  27                           X4      TEMP    Z3
!  28                           Y4              CSID 4
!  29                           Z4              X4
!  30                           TEMP            Y4
!  31                                           Z4
!  32                                           TEMP
! **********************************************************************
 
 
 INTEGER, INTENT(IN OUT)                  :: ntype
 DIMENSION SAVE(32)
 
!     ********FOLLOWING BLOCK IS SET AT MINIMUM LENGTH REQUIRED BY
!     ********THIS ROUTINE.....
 
 COMMON /sdr2x5/  ecpt(100) ,    ph1out(176)
 
 
!     THIS SUBROUTINE INCORPORATES TRIA1, QUAD1, TRIA2, QUAD2
 
!              NTYPE = 1  IMPLIES STRIA1
!              NTYPE = 2  IMPLIES STRIA2
!              NTYPE = 3  IMPLIES SQUAD1
!              NTYPE = 4  IMPLIES SQUAD2
 
!     SAVE THE INCOMING ECPT
 
 DO  i=1,32
   SAVE(i) = ecpt(i)
 END DO
 
!     TRANSFER TO OPERATIONS DESIRED
 
!              STRIA1    STRIA2    SQUAD1    SQUAD2
 SELECT CASE ( ntype )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 150
   CASE (    4)
     GO TO 230
 END SELECT
 
!     **************
! 100 *** STRIA1 ***
!     **************
 
!     SET UP ECPT FOR CALL TO STRME1(0), FIRST CHECK T1 FOR ZERO.
 20 IF( SAVE(7) == 0.0E0 ) GO TO 50
 DO  i=9,21
   ecpt(i) = SAVE(i + 6)
 END DO
 
 CALL strme1(0)
 
!     MOVE OUTPUT FROM STRME1 DOWN TO NEAR BOTTOM OF PH1OUT
!     WORDS (1 THRU 36) DOWN TO (102 THRU 137)
 
 DO  i=1,36
   ph1out(i+101)  = ph1out(i)
 END DO
 GO TO 60
 
 50 ph1out(102) = ecpt(1)
 ph1out(103) = 0.0E0
 
! 150 SET UP ECPT FOR CALL TO STQPL1(3), FIRST CHECK I AND T2 EQUAL ZERO
 60 IF( SAVE(9) == 0.0E0 ) GO TO 90
 DO  i=1,5
   ecpt(i) = SAVE(i)
 END DO
 DO  i=6,25
   ecpt(i) = SAVE(i + 2)
 END DO
 
 CALL strpl1
 RETURN
 
 90 ph1out(1) = ecpt(1)
 ph1out(2) = 0.0E0
 RETURN
 
!     **************
! 200 *** STRIA2 ***
!     **************
 100 IF( SAVE(7) == 0.0E0 ) GO TO 140
!     SET UP ECPT FOR CALL TO STRME1(0)
 
!      ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
 CALL strme1( 0 )
 
!     MOVE OUTPUT FROM STRME1 DOWN TO NEAR BOTTOM OF PH1OUT
!     WORDS (1 THRU 36) DOWN TO (102 THRU 137)
 
 DO  i=1,36
   ph1out(i+101)  = ph1out(i)
 END DO
 
!     SET UP ECPT FOR CALL TO STQPL1(3)
 
 DO  i=1,6
   ecpt(i) = SAVE(i)
 END DO
 ecpt(7) = SAVE(7) ** 3  / 12.0E0
 ecpt(8) = SAVE(6)
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 ecpt(11) =-SAVE(7)/2.0E0
 ecpt(12) = -ecpt(11)
 DO  i=13,25
   ecpt(i) = SAVE(i - 4)
 END DO
 
 CALL strpl1
 RETURN
 
 140 ph1out(  1) = ecpt(1)
 ph1out(  2) = 0.0E0
 ph1out(102) = ecpt(1)
 ph1out(103) = 0.0E0
 RETURN
 
!     **************
! 300 *** SQUAD1 ***
!     **************
 
 150 IF(SAVE(8) == 0.0E0)GO TO 180
 
!     SET UP ECPT FOR CALL TO SQDME1
 
 ecpt(9) = SAVE(13)
 DO  i=10,26
   ecpt(i) = SAVE(i+6)
 END DO
 
 CALL sqdme1
 
!     MOVE OUTPUT FROM SQDME1 DOWN TO BOTTOM OF PH1OUT
!     WORDS (1 THRU 45) DOWN TO (132 THRU 176)
 
 DO  i=1,45
   ph1out(i+131)   = ph1out(i)
 END DO
 
 GO TO 190
 180 ph1out(132) = ecpt(1)
 ph1out(133) = 0.0E0
 
 190 IF( SAVE(10) == 0.0E0 ) GO TO 220
 
!     SET UP ECPT FOR CALL TO STQPL1(4)
 
 DO  i=1,6
   ecpt(i) = SAVE(i)
 END DO
 DO  i=7,30
   ecpt(i) = SAVE(i + 2)
 END DO
 
 CALL sqdpl1
 RETURN
 
 220 ph1out(1) = ecpt(1)
 ph1out(2) = 0.0E0
 RETURN
 
!     **************
! 400 *** SQUAD2 ***
!     **************
 
 230 IF( SAVE(8) == 0.0E0 ) GO TO 270
 
!     SET UP ECPT FOR CALL TO SQDME1
 
!      ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
 CALL sqdme1
!     MOVE OUTPUT FROM SQDME1 DOWN TO BOTTOM OF PH1OUT
!     WORDS (1 THRU 45) DOWN TO (132 THRU 176)
 
 DO  i=1,45
   ph1out(i+131)   = ph1out(i)
 END DO
 
 
!     SET UP ECPT FOR CALL TO STQPL1(4)
 
 DO  i=1,7
   ecpt(i) = SAVE(i)
 END DO
 ecpt(8) = SAVE(8) **3 / 12.0E0
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 ecpt(11)= SAVE(9)
 ecpt(12) =-SAVE(8)/2.0E0
 ecpt(13) =-ecpt(12)
 DO  i=14,30
   ecpt(i) = SAVE(i - 4)
 END DO
 
 CALL sqdpl1
 
 RETURN
 
 270 ph1out(1) = ecpt(1)
 ph1out(2) = 0.0E0
 ph1out(132) = ecpt(1)
 ph1out(133) = 0.0E0
 RETURN
END SUBROUTINE strqd1
