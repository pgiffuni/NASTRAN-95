SUBROUTINE pktrqd (ntype)
!  THIS ROUTINE SETS UP THE ECPT FOR COMBINATION ELEMENTS IN PLA4
 
 
! **********************************************************************
 
!     8/18/67         E C P T     L I S T I N G
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
 COMMON /pla4es/ ecpt(100)
 COMMON /pla42d/ dummy(600)
 EQUIVALENCE (SAVE(1),ecpt(50))
 
!     THIS SUBROUTINE INCORPORATES TRIA1, QUAD1, TRIA2, QUAD2
 
!              NTYPE = 1  IMPLIES KTRIA1
!              NTYPE = 2  IMPLIES KTRIA2
!              NTYPE = 3  IMPLIES KQUAD1
!              NTYPE = 4  IMPLIES KQUAD2
 
 
!     THE SAVED ECPT IS EQUIVALENCED TO ECPT(50)
 
!     SAVE THE INCOMING ECPT
 
 DO  i=1,32
   SAVE(i) = ecpt(i)
 END DO
 
!     TRANSFER TO OPERATIONS DESIRED
 
!              KTRIA1    KTRIA2    KQUAD1    KQUAD2
 SELECT CASE ( ntype )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 70
   CASE (    3)
     GO TO 100
   CASE (    4)
     GO TO 150
 END SELECT
 
!     *** KTRIA1 ***
 
!     SET UP ECPT FOR CALL TO PKTRMS
 20 IF( SAVE(7) == 0.0E0 ) GO TO 40
 DO  i=9,21
   ecpt(i) = SAVE(i + 6)
 END DO
 
 CALL pktrms(0)
 
!     SET UP CALL TO PKTRPL
 40 IF( SAVE(9) == 0.0E0 ) RETURN
 DO  i=1,5
   ecpt(i) = SAVE(i)
 END DO
 DO  i=6,25
   ecpt(i) = SAVE(i + 2)
 END DO
 
 CALL pktrpl
 RETURN
 
!     *** KTRIA2 ***
 
 70 IF( SAVE(7) == 0.0E0 ) RETURN
!     SET UP CALL TO PKTRMS
 
!      ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
 CALL pktrms(0)
 
!     SET UP CALL TO PKTRPL
 
 DO  i=1,6
   ecpt(i) = SAVE(i)
 END DO
 ecpt(7) = SAVE(7) ** 3  / 12.0E0
 ecpt(8) = SAVE(6)
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 DO  i=13,25
   ecpt(i) = SAVE(i - 4)
 END DO
 
 CALL pktrpl
 RETURN
 
!     *** KQUAD1 ***
 
 100 IF(SAVE(8) == 0.0E0)GO TO 120
 
!     SET UP CALL TO PKQDMS
 
 ecpt(9) = SAVE(13)
 DO  i=10,26
   ecpt(i) = SAVE(i+6)
 END DO
 
 CALL pkqdms
 
 120 IF( SAVE(10) == 0.0E0 ) RETURN
 
!     SET UP CALL TO PKQDPL
 
 DO  i=1,6
   ecpt(i) = SAVE(i)
 END DO
 DO  i=7,30
   ecpt(i) = SAVE(i + 2)
 END DO
 
 CALL pkqdpl
 RETURN
 
!     *** KQUAD2 ***
 
 150 IF( SAVE(8) == 0.0E0 ) RETURN
 
!     SET UP CALL TO PKQDMS
 
!      ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
 CALL pkqdms
 
!     SET UP CALL TO PKQDPL
 
 DO  i=1,7
   ecpt(i) = SAVE(i)
 END DO
 ecpt(8) = SAVE(8) **3 / 12.0E0
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 ecpt(11)= SAVE(9)
 DO  i=14,30
   ecpt(i) = SAVE(i - 4)
 END DO
 
 CALL pkqdpl
 
 RETURN
END SUBROUTINE pktrqd
