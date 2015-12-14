SUBROUTINE ktriqd (ntype)
     
 
!     8/18/67         E C P T     L I S T I N G
 
! ECPT  TRMEM   QDMEM   TRPLT   QDPLT   TRIA1   QUAD1   TRIA2   QUAD2
! ***** ******* ******* ******* ******* ******* ******* ******* ********
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
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ntype
 LOGICAL :: heat
 INTEGER :: scr4,iecpt(4),bcd(2,4),bgpdt(4)
 DIMENSION       SAVE(32)
 COMMON /BLANK / skip(16),volume,surfac
 COMMON /matout/ dum(6),rho
 COMMON /sma1ht/ heat
 COMMON /sma1et/ ecpt(100)
 COMMON /sma1dp/ dummy(600)
 EQUIVALENCE     (SAVE(1),ecpt(50)),(ecpt(1),iecpt(1))
 DATA    bcd   / 4HCTRI,2HA1,4HCTRI,2HA2,4HCQUA,2HD1,4HCQUA,2HD2 /
 DATA    OLD   , kount,ngpt / 0.0,   2*0     /
 DATA    scr4  , bgpdt/  304, 15, 9, 16, 10  /
 
!     THIS SUBROUTINE INCORPORATES TRIA1, QUAD1, TRIA2, QUAD2
 
!             NTYPE = 1  IMPLIES KTRIA1
!             NTYPE = 2  IMPLIES KTRIA2
!             NTYPE = 3  IMPLIES KQUAD1
!             NTYPE = 4  IMPLIES KQUAD2
 
!     CALLS FROM THIS ROUTINE CAN BE MADE TO
 
!            KTRMEM - TRIANGULAR MEMBRANE ROUTINE
!            KQDMEM - QUADRILATERAL MEMBRANE ROUTINE
!            KTQPLT - TRIANGULAR OR QUADRILATERAL PLATE ROUTINE
!            QDMM1X - HIGHER LEVEL QUADRIALATER MEMBRANE ROUTINE
 
!     ALL INSERTIONS OF 6X6 ELEMENT STIFFNESS MATRICES ARE HANDLED BY
!     THE ABOVE ROUTINES.
 
 
!     THE SAVED ECPT IS EQUIVALENCED TO ECPT(50)
 
!     SAVE THE INCOMING ECPT
 
 DO  i = 1,32
   SAVE(i) = ecpt(i)
 END DO
 
!     TRANSFER TO OPERATIONS DESIRED
 
!           KTRIA1 KTRIA2 KQUAD1 KQUAD2
 SELECT CASE ( ntype )
   CASE (    1)
     GO TO   20
   CASE (    2)
     GO TO     70
   CASE (    3)
     GO TO    100
   CASE (    4)
     GO TO    150
 END SELECT
 
!     *** KTRIA1 ***
 
!     SET UP ECPT FOR CALL TO KTRMEM (0), FIRST CHECK T1 FOR ZERO.
 
 20 IF (SAVE(7) == 0.0) GO TO 40
 DO  i = 9,21
   ecpt(i) = SAVE(i+6)
 END DO
 
 CALL ktrmem (0)
 
!     SET UP ECPT FOR CALL TO TQPLT(3), FIRST CHECK I AND T2 EQUAL ZERO.
 
 40 IF (SAVE(9) == 0.0) GO TO 200
 DO  i = 1,5
   ecpt(i) = SAVE(i)
 END DO
 DO  i = 6,25
   ecpt(i) = SAVE(i+2)
 END DO
 
 IF (.NOT.heat) CALL ktrplt
 GO TO 200
 
!     *** KTRIA2 ***
 
 70 IF (SAVE(7) == 0.0) GO TO 200
 
!     SET UP ECPT FOR CALL TO KTRMEM (0)
 
!     ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
 CALL ktrmem (0)
 
!     SET UP ECPT FOR CALL TO KTQPLT (3)
 
 DO  i = 1,6
   ecpt(i) = SAVE(i)
 END DO
 ecpt(7) = SAVE(7)**3/12.0
 ecpt(8) = SAVE(6)
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 DO  i = 13,25
   ecpt(i) = SAVE(i-4)
 END DO
 
 IF (.NOT.heat) CALL ktrplt
 GO TO 200
 
!     *** KQUAD1 ***
 
 100 IF (SAVE(8) == 0.0) GO TO 120
 
!     SET UP ECPT FOR CALL TO KQDMEM
 
 ecpt(9) = SAVE(13)
 DO  i = 10,26
   ecpt(i) = SAVE(i+6)
 END DO
 
 CALL kqdmem
 
 120 IF (SAVE(10) == 0.0) GO TO 200
 
!     SET UP ECPT FOR CALL TO KTQPLT (4)
 
 DO  i = 1,6
   ecpt(i) = SAVE(i)
 END DO
 DO  i = 7,30
   ecpt(i) = SAVE(i+2)
 END DO
 
 IF (.NOT.heat) CALL kqdplt
 GO TO 200
 
!     *** KQUAD2 ***
 
 150 IF (SAVE(8) == 0.0) GO TO 200
 
!     SET UP ECPT FOR CALL TO KQDMEM
!     (WHICH HAS WEAK STIFFNESS MATRIX FORMULATION)
!     OR
!     SET UP ECPT FOR CALL TO QDMM1D/S (BETTER MEMBRANE FORMALATION)
!     THE PROBLEM HERE IS THAT KTRIQD AND KQDMEM ARE EMGOLD ELEMENTS
!     WHILE QDMM1D/S ARE EMGPRO NEW ELEMENTS.
!     TO SOLVE THIS PROPLEM, WE NEED A S.P./D.P. QDMM1X ELEMENT ROUTINE
!     THAT USES QDMM1D/S FORMULATION WITH EMGOLD/SMA1B TECHNIQUE.
 
!     ECPT IS OK AS DELIVERED TO THIS ROUTINE
 
!     CALL QDMM1X
!     (QDMM1X IS INCOMPLETE AS OF 3/92. GO BACK TO KQDMEM)
 
 CALL kqdmem
 
!     SET UP ECPT FOR CALL TO KTQPLT (4)
 
 DO  i = 1,7
   ecpt(i) = SAVE(i)
 END DO
 ecpt(8) = SAVE(8)**3/12.0
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 ecpt(11)= SAVE(9)
 DO  i = 14,30
   ecpt(i) = SAVE(i-4)
 END DO
 
 IF (.NOT. heat) CALL kqdplt
 
 
!     SAVE ELEMENT NAME, ID, THICKNESS, DENSITY, NO. OF GRID POINTS,
!     AND GRID PT DATA IF USER REQUESTED VOLUME AND AREA COMPUTATION
 
 200 IF (heat .OR. (volume <= 0.0 .AND. surfac <= 0.0)) GO TO 220
 IF (SAVE(1) == OLD) GO TO 210
 OLD  = SAVE(1)
 ngpt = 3
 IF (ntype >= 3) ngpt = 4
 kount = 0
 210 kount = kount + 1
 IF (kount < ngpt) GO TO 220
 ecpt(2) = SAVE(7)
 ecpt(3) = rho
 iecpt(4)= ngpt
 i = bgpdt(ntype)
 k = ngpt*4
 IF (ntype >= 3) ecpt(2) = SAVE(8)
 CALL WRITE (scr4,bcd(1,ntype),2,0)
 CALL WRITE (scr4,ecpt(1),4,0)
 CALL WRITE (scr4,SAVE(2),ngpt,0)
 CALL WRITE (scr4,SAVE(i),k,1)
 220 RETURN
END SUBROUTINE ktriqd
