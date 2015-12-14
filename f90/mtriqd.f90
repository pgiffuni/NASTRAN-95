SUBROUTINE mtriqd (ntype)
     
 
!     8/18/67           E C P T     L I S T I N G
 
!     ECPT    TRPLT     TRIA1     TRIA2     QDPLT     QUAD1     QUAD2
!     *****************************************************************
 
!      1     ELEM ID   ELEM ID   ELEM ID   ELEM ID   ELEM ID   ELEM ID
!      2     GRID A    GRID A    GRID A    GRID A    GRID A    GRID A
!      3     GRID B    GRID B    GRID B    GRID B    GRID B    GRID B
!      4     GRID C    GRID C    GRID C    GRID C    GRID C    GRID C
!      5     THETA     THETA     THETA     GRID D    GRID D    GRID D
!      6     MATID1    MATID1    MAT ID    THETA     THETA     THETA
!      7     I         T1        T         MATID1    MATID1    MAT ID
!      8     MATID2    MATID2    NS MASS   I         T1        T
!      9     T2        I         CSID 1    MATID2    MATID2    NS MASS
!     10     NS MASS   MATID3    X1        T2        I         CSID 1
!     11     Z1        T2        Y1        NS MASS   MATID3    X1
!     12     Z2        NS MASS   Z1        Z1        T2        Y1
!     13     CSID 1    Z1        CSID 3    Z2        NS MASS   Z1
!     14     X1        Z2        X2        CSID 1    Z1        CSID 2
!     15     Y1        CSID 1    Y2        X1        Z2        X2
!     16     Z1        X1        Z2        Y1        CSID 1    Y2
!     17     CSID 2    Y1        CSID 3    Z1        X1        Z2
!     18     X2        Z1        X3        CSID 2    Y1        CSID 3
!     19     Y2        CSID 2    Y3        X2        Z1        X3
!     20     Z2        X2        Z3        Y2        CSID 2    Y3
!     21     CSID 3    Y2        TEMP      Z2        X2        Z3
!     22     X3        Z2                  CSID 3    Y2        CSID 4
!     23     Y3        CSID 3              X3        Z2        X4
!     24     Z3        X3                  Y3        CSID 3    Y4
!     25     TEMP      Y3                  Z3        X3        Z4
!     26               Z3                  CSID 4    Y3        TEMP
!     27               TEMP                X4        Z3
!     28                                   Y4        CSID 4
!     29                                   Z4        X4
!     30                                   TEMP      Y4
!     31                                             Z4
!     32                                             TEMP
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ntype
 LOGICAL :: heat
 DIMENSION       SAVE(32),isave(32)
 COMMON /sma2et/ ecpt(100)
 COMMON /sma2ht/ heat
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ rho
 EQUIVALENCE     (SAVE(1),isave(1),ecpt(50))
 
!     THIS SUBROUTINE INCORPORATES TRIA1, QUAD1, TRIA2, QUAD2
 
!             NTYPE = 1  IMPLIES MTRIA1
!             NTYPE = 2  IMPLIES MTRIA2
!             NTYPE = 3  IMPLIES MQUAD1
!             NTYPE = 4  IMPLIES MQUAD2
 
!     IF (I . EQ. 0) THEN COMPUTE UNCOUPLED MASS
 
!     CALL MASSTQ (NARG)
!          WHERE   NARG = 5 FOR TRIA1
!                  NARG = 4 FOR TRIA2
!                  NARG = 2 FOR QUAD1
!                  NARG = 1 FOR QUAD2
 
!     CALLS FROM THIS ROUTINE CAN BE MADE TO
 
!            MTRPLT - TRIANGULAR PLATE ROUTINE
!            MQDPLT - QUADRILATERAL PLATE ROUTINE
!            MASSTQ - UNCOUPLED MASS COMBINATION ELEMENT ROUTINE
 
!     ALL INSERTIONS OF 6X6 ELEMENT MASS MATRICES ARE HANDLED BY
!     THE ABOVE ROUTINES.
 
!     THE SAVED ECPT IS EQUIVALENCED TO ECPT(50)
 
 
!     SAVE THE INCOMING ECPT
 
 inflag = 4
 DO  i = 1,32
   SAVE(i) = ecpt(i)
 END DO
 
!     TRANSFER TO OPERATIONS DESIRED
 
!            MTRIA1 MTRIA2 MQUAD1 MQUAD2
 SELECT CASE ( ntype )
   CASE (    1)
     GO TO    20
   CASE (    2)
     GO TO     60
   CASE (    3)
     GO TO    100
   CASE (    4)
     GO TO    150
 END SELECT
 
!     *** MTRIA1 ***
 
!     SET UP ECPT FOR CALL TO MTRPLT.  FIRST CHECK I EQUAL ZERO
 
 20 IF (SAVE(9) /= 0.0) GO TO 30
 narg = 5
 CALL masstq (narg)
 GO TO 200
 
 30 DO  i = 1,5
   ecpt(i) = SAVE(i)
 END DO
 DO  i = 6,25
   ecpt(i) = SAVE(i+2)
 END DO
 matid   = isave(6)
 IF (SAVE(7) == 0.0) GO TO 54
 CALL mat (ecpt(1))
 ecpt(10) = SAVE(12) + rho*SAVE(7)
 
 GO TO 56
 54 ecpt(10) = SAVE(12)
 56 IF (.NOT.heat) CALL mtrplt
 GO TO 200
 
!     *** MTRIA2 ***
 
!     SET UP ECPT FOR CALL TO MTRPLT
 
 60 IF (SAVE(7) /= 0.0) GO TO 70
 narg = 4
 CALL masstq (narg)
 GO TO 200
 
 70 DO  i = 1,6
   ecpt(i) = SAVE(i)
 END DO
 ecpt(7) = SAVE(7)**3/12.0
 ecpt(8) = SAVE(6)
 ecpt(9) = SAVE(7)
 matid   = isave(6)
 CALL mat (ecpt(1))
 ecpt(10) = SAVE(8) + rho*SAVE(7)
 DO  i = 13,25
   ecpt(i) = SAVE(i-4)
 END DO
 
 IF (.NOT. heat) CALL mtrplt
 GO TO 200
 
!     *** MQUAD1 ***
 
!      SET UP ECPT FOR CALL TO MQDPLT.  FIRST CHECK I EQUAL ZERO
 
 100 IF (SAVE(10) /= 0.0) GO TO 110
 narg = 2
 CALL masstq (narg)
 GO TO 200
 
 110 DO  i = 1,6
   ecpt(i) = SAVE(i)
 END DO
 DO  i = 7,30
   ecpt(i) = SAVE(i+2)
 END DO
 matid = isave(7)
 IF (SAVE(8) == 0.0) GO TO 144
 CALL mat (ecpt(1))
 ecpt(11) = SAVE(13) + rho*SAVE(8)
 
 GO TO 146
 144 ecpt(11) = SAVE(13)
 146 IF (.NOT.heat) CALL mqdplt
 GO TO 200
 
!     *** MQUAD2 ***
 
!     SET UP ECPT FOR CALL TO MQDPLT
 
 150 IF (SAVE(8) /= 0.0) GO TO 160
 narg = 1
 CALL masstq (narg)
 GO TO 200
 
 160 DO  i = 1,7
   ecpt(i) = SAVE(i)
 END DO
 ecpt(8) = SAVE(8)**3/12.0
 ecpt(9) = SAVE(7)
 ecpt(10)= SAVE(8)
 matid   = isave(7)
 CALL mat (ecpt(1))
 ecpt(11) = SAVE(9) + rho*SAVE(8)
 DO  i = 14,30
   ecpt(i) = SAVE(i-4)
 END DO
 
 IF (.NOT. heat) CALL mqdplt
 200 RETURN
END SUBROUTINE mtriqd
