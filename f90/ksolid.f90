SUBROUTINE ksolid (itype)
     
!     IOPT = 1 IMPLIES WEDGE - 3 TETRAHEDRONS
!     IOPT = 2 IMPLIES HEXA(6-SIDED-SOLID) 5  TETRAHEDRONS
!     IOPT = 3 IMPLIES HEXA(6-SIDED-SOLID) 10 TETRAHEDRONS
 
!     ECPT        TETRA          WEDGE          HEXA
!     ------------------------------------------------
!     ECPT( 1) =  EL ID          EL ID          EL ID
!     ECPT( 2) =  MAT-ID         MAT-ID         MAT-ID
!     ECPT( 3) =  GRID-1         GRID-1         GRID-1
!     ECPT( 4) =  GRID-2         GRID-2         GRID-2
!     ECPT( 5) =  GRID-3         GRID-3         GRID-3
!     ECPT( 6) =  GRID-4         GRID-4         GRID-4
!     ECPT( 7) =  CSID-1         GRID-5         GRID-5
!     ECPT( 8) =  X1             GRID-6         GRID-6
!     ECPT( 9) =  Y1             CSID-1         GRID-7
!     ECPT(10) =  Z1             X1             GRID-8
!     ECPT(11) =  CSID-2         Y1             CSID-1
!     ECPT(12) =  X2             Z1             X1
!     ECPT(13) =  Y2             CSID-2         Y1
!     ECPT(14) =  Z2             X2             Z1
!     ECPT(15) =  CSID-3         Y2             CSID-2
!     ECPT(16) =  X3             Z2             X2
!     ECPT(17) =  Y3             CSID-3         Y2
!     ECPT(18) =  Z3             X3             Z2
!     ECPT(19) =  CSID-4         Y3             CSID-3
!     ECPT(20) =  X4             Z3             X3
!     ECPT(21) =  Y4             CSID-4         Y3
!     ECPT(22) =  Z4             X4             Z3
!     ECPT(23) =  EL-TEM         Y4             CSID-4
!     ECPT(24)                   Z4             X4
!     ECPT(25)                   CSID-5         Y4
!     ECPT(26)                   X5             Z4
!     ECPT(27)                   Y5             CSID-5
!     ECPT(28)                   Z5             X5
!     ECPT(29)                   CSID-6         Y5
!     ECPT(30)                   X6             Z5
!     ECPT(31)                   Y6             CSID-6
!     ECPT(32)                   Z6             X6
!     ECPT(33)                   ELTEMP         Y6
!     ECPT(34)                                  Z6
!     ECPT(35)                                  CSID-7
!     ECPT(36)                                  X7
!     ECPT(37)                                  Y7
!     ECPT(38)
!     ECPT(39)                                  CSID-8
!     ECPT(40)                                  X8
!     ECPT(41)                                  Y8
!     ECPT(42)                                  Z8
!     ECPT(43)                                  EL-TEMP
 
!     MAP FOR WEDGE  M(I,J)  I = TETRAHEDRON, J = GRID POINT
 
 
 INTEGER, INTENT(IN)                      :: itype
 LOGICAL :: nogo
 INTEGER :: necpt(52),out,m(22,4)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,out,nogo
 COMMON /sma1et/ ecpt(100)
 COMMON /sma1cl/ iopt4,k4ggsw,npvt,iskp(19),nogoo
 COMMON /sma1dp/ r12(3),r13(3),r(3),rxr(3),r24(3)
 EQUIVALENCE     (necpt(1),ecpt(1))
 DATA    m( 1,1),m( 1,2),m( 1,3),m( 1,4) / 1   ,2   ,3   ,4 /
 DATA    m( 2,1),m( 2,2),m( 2,3),m( 2,4) / 1   ,2   ,3   ,5 /
 DATA    m( 3,1),m( 3,2),m( 3,3),m( 3,4) / 1   ,2   ,3   ,6 /
 DATA    m( 4,1),m( 4,2),m( 4,3),m( 4,4) / 1   ,4   ,5   ,6 /
 DATA    m( 5,1),m( 5,2),m( 5,3),m( 5,4) / 2   ,4   ,5   ,6 /
 DATA    m( 6,1),m( 6,2),m( 6,3),m( 6,4) / 3   ,4   ,5   ,6 /
 DATA    m( 7,1),m( 7,2),m( 7,3),m( 7,4) / 2   ,1   ,4   ,6 /
 DATA    m( 8,1),m( 8,2),m( 8,3),m( 8,4) / 2   ,3   ,4   ,6 /
 DATA    m( 9,1),m( 9,2),m( 9,3),m( 9,4) / 1   ,3   ,4   ,5 /
 DATA    m(10,1),m(10,2),m(10,3),m(10,4) / 2   ,3   ,4   ,5 /
 DATA    m(11,1),m(11,2),m(11,3),m(11,4) / 3   ,1   ,5   ,6 /
 DATA    m(12,1),m(12,2),m(12,3),m(12,4) / 2   ,1   ,5   ,6 /
 
!     MAP FOR HEXA-SOLID (5 OR 10 TETRAHEDRONS)
 
 DATA    m(13,1),m(13,2),m(13,3),m(13,4) / 1   ,2   ,3   ,6 /
 DATA    m(14,1),m(14,2),m(14,3),m(14,4) / 1   ,3   ,4   ,8 /
 DATA    m(15,1),m(15,2),m(15,3),m(15,4) / 1   ,3   ,8   ,6 /
 DATA    m(16,1),m(16,2),m(16,3),m(16,4) / 1   ,5   ,6   ,8 /
 DATA    m(17,1),m(17,2),m(17,3),m(17,4) / 3   ,6   ,7   ,8 /
 DATA    m(18,1),m(18,2),m(18,3),m(18,4) / 2   ,3   ,4   ,7 /
 DATA    m(19,1),m(19,2),m(19,3),m(19,4) / 1   ,2   ,4   ,5 /
 DATA    m(20,1),m(20,2),m(20,3),m(20,4) / 2   ,4   ,5   ,7 /
 DATA    m(21,1),m(21,2),m(21,3),m(21,4) / 2   ,5   ,6   ,7 /
 DATA    m(22,1),m(22,2),m(22,3),m(22,4) / 4   ,5   ,7   ,8 /
 DATA    idelem / 0 /
 
!     BRANCH ON ELEMENT TYPE
 
 igflag = 0
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
 
!     COME HERE FOR WEDGE COMPUTATIONS.
!     KTETRA IS CALLED 3 TIMES BASED ON WEDGE MAPPING MATRIX.
 
 1000 itet   = 1
 ntet   = 12
 itemp  = 33
 ngrids = 6
 iopt   = 0
 
!     BASE CROSS PRODUCT
 
 IF (necpt(1) == idelem) GO TO 1951
 idelem = necpt(1)
 igflag = 1
 r12(1) = ecpt(14) - ecpt(10)
 r12(2) = ecpt(15) - ecpt(11)
 r12(3) = ecpt(16) - ecpt(12)
 r13(1) = ecpt(18) - ecpt(10)
 r13(2) = ecpt(19) - ecpt(11)
 r13(3) = ecpt(20) - ecpt(12)
 CALL saxb (r12,r13,rxr)
 
!     IN THE ABOVE, THE WEDGE IS NUMBERED 1,2,3 COUNTERCLOCKWISE AT THE
!     BASE AND 4,5,6 COUNTER CLOCKWISE AT THE TOP. (LOOKING DOWN ON WED)
 
 r12(1) = ecpt(26) - ecpt(22)
 r12(2) = ecpt(27) - ecpt(23)
 r12(3) = ecpt(28) - ecpt(24)
 r13(1) = ecpt(30) - ecpt(22)
 r13(2) = ecpt(31) - ecpt(23)
 r13(3) = ecpt(32) - ecpt(24)
 CALL saxb (r12,r13,r)
 
 IF (sadotb(r,rxr) > 0.0) THEN
   GO TO  1950
 END IF
 
!     ERROR CONDITION - BAD GEOMETRY
 
 1800 WRITE  (out,1900) ufm,necpt(1)
 1900 FORMAT (a23,' 4001, ELEMENT',i10,' HAS BAD GEOMETRY.')
 nogoo = 1
 RETURN
 
!     PLANER CHECKS FOR WEDGE
 
 1950 CALL kpltst (ecpt(10),ecpt(14),ecpt(26),ecpt(22))
 CALL kpltst (ecpt(10),ecpt(22),ecpt(30),ecpt(18))
 CALL kpltst (ecpt(14),ecpt(18),ecpt(30),ecpt(26))
 1951 IF (nogoo == 1) RETURN
 GO TO 6000
 
!     COME HERE FOR 5-TETRAHEDRON 6-SIDED SOLID
 
 2000 itet   = 13
 ntet   = 17
 itemp  = 43
 ngrids = 8
 iopt   = 0
 GO TO 3500
 
!     COME HERE FOR 10-TETRAHEDRON 6-SIDED SOLID
 
 3000 itet   = 13
 ntet   = 22
 itemp  = 43
 ngrids = 8
 iopt   = 1
 
!     CHECK GEOMETRY OF 6-SIDED SOLID AT THIS POINT
 
 3500 IF (necpt(1) == idelem) GO TO 2951
 idelem = necpt(1)
 igflag = 1
 r13(1) = ecpt(20) - ecpt(12)
 r13(2) = ecpt(21) - ecpt(13)
 r13(3) = ecpt(22) - ecpt(14)
 r24(1) = ecpt(24) - ecpt(16)
 r24(2) = ecpt(25) - ecpt(17)
 r24(3) = ecpt(26) - ecpt(18)
 CALL saxb (r13,r24,rxr)
 
 r12(1) = ecpt(36) - ecpt(28)
 r12(2) = ecpt(37) - ecpt(29)
 r12(3) = ecpt(38) - ecpt(30)
 r13(1) = ecpt(40) - ecpt(32)
 r13(2) = ecpt(41) - ecpt(33)
 r13(3) = ecpt(42) - ecpt(34)
 CALL saxb (r12,r13,r)
 
 IF (sadotb(rxr,r) > 0.0) THEN
   GO TO  2950
 ELSE
   GO TO  1800
 END IF
 
!     PLANER CHECKS FOR HEXA-5 OR HEXA-10
 
 2950 CALL kpltst (ecpt(12),ecpt(16),ecpt(20),ecpt(24))
 CALL kpltst (ecpt(12),ecpt(16),ecpt(32),ecpt(28))
 CALL kpltst (ecpt(16),ecpt(20),ecpt(36),ecpt(32))
 CALL kpltst (ecpt(20),ecpt(24),ecpt(40),ecpt(36))
 CALL kpltst (ecpt(24),ecpt(12),ecpt(28),ecpt(40))
 CALL kpltst (ecpt(28),ecpt(32),ecpt(36),ecpt(40))
 2951 IF (nogoo == 1) RETURN
 GO TO 6000
 
!     AT THIS POINT ALL CHECKS HAVE BEEN MADE. NOW FORM THE ECPT FOR
!     EACH TETRAHEDRON AND CALL KTETRA(IOPT). IOPT = 1 IMPLIES TO COMPUT
!     HALF STIFFNESS. IOPT = 0 IMPLIES COMPUTE FULL STIFFNESS.
 
 6000 DO  j = 1,50
   ecpt(j+50) = ecpt(j)
 END DO
 
!     FILL MAT ID AND EL TEMP
 
 necpt( 2) = necpt(52)
 necpt(23) = necpt(itemp+50)
 jtype     = itype
 DO  i = itet,ntet
   IF (i == ntet) jtype = -itype
   IF (itype == 1) iopt = i + 10
   
!     FILL IN GRID SIL-S AND COORDINATE SETS
   
   DO  j = 1,4
     kpoint = m(i,j)
     necpt(j+2) = necpt(kpoint+52)
     kpoint = 4*kpoint + ngrids - 3
     jpoint = 4*j + 2
     necpt(jpoint+1) = necpt(kpoint+52)
     necpt(jpoint+2) = necpt(kpoint+53)
     necpt(jpoint+3) = necpt(kpoint+54)
     necpt(jpoint+4) = necpt(kpoint+55)
   END DO
   
!     BUMP IOPT IF GEOMETRY TESTS ARE TO BE MADE
   
   IF (igflag == 1) iopt = iopt + 100
   CALL ktetra (iopt,jtype)
 END DO
 
!     ALL THROUGH
 
 RETURN
END SUBROUTINE ksolid
