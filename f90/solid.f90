SUBROUTINE solid( temps, pg, itype )
!*****
!  ELEMENT THERMAL LOAD GENERATOR FOR THE WEDGE, HEXA1, AND HEXA2
 
!  ITYPE = 1 IMPLIES WEDGE - 3 TETRAHEDRONS
 
!  ITYPE = 2 IMPLIES HEXA(6-SIDED-SOLID) 5 TETRAHEDRONS
 
!  ITYPE = 3 IMPLIES HEXA(6-SIDED-SOLID) 10 TETRAHEDRONS
 
!*****
 
 REAL, INTENT(IN)                         :: temps(8)
 REAL, INTENT(IN OUT)                     :: pg(6)
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER :: necpt(52)           ,m(13,4)
 
 REAL :: tmps(4)
 
 COMMON/trimex/ecpt(100)
 
 EQUIVALENCE( necpt(1), ecpt(1) )
!*****
 
!  E C P T     TETRA          WEDGE          HEXA
!  -----------------------------------------------
!  ECPT( 1) =  EL ID          EL ID          EL ID
!  ECPT( 2) =  MAT-ID         MAT-ID         MAT-ID
!  ECPT( 3) =  GRID-1         GRID-1         GRID-1
!  ECPT( 4) =  GRID-2         GRID-2         GRID-2
!  ECPT( 5) =  GRID-3         GRID-3         GRID-3
!  ECPT( 6) =  GRID-4         GRID-4         GRID-4
!  ECPT( 7) =  CSID-1         GRID-5         GRID-5
!  ECPT( 8) =  X1             GRID-6         GRID-6
!  ECPT( 9) =  Y1             CSID-1         GRID-7
!  ECPT(10) =  Z1             X1             GRID-8
!  ECPT(11) =  CSID-2         Y1             CSID-1
!  ECPT(12) =  X2             Z1             X1
!  ECPT(13) =  Y2             CSID-2         Y1
!  ECPT(14) =  Z2             X2             Z1
!  ECPT(15) =  CSID-3         Y2             CSID-2
!  ECPT(16) =  X3             Z2             X2
!  ECPT(17) =  Y3             CSID-3         Y2
!  ECPT(18) =  Z3             X3             Z2
!  ECPT(19) =  CSID-4         Y3             CSID-3
!  ECPT(20) =  X4             Z3             X3
!  ECPT(21) =  Y4             CSID-4         Y3
!  ECPT(22) =  Z4             X4             Z3
!  ECPT(23) =  EL-TEM         Y4             CSID-4
!  ECPT(24)                   Z4             X4
!  ECPT(25)                   CSID-5         Y4
!  ECPT(26)                   X5             Z4
!  ECPT(27)                   Y5             CSID-5
!  ECPT(28)                   Z5             X5
!  ECPT(29)                   CSID-6         Y5
!  ECPT(30)                   X6             Z5
!  ECPT(31)                   Y6             CSID-6
!  ECPT(32)                   Z6             X6
!  ECPT(33)                   ELTEMP         Y6
!  ECPT(34)                                  Z6
!  ECPT(35)                                  CSID-7
!  ECPT(36)                                  X7
!  ECPT(37)                                  Y7
!  ECPT(38)
!  ECPT(39)                                  CSID-8
!  ECPT(40)                                  X8
!  ECPT(41)                                  Y8
!  ECPT(42)                                  Z8
!  ECPT(43)                                  EL-TEMP
!*****
 
!*****
!  MAP FOR WEDGE  M(I,J)  I=TETRAHEDRON, J=GRID POINT
!*****
 DATA m( 1,1),m( 1,2),m( 1,3),m( 1,4) / 1   ,2   ,3   ,6 /
 DATA m( 2,1),m( 2,2),m( 2,3),m( 2,4) / 1   ,2   ,6   ,5 /
 DATA m( 3,1),m( 3,2),m( 3,3),m( 3,4) / 1   ,4   ,5   ,6 /
!*****
!  MAP FOR HEXA-SOLID (5 OR 10 TETRAHEDRONS)
!*****
 DATA m( 4,1),m( 4,2),m( 4,3),m( 4,4) / 1   ,2   ,3   ,6 /
 DATA m( 5,1),m( 5,2),m( 5,3),m( 5,4) / 1   ,3   ,4   ,8 /
 DATA m( 6,1),m( 6,2),m( 6,3),m( 6,4) / 1   ,3   ,8   ,6 /
 DATA m( 7,1),m( 7,2),m( 7,3),m( 7,4) / 1   ,5   ,6   ,8 /
 DATA m( 8,1),m( 8,2),m( 8,3),m( 8,4) / 3   ,6   ,7   ,8 /
 DATA m( 9,1),m( 9,2),m( 9,3),m( 9,4) / 2   ,3   ,4   ,7 /
 DATA m(10,1),m(10,2),m(10,3),m(10,4) / 1   ,2   ,4   ,5 /
 DATA m(11,1),m(11,2),m(11,3),m(11,4) / 2   ,4   ,5   ,7 /
 DATA m(12,1),m(12,2),m(12,3),m(12,4) / 2   ,5   ,6   ,7 /
 DATA m(13,1),m(13,2),m(13,3),m(13,4) / 4   ,5   ,7   ,8 /
!*****
!  BRANCH ON ELEMENT TYPE
!*****
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
!*****
!  COME HERE FOR WEDGE COMPUTATIONS.
!  KTETRA IS CALLED 3 TIMES BASED ON WEDGE MAPPING MATRIX.
!*****
 1000 itet = 1
 ntet = 3
 itemp= 33
 ngrids = 6
 iopt = 0
 GO TO 6000
!*****
!  COME HERE FOR 5-TETRAHEDRON 6-SIDED SOLID
!*****
 2000 itet = 4
 ntet = 8
 itemp= 43
 ngrids = 8
 iopt = 0
 GO TO 6000
!*****
!  COME HERE FOR 10-TETRAHEDRON 6-SIDED SOLID
!*****
 3000 itet = 4
 ntet =13
 itemp= 43
 ngrids = 8
 iopt = 1
 GO TO 6000
 6000 DO  j = 1,50
   ecpt(j+50) = ecpt(j)
 END DO
 
!  FILL MAT ID AND EL TEMP
 
 necpt(2) = necpt(52)
 necpt(23) = necpt (itemp+50)
 DO  i = itet,ntet
   
!     FILL IN GRID SIL-S AND COORDINATE SETS
   
   DO  j = 1,4
     kpoint = m(i,j)
     tmps(j) = temps(kpoint)
     necpt(j+2) = necpt(kpoint+52)
     kpoint = 4*kpoint + ngrids - 3
     jpoint = 4*j + 2
     necpt(jpoint+1) = necpt(kpoint+52)
     necpt(jpoint+2) = necpt(kpoint+53)
     necpt(jpoint+3) = necpt(kpoint+54)
     necpt(jpoint+4) = necpt(kpoint+55)
   END DO
   CALL tetra( tmps(1), pg(1), iopt )
 END DO
!*****
!  ALL THROUGH
!*****
 RETURN
END SUBROUTINE solid
