SUBROUTINE mrge (list, n, string, m)
!*****
! MRGE IS A MERGE ROUTINE. GIVEN A SORTED LIST AND A SORTED STRING,
! MRGE ADDS THE ENTRIES IN THE STRING TO THE LIST IN THEIR APPROPRIATE
! POSITIONS.  DUPLICATES ARE DISCARDED.
 
!  ARGUMENTS
 
!     LIST   --- THE ARRAY CONTAINING THE SORTED LIST
!     N      --- THE NUMBER OF TERMS BEFORE AND AFTER THE MERGE
!     STRING --- THE ARRAY CONTAINING THE SORTED STRING
!     M      --- THE NUMBER OF TERMS IN THE STRING
 
!*****
 
 
 INTEGER, INTENT(OUT)                     :: list(1)
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN)                      :: string(1)
 INTEGER, INTENT(IN)                      :: m
 
 
! LOCATE THE POSITION IN THE LIST OF THE FIRST TERM IN THE STRING
 
 kk = 1
 id = string(kk)
 CALL bisloc (*12, id, list, 1, n, k)
 kstart = MIN0( k+1, n )
 k2 = 2
 GO TO 13
 12 kstart = MAX0( 1, k-1 )
 k2 = 1
 
! CREATE A HOLE IN THE LIST BY MOVING THE END OF THE LIST.
 
 13 k = n
 14 list(k+m) = list(k)
 k = k - 1
 IF( k >= kstart ) GO TO 14
 k1 = kstart + m
 nm = n + m
 k = kstart
 
! NOW ADD TO THE LIST BY MERGING FROM THE TWO STRINGS
 
 16 IF( k1 > nm ) GO TO 60
 IF( k2 > m  ) GO TO 50
 IF (list(k1) - string(k2) < 0) THEN
   GO TO    20
 ELSE IF (list(k1) - string(k2) == 0) THEN
   GO TO    40
 ELSE
   GO TO    30
 END IF
 
!    CHOOSE TERM FROM OLD LIST
 
 20 list(k) = list(k1)
 k1 = k1 + 1
 k  = k  + 1
 GO TO 16
 
!    CHOOSE TERM FROM STRING
 
 30 list(k) = string(k2)
 k2 = k2 + 1
 k  = k  + 1
 GO TO 16
 
!    DUPLICATES -- DISCARD TERM FROM STRING
 
 40 k2 = k2 + 1
 GO TO 20
 
!    STRING EXHAUSTED -- COMPLETE LIST FROM OLD LIST
 
 50 DO  kx=k1,nm
   list(k) = list(kx)
   k = k + 1
 END DO
 GO TO 68
 
!    OLD LIST EXHAUSTED -- COMPLETE LIST FROM STRING
 
 60 IF( k2 > m ) GO TO 68
 DO  kx=k2,m
   list(k) = string(kx)
   k = k + 1
 END DO
 
! RETURN NEW NUMBER OF TERMS IN LIST.
 
 68 n = k - 1
 RETURN
END SUBROUTINE mrge
