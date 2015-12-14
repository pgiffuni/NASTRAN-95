SUBROUTINE bislc2 (*,id,aa,nc,nr,loc)
!-----
!     BINARY SEARCH ROUTINE - LOCATE ID POSTION IN AA
!     SEARCH BY FIRST 2 WORDS (ROWS) OF ENTRIES.
 
!     ID  = TARGET WORD SEARCH, 2 BCD-WORDS
!     AA  = A (NR X NC) TABLE TO SEARCH FOR ID.
!     NR  = SIZE   OF ENTRIES (ROW   ) IN THE AA.
!     NC  = NUMBER OF ENTRIES (COLUMN) IN THE AA.
!     LOC = POINTER RETURNED, OF NC LOCATION
 
!     NONSTANDARD RETURN IN THE EVENT OF NO MATCH.
 
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: id(2)
 INTEGER, INTENT(IN OUT)                  :: aa(nr,nc)
 INTEGER, INTENT(IN)                      :: nc
 INTEGER, INTENT(IN OUT)                  :: nr
 INTEGER, INTENT(OUT)                     :: loc
 
 
 klo = 1
 khi = nc
 10 k   = (klo+khi+1)/2
 20 IF (id(1) - aa(1,k) < 0) THEN
   GO TO    30
 ELSE IF (id(1) - aa(1,k) == 0) THEN
   GO TO    25
 ELSE
   GO TO    40
 END IF
 25 IF (id(2) - aa(2,k) < 0) THEN
   GO TO    30
 ELSE IF (id(2) - aa(2,k) == 0) THEN
   GO TO    90
 ELSE
   GO TO    40
 END IF
 30 khi = k
 GO TO 50
 40 klo = k
 50 IF (khi-klo -1 < 0) THEN
   GO TO   100
 ELSE IF (khi-klo -1 == 0) THEN
   GO TO    60
 ELSE
   GO TO    10
 END IF
 60 IF (k == klo)  GO TO 70
 k   = klo
 GO TO 80
 70 k   = khi
 80 klo = khi
 GO TO 20
 90 loc = k
 RETURN
 100 RETURN 1
END SUBROUTINE bislc2
