SUBROUTINE bisloc (*,id,arr,LEN,kn,jloc)
!-----
!     BINARY SEARCH - LOCATE KEY WORD 'ID' IN ARRAY 'ARR', 1ST ENTRY
!     IF FOUND, 'JLOC' IS THE MATCHED POSITION IN 'ARR'
!     IF NOT FOUND, NON-STANDARD RETURN
!                                                     I.E.
!     ID  = KEY WORD TO MATCH IN ARR.      MATCH AGAINST 1ST COL OF ARR
!     ARR = ARRAY TO SEARCH.                          ARR(ROW,COL)
!     LEN = LENGTH OF EACH ENTRY IN ARRAY.            LEN=ROW
!     KN  = NUMBER OF ENTRIES IN THE ARR.             KN =COL
!     JLOC= POINTER RETURNED - FIRST WORD OF ENTRY.   MATCHED ROW
!-----
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: id
 INTEGER, INTENT(IN)                      :: arr(1)
 INTEGER, INTENT(IN)                      :: LEN
 INTEGER, INTENT(IN)                      :: kn
 INTEGER, INTENT(OUT)                     :: jloc
 
 DATA     iswtch / 16 /
 
 jj = LEN - 1
 IF (kn < iswtch) GO TO 120
 klo = 1
 khi = kn
 10 k   = (klo+khi+1)/2
 20 j   = k*LEN - jj
 IF (id-arr(j) < 0) THEN
   GO TO    30
 ELSE IF (id-arr(j) == 0) THEN
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
 60 IF (k == klo) GO TO 70
 k   = klo
 GO TO 80
 70 k   = khi
 80 klo = khi
 GO TO 20
 90 jloc = j
 RETURN
 100 jloc = khi*LEN - jj
 j    = kn *LEN - jj
 IF (id > arr(j)) jloc = jloc + LEN
 110 RETURN 1
 
!     SEQUENTIAL SEARCH MORE EFFICIENT
 
 120 khi = kn*LEN - jj
 DO  j = 1,khi,LEN
   IF (arr(j)-id < 0.0) THEN
     GO TO   130
   ELSE IF (arr(j)-id == 0.0) THEN
     GO TO    90
   ELSE
     GO TO   140
   END IF
 END DO
 jloc = khi + LEN
 GO TO 110
 140 jloc = j
 GO TO 110
END SUBROUTINE bisloc
