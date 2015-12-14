SUBROUTINE fvrs1e (a,k,n)
     
!     PURPOSE
!       TO SORT THE ELEMENTS OF A REAL*4 VECTOR, A, INTO ASCENDING
!       ORDER AND TO CONSTRUCT AN INTEGER*4 VECTOR, K, WHICH INDICATES
!       HOW THE ELEMENTS OF A HAVE BEEN REARRANGED.
 
!     USAGE
!       CALL FVRS1E(A,K,N)
 
!     DESCRIPTION OF PARAMETERS
!       A - REAL*4 VECTOR.
!              ON INPUT - A CONTAINS THE NUMBERS TO BE SORTED.
!              ON OUTPUT- A CONTAINS THE NUMBERS IN ASCENDING ORDER.
!       K - OUTPUT VECTOR CONTAINING INTERCHANGE INFORMATION, I.E.,
!           THE NUMBER IN A(K(I)) (OF THE INPUT A) HAS BEEN MOVED TO
!           A(I).
!       N - LENGTH OF A AND K.
 
!     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!       NONE
 
!     REMARKS
!       THE K-VECTOR CAN BE USED IN CONJUNCTION WITH SUBROUTINE FVRS1E
!       TO REARRANGE OTHER VECTORS IN THE SAME WAY THAT THE A-VECTOR
!       HAS BEEN REARRANGED.
 
!     METHOD
!       THIS ROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE,
!       'SHELLSORT', ALGORITHM 201, 'COLLECTED ALGORITHMS FROM CACM',
!       BY J. BOOTHROYD.
 
 
 
 REAL, INTENT(IN OUT)                     :: a(1)
 INTEGER, INTENT(OUT)                     :: k(1)
 INTEGER, INTENT(IN)                      :: n
 
 
 DO  ikl =1,n
   k(ikl) = ikl
 END DO
 i = 1
 1 i = i+i
 IF(i-n < 0) THEN
   GO TO     1
 ELSE IF (i-n == 0) THEN
   GO TO     2
 END IF
 7 i = i/2
 2 CONTINUE
 m = 2*i-1
 5 CONTINUE
 m = m/2
 k1 = n-m
 DO  j=1,k1
   i = j
   3 ipm = i+m
   aipm = a(ipm)
   IF(aipm >= a(i)) GO TO 4
   w = a(i)
   kw = k(i)
   a(i) = aipm
   k(i) = k(ipm)
   a(ipm) = w
   k(ipm) = kw
   i = i-m
   IF(i >= 1) GO TO 3
   4 CONTINUE
 END DO
 IF(m > 1) GO TO 5
 RETURN
END SUBROUTINE fvrs1e
