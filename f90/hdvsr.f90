SUBROUTINE hdvsr(a,la,ir)
!   IMSL ROUTINE NAME   - HDVSR
 
!-----------------------------------------------------------------------
 
!   COMPUTER            - CDC/SINGLE
 
!   LATEST REVISION     - JANUARY 1, 1978
 
!   PURPOSE             - SORTING OF ARRAYS BY ALGEBRAIC VALUE -
!                           PERMUTATIONS RETURNED
 
!   USAGE               - CALL HDVSR (A,LA,IR)
 
!   ARGUMENTS    A      - ON INPUT, A CONTAINS THE ARRAY TO BE SORTED.
!                         ON OUTPUT, A CONTAINS THE SORTED ARRAY.
!                LA     - INPUT VARIABLE CONTAINING THE NUMBER OF
!                           ELEMENTS IN THE ARRAY TO BE SORTED.
!                IR     - VECTOR OF LENGTH LA.
!                         ON INPUT, IR CONTAINS THE INTEGER VALUES
!                           1,2,...,LA. SEE REMARKS.
!                         ON OUTPUT, IR CONTAINS A RECORD OF THE
!                           PERMUTATIONS MADE ON THE VECTOR A.
 
!   PRECISION/HARDWARE  - SINGLE/ALL
 
!   REQD. IMSL ROUTINES - NONE REQUIRED
 
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
 
!   REMARKS      THE VECTOR IR MUST BE INITIALIZED BEFORE ENTERING
!                HDVSR.  ORDINARILY, IR(1)=1, IR(2)=2, ...,
!                IR(LA)=LA.  FOR WIDER APPLICABILITY, ANY INTEGER
!                THAT IS TO BE ASSOCIATED WITH A(I) FOR I=1,2,...,LA
!                MAY BE ENTERED INTO IR(I).
 
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
 
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE.  NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
 
!-----------------------------------------------------------------------
 
 
 REAL, INTENT(IN OUT)                     :: a(1)
 INTEGER, INTENT(IN)                      :: la
 INTEGER, INTENT(IN OUT)                  :: ir(1)
 
!                                  SPECIFICATIONS FOR ARGUMENTS
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
 INTEGER :: iu(21),il(21),i,m,j,k,ij,it,l,itt
 REAL :: t,tt,r
!                                  FIRST EXECUTABLE STATEMENT
 IF (la <= 0) RETURN
 m = 1
 i = 1
 j = la
 r = .375
 5 IF (i == j) GO TO 45
 IF (r > .5898437) GO TO 10
 r = r+3.90625E-2
 GO TO 15
 10 r = r-.21875
 15 k = i
!                                  SELECT A CENTRAL ELEMENT OF THE
!                                  ARRAY AND SAVE IT IN LOCATION T
 ij = i+(j-i)*r
 t = a(ij)
 it = ir(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
 IF (a(i) <= t) GO TO 20
 a(ij) = a(i)
 a(i) = t
 t = a(ij)
 ir(ij) = ir(i)
 ir(i) = it
 it = ir(ij)
 20 l = j
!                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
!                                  T, INTERCHANGE WITH T
 IF (a(j) >= t) GO TO 30
 a(ij) = a(j)
 a(j) = t
 t = a(ij)
 ir(ij) = ir(j)
 ir(j) = it
 it = ir(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
 IF (a(i) <= t) GO TO 30
 a(ij) = a(i)
 a(i) = t
 t = a(ij)
 ir(ij) = ir(i)
 ir(i) = it
 it = ir(ij)
 GO TO 30
 25 IF (a(l) == a(k)) GO TO 30
 tt = a(l)
 a(l) = a(k)
 a(k) = tt
 itt = ir(l)
 ir(l) = ir(k)
 ir(k) = itt
!                                  FIND AN ELEMENT IN THE SECOND HALF OF
!                                  THE ARRAY WHICH IS SMALLER THAN T
 30 l = l-1
 IF (a(l) > t) GO TO 30
!                                  FIND AN ELEMENT IN THE FIRST HALF OF
!                                  THE ARRAY WHICH IS GREATER THAN T
 35 k = k+1
 IF (a(k) < t) GO TO 35
!                                  INTERCHANGE THESE ELEMENTS
 IF (k <= l) GO TO 25
!                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
!                                  THE ARRAY YET TO BE SORTED
!                                  THE ARRAY YET TO BE SORTED
 IF (l-i <= j-k) GO TO 40
 il(m) = i
 iu(m) = l
 i = k
 m = m+1
 GO TO 50
 40 il(m) = k
 iu(m) = j
 j = l
 m = m+1
 GO TO 50
!                                  BEGIN AGAIN ON ANOTHER PORTION OF
!                                  THE UNSORTED ARRAY
 45 m = m-1
 IF (m == 0) RETURN
 i = il(m)
 j = iu(m)
 50 IF (j-i >= 11) GO TO 15
 IF (i == 1) GO TO 5
 i = i-1
 55 i = i+1
 IF (i == j) GO TO 45
 t = a(i+1)
 it = ir(i+1)
 IF (a(i) <= t) GO TO 55
 k = i
 60 a(k+1) = a(k)
 ir(k+1) = ir(k)
 k = k-1
 IF (t < a(k)) GO TO 60
 a(k+1) = t
 ir(k+1) = it
 GO TO 55
END SUBROUTINE hdvsr
