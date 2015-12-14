SUBROUTINE na1 2 IF (*,a,n,b,INT)
     
!     VAX, IBM AND UNIVAC VERSION (CHARACTER FUNCTION PROCESSING)
!     ===========================
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: a(1)
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(OUT)                        :: b
 INTEGER, INTENT(IN OUT)                  :: INT
 COMMON /xreadx/nout
 
 CHARACTER (LEN=1) :: bk,    pt,    tj,   t(24), c(1), num(10)
 CHARACTER (LEN=12) :: temp,  next,  blnk
 EQUIVALENCE    (temp,t(1)), (next,t(13)), (i,xi)
 DATA           bk,    pt,    blnk  / ' ', '.', '            ' /
 DATA           num / '0','1','2','3','4','5','6','7','8','9'  /
 
!     ARRAY A, IN NA1 BCD WORDS (OR C IN CHARACTERS), IS DECODED TO
!     AN INTEGER OR TO A F.P. NUMBER IN B.
!     INT SHOULD BE SET TO +1 IF CALLER IS EXPECTING B TO BE AN INTEGER,
!     OR SET TO -1 IF B IS TO BE A F.P. NUMBER.   SET INT TO ZERO IF
!     CALLER IS NOT SURE.  IN THIS LAST CASE, INT WILL BE SET TO +1 OR
!     -1 BY NA12IF/NK12IF ACCORDING TO THE INPUT DATA TYPE.
!     THESE ROUTINES HANDLE UP TO 12 DIGITS INPUT DATA (N .LE. 12)
!     (NO SYSTEM ENCODE/DECODE FUNCTIONS ARE USED)
 
!     ENTRY POINTS   NA1 2 IF  (BCD-INTEGER/FP VERSION)
!                    NK1 2 IF  (CHARACTER-INTEGER/FP VERSION)
 
!     WRITTEN BY G.CHAN/SPERRY IN AUG. 1985
!     PARTICULARLY FOR XREAD ROUTINE, IN SUPPORT OF ITS NEW FREE-FIELD
!     INPUT FORMAT.  THIS SUBROUTINE IS MACHINE INDEPENDENT
 
 IF (n > 12) GO TO 150
 CALL b2k (a,temp,n)
 GO TO 20
 
 ENTRY nk1 2 IF (*,c,n,b,INT)
!     ****************************
 
 IF (n > 12) GO TO 150
 DO  i=1,n
   t(i)=c(i)
 END DO
 
 20   IF (INT >= 1) GO TO 110
 25   nt=1
 k =24
 j =n
 next=blnk
 DO  i=1,12
   IF (i    >  n) GO TO 30
   tj=t(j)
   IF (tj == bk) GO TO 50
   IF (tj == pt) nt=nt-2
   t(k)=tj
   GO TO 40
   30   t(k)=bk
   40   k=k-1
   50   j=j-1
 END DO
 
 IF (nt < -1 .OR. INT*nt < 0) GO TO 170
 IF (INT == 0) INT=nt
 IF (INT) 60,170,80
 60   READ (next,70) b
 70   FORMAT (f12.0)
 RETURN
 80   READ (next,90) i
 90   FORMAT (i12)
 100  b=xi
 RETURN
 
!     QUICK WAY TO GET THE INTEGER
 
 110  i=0
 j=0
 120  j=j+1
 IF (j > n) GO TO 100
 tj=t(j)
 IF (tj == bk) GO TO 120
 DO  k=1,10
   IF (tj == num(k)) GO TO 140
 END DO
 GO TO 25
 140  i=i*10 + k-1
 GO TO 120
 
 150  b=0.
 WRITE (nout,160) n
 160  FORMAT (5X,'*** N.GT.12/NA12IF',i6)
 170  RETURN 1
END SUBROUTINE na1 2 IF
