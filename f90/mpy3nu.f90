SUBROUTINE mpy3nu (iz)
     
!     CALCULATES NEXT TIME USED FOR INDIVIDUAL COLUMNS OF B OR FOR ROWS
!     CORRESPONDING TO NON-ZERO TERMS IN COLUMN OF A.
 
 
 INTEGER, INTENT(IN)                      :: iz(1)
 INTEGER :: zpntrs
 DIMENSION  NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,nout
 COMMON /mpy3cp/ itrl,icore,n,ncb,dum1(4),zpntrs(22),laend, dum2(8),j,id,ntbu
 EQUIVALENCE     (ipoint,zpntrs(3)),(iacols,zpntrs(5))
 DATA    NAME  / 4HMPY3,4HNU   /
 
!     CALCULATION BY SEARCH THROUGH ROW OF A IN QUESTION.
 
 lp = ipoint + id - 1
 l1 = iz(lp)
 IF (l1 ==   0) GO TO 60
 IF (id == ncb) GO TO 20
 ll = id + 1
 DO  l = ll,ncb
   lp = lp + 1
   IF (iz(lp) == 0) CYCLE
   l2 = iz(lp) - 1
   GO TO 30
 END DO
 20 l2  = laend
 30 lac = iacols + l1 - 2
 DO  l = l1,l2
   lac = lac + 1
   IF (j < iz(lac)) GO TO 50
 END DO
 ntbu = 99999999
 GO TO 80
 50 ntbu = iz(lac)
 GO TO 80
 
!    ERROR MESSAGE.
 
 60 WRITE  (nout,70) ufm
 70 FORMAT (a23,' 6557, UNEXPECTED NULL COLUMN OF A(T) ENCOUNTERED.')
 CALL mesage (-37,0,NAME)
 
 80 RETURN
END SUBROUTINE mpy3nu
