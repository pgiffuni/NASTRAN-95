SUBROUTINE bug (NAME,loc,buf,nwds)
     
!     THIS ROUTINE PRINTS NAME,LOC, AND CONTENT OF BUF ARRAY
!     E.G.   CALL BUG ('SUBR ABC',105,CORE(1),120)
!     LIMITED TO 5000 LINES EACH CALL,  14 VALUES PER LINE
 
!     (THIS ROUTINE REPLACES THE OLD ONE IN NASTRAN)
!     WRITTEN BY G.CHAN/SPERRY     MARCH 1986
 
 
 REAL, INTENT(IN OUT)                     :: NAME(3)
 INTEGER, INTENT(IN OUT)                  :: loc
 REAL, INTENT(IN OUT)                     :: buf(1)
 INTEGER, INTENT(IN OUT)                  :: nwds
 
 CHARACTER (LEN=4) :: a(28),    xloc,     BLANK
 CHARACTER (LEN=8) :: b(14),    zero,     ERR
 COMMON /system/ ibuf,     nout
 EQUIVALENCE     (a(1),b(1))
 DATA    line,   nwpl,     limit              /  &
     0,      14,       5000               /
 DATA    zero,   BLANK,    xloc,     ERR      /  &
     ' 00 ', '    ',   'LOC',    '(ERR)'  /
 
 CALL sswtch (20,l)
 IF (l == 0) RETURN
 GO TO 5
 
 ENTRY bug1 (NAME,loc,buf,nwds)
!     ==============================
 
 5    IF (nwds < 0) RETURN
 l = 2
 i = 0
 CALL a42k8 (NAME(1),NAME(2),b(1))
 CALL int2k8 (*20,loc,a(3))
 a(4) = a(3)
 a(3) = xloc
 
 10   IF (i >= nwds) GO TO 60
 15   i = i + 1
 l = l + 1
 j = numtyp(buf(i)) + 1
 SELECT CASE ( j )
   CASE (    1)
     GO TO   25
   CASE (    2)
     GO TO  30
   CASE (    3)
     GO TO   35
   CASE (    4)
     GO TO  40
 END SELECT
!            ZERO,INT,REAL,BCD
 20   b(l) = ERR
 GO TO 55
 25   b(l) = zero
 GO TO 55
 30   CALL int2k8 (*20,buf(i),b(l))
 GO TO 55
 35   CALL fp2k8  (*20,buf(i),b(l))
 GO TO 55
 40   CALL a42k8 (buf(i),buf(i+1),b(l))
 IF (numtyp(buf(i+1)) /= 3) GO TO 45
 i = i + 1
 GO TO 50
 45   a(l*2) = BLANK
 50   IF (i >= nwds) GO TO 60
 55   IF (l < nwpl) GO TO 10
 60   IF (l > 0) WRITE (nout,65) (b(j),j=1,l)
 65   FORMAT (2X,14(a8,1X))
 line = line + 1
 IF (line > limit) GO TO 70
 l = 0
 IF (i < nwds) GO TO 15
 RETURN
 
 70   WRITE  (nout,75) limit
 75   FORMAT (/2X,'PRINT LINES IN BUG EXCEEDS LIMIT OF',i6)
 RETURN
END SUBROUTINE bug
