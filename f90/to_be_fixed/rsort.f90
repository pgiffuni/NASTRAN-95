SUBROUTINE rsort (nwds,keywx,l,nx)
     
!     RSORT SORTS REAL NUMBERS IN L(NWDS,NCOL)
!           WHERE NCOL = IABS(NX)/NWDS
 
!     IF KEYWX .LT. 0 SORT BY ABSOLUTE VALUE
!     IF NX    .LT. 0 SORT IN DECREASING SEQUENCE
 
!     COMMENTS FROM G.C./UNISYS
!     THIS ROUTINE IS INEFFICIENT FOR LARGE ARRAY OF L
 
 
 INTEGER, INTENT(IN)                      :: nwds
 INTEGER, INTENT(IN OUT)                  :: keywx
 REAL, INTENT(IN OUT)                     :: l(1)
 INTEGER, INTENT(IN OUT)                  :: nx
 LOGICAL :: mag,bck
 REAL :: temp(50)
 INTEGER :: nam(2)
 COMMON /system/ ibuf,nout
 DATA    nam   / 4HRSOR, 4HT    /
 
 IF (nwds <= 50) GO TO 30
 WRITE  (nout,20)
 20 FORMAT (' *** ARRAY TEMP OF 50 EXCEEDED')
 CALL mesage (-37,0,nam)
 
 30 mag = .false.
 bck = .false.
 IF (keywx < 0) mag = .true.
 IF (nx    < 0) bck = .true.
 keywd= IABS(keywx)
 nnn  = IABS(nx)
 iii  = nwds+keywd
 ia   = nwds-keywd
 IF (nnn-nwds-nwds < 0) GO TO 150
 DO  i = iii,nnn,nwds
   jj = i-nwds
   IF (bck) GO TO 40
   IF (mag) IF (ABS(l(i))-ABS(l(jj))) 50,140,140
   IF (l(i)-l(jj) < 0) THEN
     GO TO    50
   ELSE
     GO TO   140
   END IF
   40 IF (mag) IF (ABS(l(jj))-ABS(l(i))) 50,140,140
   IF (l(jj)-l(i) < 0) THEN
     GO TO    50
   ELSE
     GO TO   140
   END IF
   50 jj = jj-nwds
   IF (jj <= 0) GO TO 70
   IF (bck) GO TO 60
   IF (mag) IF (ABS(l(i))-ABS(l(jj))) 50,80,80
   IF (l(i) - l(jj) < 0) THEN
     GO TO    50
   ELSE
     GO TO    80
   END IF
   60 IF (mag) IF (ABS(l(jj))-ABS(l(i))) 50,80,80
   IF (l(jj)-l(i) < 0) THEN
     GO TO    50
   ELSE
     GO TO    80
   END IF
   70 jj = nwds
   GO TO 90
   80 jj = jj+ia+nwds
   90 ii = i-keywd
   DO  j = 1,nwds
     ii = ii+1
     temp(j) = l(ii)
   END DO
   110 iia = ii-nwds
   l(ii) = l(iia)
   ii = ii-1
   IF (ii-jj > 0) THEN
     GO TO   110
   END IF
   120 ii = ii-nwds
   DO  j = 1,nwds
     ii = ii+1
     l(ii) = temp(j)
   END DO
 END DO
 
 150 RETURN
END SUBROUTINE rsort
