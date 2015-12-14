SUBROUTINE stack (ideg,NEW,ild,iw)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     STACK POINTS OF ZERO DEGREE AT END OF THE NUMBERING.
!     IW IS SCRATCH STORAGE.
 
 
 INTEGER, INTENT(IN OUT)                  :: ideg(1)
 INTEGER, INTENT(IN OUT)                  :: NEW(1)
 INTEGER, INTENT(IN OUT)                  :: ild(1)
 INTEGER, INTENT(OUT)                     :: iw(1)
 
 COMMON /bands / nn
 COMMON /bandd / dum(5),   kt
 
 kt  = 0
 nn1 = nn - 1
 
!     LIST POINTS OF ZERO DEGREE AND INCREMENT COUNTER KT.
 
 DO  i = 1,nn
   IF (ideg(i) > 0) CYCLE
   kt = kt + 1
   iw(kt) = ild(i)
 END DO
 IF (kt <= 0) GO TO 80
 
!     SORT LIST OF RENUMBERED NUMBERS TO BE STACKED.
 
 IF (kt <= 1) GO TO 40
 kt1 = kt-1
 DO  i = 1,kt1
   k = kt - i
   kflag = 0
   DO  j = 1,k
     j1 = j + 1
     IF (iw(j) <= iw(j1))  CYCLE
     kflag = 1
     l = iw(j)
     iw(j ) = iw(j1)
     iw(j1) = l
   END DO
   IF (kflag == 0) EXIT
 END DO
 40 CONTINUE
 
!     STACK POINTS OF ZERO DEGREE AT END OF NEW.
 
 DO  l = 1,kt
   i = iw(l) - l + 1
   k = NEW(i)
   IF (i >= nn) GO TO 60
   DO  j = i,nn1
     NEW(j) = NEW(j+1)
   END DO
   60 NEW(nn) = k
 END DO
 
!     CORRECT ILD, THE INVERSE OF NEW.
 
 80 DO  i = 1,nn
   k = NEW(i)
   ild(k) = i
 END DO
 RETURN
END SUBROUTINE stack
