SUBROUTINE scat (kg,ncon,inv,ii3,norig)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS ROUTINE USES SCATTER SORT TECHNIQUES FOR EACH GRID POINT
!     ENCOUNTERED TO DETERMINE WHETHER OR NOT THE POINT HAS BEEN SEEN
!     BEFORE.   IF NOT, INV, NORIG, AND NN ARE UPDATED.
 
!     INV(I,1) CONTAINS AN ORIGINAL GRID POINT NUMBER
!     INV(I,2) CONTAINS THE INTERNAL NUMBER ASSIGNED TO IT (BEFORE SORT)
 
 
 INTEGER, INTENT(IN OUT)                  :: kg(1)
 INTEGER, INTENT(IN)                      :: ncon
 INTEGER, INTENT(OUT)                     :: inv(ii3,2)
 INTEGER, INTENT(IN OUT)                  :: ii3
 INTEGER, INTENT(OUT)                     :: norig(1)
 
 COMMON /bandb / dum3b(3), ngrid
 COMMON /bands / nn,       dum3(3),  maxgrd,   maxdeg,   kmod
 COMMON /system/ isys,     nout
 
 IF (ncon < 1) RETURN
 DO  i = 1,ncon
   nold = kg(i)
   IF (nold == 0) CYCLE
   loc = nold - 1
   20 loc = MOD(loc,kmod) + 1
   IF (inv(loc,1) /= 0) GO TO 30
   inv(loc,1) = nold
   nn = nn + 1
   IF (nn > maxgrd) GO TO 60
   norig(nn) = nold
   inv(loc,2) = nn
   GO TO 40
   30 IF (inv(loc,1) /= nold) GO TO 20
   40 kg(i) = inv(loc,2)
 END DO
 RETURN
 
 60 ngrid = -1
 RETURN
END SUBROUTINE scat
