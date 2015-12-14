LOGICAL FUNCTION ifpdco (ic)
     
!     DECODE D.O.F. INTO LL SPACE.
!     RETURN WITH IFPDCO=.TRUE. IF ERROR ENCOUNTERED
!     FOR EXAMPLE - GIVEN IC=124, THEN
!                   LL(1)=1,   LL(2)=2,  LL(4)=4, LL(3)=LL(5)=LL(6)=0
!                   GC(1)=124, GC(2)=12, GC(3)=1, GC(4)=GC(5),GC(6)=0
!                   IFPDCO=.FALSE.
 
 
 INTEGER, INTENT(IN)                      :: ic
 INTEGER :: dg,gc
 COMMON /ifpdta/ dum(521),gc(7),ll(6)
 COMMON /system/ idummy(55),ithrml
 
 gc(1) = ic
 DO  lc=1,6
   ll(lc) = 0
 END DO
 IF (ic < 0) THEN
   GO TO   120
 ELSE IF (ic == 0) THEN
   GO TO   116
 END IF
 112 DO  lc=1,6
   gc(lc+1) = gc(lc)/10
   dg = gc(lc)-10*gc(lc+1)
   IF (ithrml /= 1 .AND. dg > 6) GO TO 120
   IF (ithrml == 1 .AND. dg > 1) GO TO 120
   IF (dg == 0) GO TO 118
   IF (ll(dg) /= 0) GO TO 120
   ll(dg) = dg
 END DO
 IF (gc(7) /= 0) GO TO 120
 116 ifpdco = .false.
 RETURN
 118 IF (gc(lc) == 0) GO TO 116
 120 ifpdco = .true.
 RETURN
END FUNCTION ifpdco
