FUNCTION maxdgr (nc,ic,ideg)
     
 INTEGER, INTENT(IN OUT)                  :: nc
 INTEGER, INTENT(IN OUT)                  :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 
 COMMON /bands / nn
 
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DEGREE OF ANY NODE OF
!     COMPONENT NC IF NC.GT.0
!     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
 
 m=0
 DO  i=1,nn
   IF (nc == 0) THEN
     GO TO    50
   END IF
   40 IF (ic(i) -nc == 0) THEN
     GO TO    50
   ELSE
     GO TO   100
   END IF
   50 IF (ideg(i)-m > 0) THEN
     GO TO    60
   ELSE
     GO TO   100
   END IF
   60 m=ideg(i)
 END DO
 maxdgr=m
 RETURN
END FUNCTION maxdgr
