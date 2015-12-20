FUNCTION mindeg (nc,ic,ideg)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS FUNCTION HAS AS ITS VALUE THE MINIMUM DEGREE OF ANY NODE OF
!     COMPONENT NC IF NC.GT.0
!     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
 
 
 INTEGER, INTENT(IN OUT)                  :: nc
 INTEGER, INTENT(IN OUT)                  :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 
 COMMON /bands / nn
 
 m=600000
 DO  i=1,nn
   IF (nc == 0) THEN
     GO TO    50
   END IF
   40 IF (ic(i) -nc == 0) THEN
     GO TO    50
   ELSE
     GO TO   100
   END IF
   50 IF (m-ideg(i) > 0) THEN
     GO TO    60
   ELSE
     GO TO   100
   END IF
   60 m=ideg(i)
   100 CONTINUE
 END DO
 mindeg=m
 RETURN
END FUNCTION mindeg
