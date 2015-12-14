SUBROUTINE degree (ig,ideg,jg)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     SET UP THE IDEG ARRAY CONTAINING THE DEGREE OF EACH NODE STORED
!     IN THE IG ARRAY.
!     IDEG(I)=DEGREE OF NODE I
 
!     INTEGER          BUNPK
 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(OUT)                     :: ideg(1)
 INTEGER, INTENT(IN OUT)                  :: jg(1)
 
 COMMON /bands /  nn,       mm
 
 DO  i=1,nn
   ideg(i)=0
   CALL bunpak(ig,i,mm,jg)
   DO  j=1,mm
!     IF (BUNPK(IG,I,J)) 100,100,50
     IF (jg(j) > 0) THEN
       GO TO    50
     ELSE
       GO TO   100
     END IF
     50 ideg(i)=ideg(i)+1
   END DO
 END DO
 RETURN
END SUBROUTINE degree
