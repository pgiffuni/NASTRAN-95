SUBROUTINE diam (nc,maxdeg,nl,nodesl,idem,maxlev,ig,ic,ideg,  &
        idis,iw,icc,jg)
     
!     DETERMINE NL STARTING POINTS AND STORE IN NODESL.
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
 
 INTEGER, INTENT(IN OUT)                  :: nc
 INTEGER, INTENT(IN OUT)                  :: maxdeg
 INTEGER, INTENT(OUT)                     :: nl
 INTEGER, INTENT(OUT)                     :: nodesl(1)
 INTEGER, INTENT(IN OUT)                  :: idem
 INTEGER, INTENT(OUT)                     :: maxlev
 INTEGER, INTENT(IN)                      :: ig(1)
 INTEGER, INTENT(IN)                      :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(IN)                      :: idis(1)
 INTEGER, INTENT(IN)                      :: iw(1)
 INTEGER, INTENT(IN)                      :: icc(1)
 INTEGER, INTENT(IN OUT)                  :: jg(1)
 
 COMMON /bands / nn
 
 nl = 0
 maxlev = 600000
 DO  i = 1,nn
   IF (nc-ic(i) == 0) THEN
     GO TO    40
   ELSE
     GO TO   100
   END IF
   40 IF (maxdeg-ideg(i) < 0) THEN
     GO TO   100
   END IF
   50 md = idist(i,ml,maxlev,ig,ic,ideg,idis,iw,icc,jg)
   IF (md > 0) THEN
     GO TO    60
   ELSE
     GO TO   120
   END IF
   60 IF (ml-maxlev < 0) THEN
     GO TO    70
   ELSE IF (ml-maxlev == 0) THEN
     GO TO    80
   ELSE
     GO TO   100
   END IF
   70 maxlev = ml
   nl = 1
   nodesl(1) = i
   CYCLE
   80 IF (nl >= idem) CYCLE
   nl = nl + 1
   nodesl(nl) = i
 END DO
 RETURN
 
 120 ml = 1
 nodesl(1) = i
 maxlev = 0
 RETURN
END SUBROUTINE diam
