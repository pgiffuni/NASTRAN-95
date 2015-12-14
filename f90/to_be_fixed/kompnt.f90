FUNCTION kompnt (ig,ic,ideg,iw,icc,jg)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     THIS FUNCTION HAS AS ITS VALUE THE NUMBER OF COMPONENTS STORED
!     IN THE CONNECTION ARRAY IG.
!     ALSO, IC AND ICC ARE SET UP.
!     IC(I) =COMPONENT INDEX FOR NODE I
!     ICC(I)=THE STARTING POSITION TO BE USED FOR LABELS IN COMPONENT I
!     THUS, ICC(I+1)-ICC(I)= THE NUMBER OF NODES IN COMPONENT I
 
!     INTEGER          BUNPK
 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(OUT)                     :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(OUT)                     :: iw(1)
 INTEGER, INTENT(OUT)                     :: icc(1)
 INTEGER, INTENT(IN)                      :: jg(1)
 
 COMMON /bands /  nn,       mm
 
 DO  i=1,nn
   icc(i)=0
   ic(i)=0
 END DO
 nc=0
 icc(1)=1
 105 DO  i=1,nn
   IF (ic(i) == 0) THEN
     GO TO   120
   ELSE
     GO TO   110
   END IF
   kompnt=nc
 END DO
 GO TO 210
 120 nc=nc+1
 ki=0
 ko=1
 iw(1)=i
 ic(i)=nc
 IF (nc-1 < 0) THEN
   GO TO   130
 END IF
 125 is=icc(nc)+1
 icc(nc+1)=is
 130 ki=ki+1
 ii=iw(ki)
 n =ideg(ii)
 IF (n == 0) THEN
   GO TO   105
 END IF
 140 CALL bunpak(ig,ii,n,jg)
 DO  i=1,n
   ia=jg(i)
   IF (ic(ia) == 0) THEN
     GO TO   150
   ELSE
     GO TO   200
   END IF
   150 ic(ia)=nc
   ko=ko+1
   iw(ko)=ia
   is=icc(nc+1)+1
   icc(nc+1)=is
 END DO
 IF (ko-ki > 0) THEN
   GO TO   130
 ELSE
   GO TO   105
 END IF
 210 RETURN
END FUNCTION kompnt
