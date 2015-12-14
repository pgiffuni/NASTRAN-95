SUBROUTINE dgree (ndstk,ndeg,iold,ibw1,ipf1,nu)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
!     DGREE COMPUTES THE DEGREE OF EACH NODE IN NDSTK AND STORES
!     IT IN THE ARRAY NDEG.  THE BANDWIDTH AND PROFILE FOR THE ORIGINAL
!     OR INPUT RENUMBERING OF THE GRAPH IS COMPUTED ALSO.
 
!     COMPUTE MAXIMUM DEGREE MM AND STORE IN IDEG.
 
!     INTEGER          BUNPK
 
 INTEGER, INTENT(IN OUT)                  :: ndstk(1)
 INTEGER, INTENT(OUT)                     :: ndeg(1)
 INTEGER, INTENT(IN)                      :: iold(1)
 INTEGER, INTENT(OUT)                     :: ibw1
 INTEGER, INTENT(OUT)                     :: ipf1
 INTEGER, INTENT(IN)                      :: nu(1)
 
 COMMON /bandg /  n,        idpth,    ideg
 COMMON /bands /  nn,       mm
 
 ibw1=0
 ipf1=0
 ideg=mm
 mm=0
 DO  i=1,n
   ndeg(i)=0
   irw=0
   CALL bunpak(ndstk,i,ideg,nu)
   DO  j=1,ideg
     itst=nu(j)
     IF (itst > 0) THEN
       GO TO    50
     ELSE
       GO TO    90
     END IF
     50 ndeg(i)=ndeg(i)+1
     idif=iold(i)-iold(itst)
     IF (irw < idif) irw=idif
     mm=MAX0(mm,j)
   END DO
   90 ipf1=ipf1+irw
   IF (irw > ibw1) ibw1=irw
 END DO
 ideg=mm
 
!     INCLUDE DIAGONAL TERMS IN BANDWIDTH AND PROFILE
 ibw1=ibw1+1
 ipf1=ipf1+n
 RETURN
END SUBROUTINE dgree
