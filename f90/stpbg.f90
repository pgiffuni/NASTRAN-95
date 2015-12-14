SUBROUTINE stpbg(bm,gm,ns,bloc,d,ca,nsize)
!     MAKES MATRICES BM AND GM FOR EACH STRIP
 
 
 REAL, INTENT(OUT)                        :: bm(4,4,ns)
 REAL, INTENT(OUT)                        :: gm(4,3,ns)
 INTEGER, INTENT(IN)                      :: ns
 REAL, INTENT(IN)                         :: bloc(1)
 REAL, INTENT(IN)                         :: d(1)
 REAL, INTENT(IN)                         :: ca(1)
 INTEGER, INTENT(IN OUT)                  :: nsize(1)
 
 DO  n=1,ns
   DO  i=1,4
     DO  j=1,4
       bm(i,j,n)=0.0
     END DO
   END DO
   DO  i=1,4
     DO  j=1,3
       gm(i,j,n)=0.0
     END DO
   END DO
   bm(1,1,n)=  bloc(n)
   bm(2,2,n)= -bloc(n)*bloc(n)
   gm(1,1,n)=-1.0/bloc(n)
   gm(2,2,n)= 1.0
   IF(nsize(n) == 2) GO TO 50
!         CONTROL SURFACE CASE
   e= ca(n) + d(n) - 1.5*bloc(n)
   bm(3,3,n)=e*bloc(n)
   bm(3,4,n)= bm(2,2,n)
   gm(3,3,n)= 1.0
   gm(4,3,n)= -e/bloc(n)
   50 CONTINUE
 END DO
 RETURN
END SUBROUTINE stpbg
