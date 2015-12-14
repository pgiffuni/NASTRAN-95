SUBROUTINE tldrs (offset,ii,trans,trans1)
     
!     &    ENTRY TLDRD (DFFSET,II,TRAND,TRAND1)
 
!     MODIFIED TO INCLUDE THE EFFECTS OF OFFSET
 
 
 REAL, INTENT(IN)                         :: offset
 INTEGER, INTENT(IN OUT)                  :: ii
 REAL, INTENT(IN)                         :: trans(1)
 REAL, INTENT(OUT)                        :: trans1(36)
 
 DOUBLE PRECISION :: trand(1),trand1(36),dffset
 
!     SINGLE PRECISION -
 
 DO  i = 1,36
   trans1(i) = 0.0
 END DO
 
 ipoint = 9*(ii-1)
 
 DO  i = 1,3
   jpoint = 6*(i-1)
   kpoint = jpoint  + 21
   lpoint = 3*(i-1) + ipoint
   
   DO  j = 1,3
     trans1(jpoint+j) = trans(lpoint+j)
     trans1(kpoint+j) = trans(lpoint+j)
   END DO
 END DO
 
 IF (offset == 0.0) GO TO 100
 trans1(4) = offset*trans(ipoint+4)
 trans1(5) = offset*trans(ipoint+5)
 trans1(6) = offset*trans(ipoint+6)
 trans1(10)=-offset*trans(ipoint+1)
 trans1(11)=-offset*trans(ipoint+2)
 trans1(12)=-offset*trans(ipoint+3)
 GO TO 100
 
 ENTRY tldrd (dffset,ii,trand,trand1)
!     ====================================
 
!     DOUBLE PRECISION -
 
 DO  i = 1,36
   trand1(i) = 0.0D0
 END DO
 
 ipoint = 9*(ii-1)
 
 DO  i = 1,3
   jpoint = 6*(i-1)
   kpoint = jpoint  + 21
   lpoint = 3*(i-1) + ipoint
   
   DO  j = 1,3
     trand1(jpoint+j) = trand(lpoint+j)
     trand1(kpoint+j) = trand(lpoint+j)
   END DO
 END DO
 
 IF (dffset == 0.0D0) GO TO 100
 trand1(4) = dffset*trand(ipoint+4)
 trand1(5) = dffset*trand(ipoint+5)
 trand1(6) = dffset*trand(ipoint+6)
 trand1(10)=-dffset*trand(ipoint+1)
 trand1(11)=-dffset*trand(ipoint+2)
 trand1(12)=-dffset*trand(ipoint+3)
 100 RETURN
END SUBROUTINE tldrs
