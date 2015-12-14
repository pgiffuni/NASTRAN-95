SUBROUTINE fvrs1d (base,base1,INDEX,nfx)
     
 
 COMPLEX, INTENT(IN OUT)                  :: base(3,nfx)
 COMPLEX, INTENT(OUT)                     :: base1(3,nfx)
 INTEGER, INTENT(IN)                      :: INDEX(nfx)
 INTEGER, INTENT(IN)                      :: nfx
 
 
 
 
 DO  i=1,nfx
   loc =INDEX(i)
   DO  l=1,3
     base1(l,i)=base(l,loc)
   END DO
 END DO
 
!-----RETURN BASE1 TO BASE
 
 DO  i=1,nfx
   DO  l=1,3
     base(l,i)=base1(l,i)
   END DO
 END DO
 RETURN
END SUBROUTINE fvrs1d
