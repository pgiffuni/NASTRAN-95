SUBROUTINE splt10( icomp, comps, nc)
     
 
 INTEGER, INTENT(IN OUT)                  :: icomp
 INTEGER, INTENT(OUT)                     :: comps(9)
 INTEGER, INTENT(OUT)                     :: nc
 
 IF( icomp == 0) icomp= 1
 ic = icomp
 nc = 0
 DO   i=1,9
   ix = ic/10
   jx = ic - 10*ix
   ic = ix
   IF ( jx == 0) GO TO 5
   nc =nc+1
   comps(nc)= jx
   5 IF( ic == 0) EXIT
 END DO
 15 IF (nc == 1) RETURN
 CALL sort (0,0,1,1, comps,nc )
 
!     REMOVE DUPLICATES
 ix= 1
 DO  i=2,nc
   IF ( comps(i) == comps(i-1)) CYCLE
   ix = ix+1
   comps(ix) = comps(i)
 END DO
 nc = ix
 RETURN
END SUBROUTINE splt10
